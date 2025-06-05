function dmrg(
        PH::ProjMPO,
        psi0::MPS,
        sweeps::Sweeps;
        degeneracy_tol=1e-10,
        maxsweeps=10,
        svd_alg=nothing,
        which_decomp="svd",
        observer=DMRGObserver(;energy_tol=1e-10, minsweeps=10),
        outputlevel=1,
        # eigsolve kwargs
        eigsolve_tol=1e-14,
        eigsolve_krylovdim=3,
        eigsolve_maxiter=1,
        eigsolve_verbosity=0,
        eigsolve_which_eigenvalue=:SR,
        ishermitian=true,
    )

    if length(psi0) == 1
        error(
        "`dmrg` currently does not support system sizes of 1. You can diagonalize the MPO tensor directly with tools like `LinearAlgebra.eigen`, `KrylovKit.eigsolve`, etc.",
        )
    end

    psi = copy(psi0)
    N = length(psi)
    if !isortho(psi) || orthocenter(psi) != 1
        psi = orthogonalize!(PH, psi, 1)
    end
    @assert isortho(psi) && orthocenter(psi) == 1

    PH = position!(PH, psi, 1)
    energy = 0.0

    for sw in 1:nsweep(sweeps)
        for ni in 1:maxsweeps
            ni_time = @elapsed begin
            maxtruncerr = 0.0

            for (b, ha) in sweepnext(N)
                @timeit_debug timer "dmrg: position!" begin
                PH = position!(PH, psi, b)
                end

                @timeit_debug timer "dmrg: psi[b]*psi[b+1]" begin
                phi = psi[b] * psi[b + 1]
                end

                @timeit_debug timer "dmrg: eigsolve" begin
                vals, vecs = eigsolve(
                    PH,
                    phi,
                    1,
                    eigsolve_which_eigenvalue;
                    ishermitian,
                    tol=eigsolve_tol,
                    krylovdim=eigsolve_krylovdim,
                    maxiter=eigsolve_maxiter,
                    verbosity=eigsolve_verbosity,
                )
                end


                energy = first(vals)
                ## Right now there is a conversion problem in CUDA.jl where `UnifiedMemory` Arrays are being converted
                ## into `DeviceMemory`. This conversion line is here temporarily to fix that problem when it arises
                ## Adapt is only called when using CUDA backend. CPU will work as implemented previously.
                ## TODO this might be the only place we really need iscu if its not fixed.
                phi = if NDTensors.iscu(phi) && NDTensors.iscu(first(vecs))
                adapt(ITensors.set_eltype(unwrap_array_type(phi), eltype(first(vecs))), first(vecs))
                else
                first(vecs)
                end

                ortho = ha == 1 ? "left" : "right"

                drho = nothing
                if noise(sweeps, sw) > 0
                @timeit_debug timer "dmrg: noiseterm" begin
                    # Use noise term when determining new MPS basis.
                    # This is used to preserve the element type of the MPS.
                    elt = real(scalartype(psi))
                    drho = elt(noise(sweeps, sw)) * noiseterm(PH, phi, ortho)
                end
                end

                @timeit_debug timer "dmrg: replacebond!" begin

                # TODO: If structure will be optimized, we need to take check the degeneracy on this structure
                maxdim_degeneracy = maxdim_with_degeneracycheck(psi, b, phi, maxdim(sweeps, sw); degeneracy_tol=degeneracy_tol)

                spec = replacebond!(
                    PH,
                    psi,
                    b,
                    phi;
                    maxdim=maxdim_degeneracy,
                    mindim=mindim(sweeps, sw),
                    cutoff=cutoff(sweeps, sw),
                    eigen_perturbation=drho,
                    ortho,
                    normalize=true,
                    which_decomp,
                    svd_alg,
                )
                end

                maxtruncerr = max(maxtruncerr, spec.truncerr)

                if outputlevel >= 2
                @printf("Sweep %d, half %d, bond (%d,%d) energy=%s\n", sw, ha, b, b + 1, energy)
                @printf(
                    "  Truncated using cutoff=%.1E maxdim=%d mindim=%d\n",
                    cutoff(sweeps, sw),
                    maxdim(sweeps, sw),
                    mindim(sweeps, sw)
                )
                @printf(
                    "  Trunc. err=%.2E, bond dimension %d\n", spec.truncerr, dim(linkind(psi, b))
                )
                flush(stdout)
                end

                sweep_is_done = (b == 1 && ha == 2)
                measure!(
                    observer;
                    energy,
                    psi,
                    projected_operator=PH,
                    bond=b,
                    sweep=sw,
                    half_sweep=ha,
                    spec,
                    outputlevel,
                    sweep_is_done,
                )
                # TODO: I don't know much about the measure! function, but it seems to be used to check for convergence
            end
            end
            # Sweep is done
            if outputlevel >= 1
                @printf(
                    "At the sweeps stage %d, after iteration %d, energy=%s  maxlinkdim=%d maxerr=%.2E time=%.3f\n",
                    sw,
                    ni,
                    energy,
                    maxlinkdim(psi),
                    maxtruncerr,
                    ni_time
                )
                flush(stdout)
            end
            isdone = checkdone!(observer; energy, psi, sweep=ni, outputlevel)
            isdone && break
        end
    end
    return (energy, psi)
end


function checkdone!(
        o::DMRGObserver; outputlevel=false, energy=nothing, psi=nothing, sweep=nothing
    )
    sweep == 1 && return false
    if length(real(energies(o))) > o.minsweeps
        outputlevel > 0 && println("Minimum number of sweeps reached, stopping sweeps stage")
        return true
    end
    if (
        abs(real(energies(o))[end]/real(energies(o))[end - 1] - 1.0) < o.etol
        # TODO: structure of the MPS can be used to check for convergence
    )
        outputlevel > 0 && println("Relative energy difference less than $(o.etol), stopping sweeps stage")
        return true
    end
    return false
end


function maxdim_with_degeneracycheck(M::MPS, b::Int, phi::ITensor, maxdim::Int; degeneracy_tol=1e-10)
    U, S, V, spec = svd(phi, inds(M[b]), cutoff=0.0)
    length(spec.eigs) <= maxdim && return maxdim
    dim = maxdim
    while dim > 1
        if abs(spec.eigs[dim+1] - spec.eigs[dim]) < degeneracy_tol
            dim -= 1
            dim == 1 && error("Degeneracy found at the lowest dimension, this is not allowed")
        else
            break
        end
    end
    return dim
end