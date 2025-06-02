using ITensors, ITensorMPS
using FlexibleDMRG: dmrg

function main()
    let
        # Create 100 spin-one indices
        N = 10
        sites = siteinds("S=1/2", N)

        # Heisenberg model
        os = OpSum()
        for j in 1:(N - 1)
            os += "Sz", j, "Sz", j + 1
            os += 0.5, "S+", j, "S-", j + 1
            os += 0.5, "S-", j, "S+", j + 1
        end
        H = MPO(os, sites)

        # Create an initial random matrix product state
        psi0 = random_mps(sites)

        # Plan to do 5 passes or 'sweeps' of DMRG,
        # setting maximum MPS internal dimensions
        # for each sweep and maximum truncation cutoff
        # used when adapting internal dimensions:
        sweeps = Sweeps([
            "maxdim" "cutoff"
            10 1E-10
            20 1E-10
            30 1E-10
            40 1E-10
        ])
        kwargs = (
            degeneracy_tol = 1e-10,
            maxiters = 10,
            svd_alg=nothing,
            which_decomp="svd",
            observer = DMRGObserver(;energy_tol=1e-10, minsweeps=10),
            outputlevel = 1,
            # eigsolve kwargs
            eigsolve_krylovdim = 20,
            eigsolve_maxiter = 1,
        )

        # Run the DMRG algorithm, returning energy
        # (dominant eigenvalue) and optimized MPS
        energy, psi = dmrg(H, psi0, sweeps; kwargs...)
        println("Final energy = $energy")
        nothing
    end
end

main()