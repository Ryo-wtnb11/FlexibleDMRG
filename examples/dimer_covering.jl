using Random
Random.seed!(123)
using ITensors, ITensorMPS
using FlexibleDMRG: dmrg, siteordering

function unique_random_pairs(N::Int)
    N % 2 == 0 || error("N must be even here.")
    sites = collect(1:N)
    shuffle!(sites)

    pairs = Vector{Tuple{Int, Int}}(undef, N÷2)
    for i in 1:N÷2
        pairs[i] = (sites[2i-1], sites[2i])
    end

    return pairs
end

function main()
    let
        # Create 100 spin-one indices
        N = 10
        sites = siteinds("S=1/2", N)

        # Heisenberg model
        os = OpSum()
        pairs = unique_random_pairs(N)
        for (i, j) in pairs
            os += "Sz", i, "Sz", j
            os += 0.5, "S+", i, "S-", j
            os += 0.5, "S-", i, "S+", j
        end

        @show pairs

        H = MPO(os, sites)

        # Create an initial random matrix product state
        psi0 = random_mps(sites)

        # Plan to do 4 steps of sweeps of DMRG. Note that it is not the same as the number of sweeps as in the original dmrg function of ITensorMPS
        sweeps = Sweeps([
            "maxdim" "cutoff"
            100 1E-10
        ])
        kwargs = (
            degeneracy_tol = 1e-10,
            outputlevel = 1,
            # eigsolve kwargs
            eigsolve_krylovdim = 20,
            eigsolve_maxiter = 2,
        )

        # Run the DMRG algorithm, returning energy
        # (dominant eigenvalue) and optimized MPS
        energy, psi = dmrg(H, psi0, sweeps; kwargs...)

        psi, H, whole_order = siteordering(H, psi, sweeps; kwargs...)
        @show whole_order

        return 0
    end
end

main()