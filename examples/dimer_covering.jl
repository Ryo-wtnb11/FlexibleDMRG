using Random
Random.seed!(123)
using ITensors, ITensorMPS
using FlexibleDMRG: siteordering, MinimumEESiteOrderProposer

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

        H = MPO(os, sites)

        # Create an initial random matrix product state
        psi0 = random_mps(sites)

        # Plan to do 1 steps of sweeps of DMRG. Note that it is not the same as the number of sweeps as in the original dmrg function of ITensorMPS
        sweeps = Sweeps([
            "maxdim" "cutoff" "noise"
            100 1E-12 1e-9
        ])
        # DMRG kwargs
        dmrg_kwargs = (
            degeneracy_tol = 1e-8,
            outputlevel = 0,
            # eigsolve kwargs
            eigsolve_krylovdim = 25,
            eigsolve_maxiter = 1
        )

        kwargs = (
            nswapsites = 3,
            proposer = MinimumEESiteOrderProposer(),
        )

        # Run the DMRG algorithm, returning energy
        # (dominant eigenvalue) and optimized MPS
        energy, psi, whole_order = siteordering(H, psi0, sweeps; kwargs..., dmrg_kwargs=dmrg_kwargs)
        @show pairs
        @show whole_order
        return 0
    end
end

main()