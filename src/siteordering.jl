using Combinatorics

function siteordering(
    H::MPO,
    psi::MPS,
    sweeps::Sweeps;
    sweeps_of_siteordering=nothing, # Sweeps for DMRG in site ordering. If nothing, then it will be set to the last sweep of the sweeps.
    nswapsites::Int=3, # that should be around 3 - 4.
    maxswapsweeps::Int=100, # Number of siteordring sweeps
    order_degeneracy_check::Bool=true, # If true, then we consider degeneracy of the site ordering when it decide to swap sites order.
    order_convergence_count::Int=5, # Number of counts where orders are the same with those of previous sweeps sequentially to detect convergence of the site ordering.
    energy_convergence_count::Int=5, # Number of counts where energy is the same with that of previous sweeps sequentially to detect convergence of the energy.
    energy_tol::Float64=1e-8, # Tolerance for energy convergence.
    mps_cutoff::Float64=1e-13, # MPS truncation cutoff .
    mpo_cutoff::Float64=1e-13, # MPO truncation cutoff for emerging of the bond dimension introduced by swaps.
    proposer::AbstractSiteOrderProposer=MinimumEESiteOrderProposer(), # Site ordering proposer
    verbose::Int=1, # If true, then it will print the energy and site ordering at each sweep.
    dmrg_kwargs=(degeneracy_tol=1e-10), # DMRG kwargs,
)
    verbose > 0 && begin
        println("Site ordering settings:")
        println("  nswapsites: $nswapsites")
        println("  maxswapsweeps: $maxswapsweeps")
        println("  order_degeneracy_check: $order_degeneracy_check")
        println("  order_convergence_count: $order_convergence_count")
        println("  energy_convergence_count: $energy_convergence_count")
        println("  energy_tol: $energy_tol")
    end

    nswapsites == 1 && error("nswapsites must be greater than 1.")
    energy, psi = dmrg(H, psi, sweeps; dmrg_kwargs...) # Initial DMRG to get the initial energy and MPS.
    verbose > 0 && println("Initial energy: $energy")
    prev_energy = energy

    N = length(siteinds(psi))
    whole_order = collect(1:N) # Initial site ordering.
    prev_whole_order = copy(whole_order) # Previous site ordering.

    order_count = 0
    energy_count = 0

    order_flag = false # Flag to check if the site ordering is converged.
    energy_flag = false # Flag to check if the energy is converged.

    if sweeps_of_siteordering === nothing
        l = length(sweeps)
        sweeps_of_siteordering = Sweeps(1,
            ["maxdim"  "mindim"   "cutoff"      "noise";
            maxdim(sweeps, l)  mindim(sweeps, l)  cutoff(sweeps, l)  noise(sweeps, l)]
        )
    end

    for ni in 1:maxswapsweeps
        for (b, ha) in sweepnext(N)
            ha == 1 && b > N - nswapsites + 1 && continue
            ha == 2 && b < nswapsites && continue

            ortho = ha == 1 ? "left" : "right"
            sites = ortho == "left" ? collect(b:b+nswapsites-1) : collect(b-nswapsites+1:b)

            orthogonalize!(psi, sites[1])

            order = swap_order(psi, sites, ortho; proposer=proposer, order_degeneracy_check=order_degeneracy_check)

            before_order = copy(whole_order)
            if ortho == "left"
                whole_order[b:b+nswapsites-1] = whole_order[order]
            else
                whole_order[b-nswapsites+1:b] = whole_order[order]
            end
            psi = swap_sites(psi, before_order, whole_order; cutoff=mps_cutoff)
        end
        H = swap_sites(H, prev_whole_order, whole_order; cutoff=mpo_cutoff)
        energy, psi = dmrg(H, psi, sweeps_of_siteordering; dmrg_kwargs...)
        verbose > 0 && println("Energy at stage $ni: $energy")
        verbose > 1 && println("Site ordering at stage $ni: $whole_order")

        ΔE = abs(1.0 - energy/prev_energy)
        energy_converged = (ΔE < energy_tol)
        if energy_converged
            energy_count += 1
            if energy_count > energy_convergence_count
                energy_flag = true
            end
        else
            energy_count = 0
            energy_flag = false
        end
        prev_energy = energy

        # site ordering
        if prev_whole_order == whole_order
            order_count += 1
            if order_count > order_convergence_count
                order_flag = true
            end
        else
            order_count = 0
            order_flag = false
        end
        prev_whole_order = copy(whole_order)

        !order_degeneracy_check && energy_flag && order_count == 0 && begin
            order_degeneracy_check = true
            verbose > 1 && println("Order degeneracy check is now enabled, since energy is converged and site ordering is fluctuating.")
        end
        energy_flag && order_flag && break
    end
    return energy, psi, whole_order
end

function swap_sites(
    mp,
    sites::Vector{Int},
    order::Vector{Int};
    cutoff::Float64=1e-13,
)
    swaps = swap_sequence(sites, order)
    for swap in swaps
        typeof(mp) <: MPS && orthogonalize!(mp, swap)
        mp = swapbondsites(mp, swap)
        truncate!(mp, cutoff=cutoff)
    end
    return mp
end


function swap_sequence(initial::Vector{Int}, target::Vector{Int})
    sequence = collect(initial)
    swaps = Int[]

    pos_map = Dict(v => i for (i, v) in enumerate(sequence))

    for tgt_pos in 1:length(target)
        v = target[tgt_pos]
        cur_pos = pos_map[v]
        while cur_pos > tgt_pos
            push!(swaps, cur_pos - 1)
            sequence[cur_pos], sequence[cur_pos-1] =
                sequence[cur_pos-1], sequence[cur_pos]

            pos_map[ sequence[cur_pos]   ] = cur_pos
            pos_map[ sequence[cur_pos-1] ] = cur_pos - 1

            cur_pos -= 1
        end
    end

    return swaps
end


function swap_order(
    psi::MPS,
    sites::Vector{Int},
    ortho::String;
    proposer::AbstractSiteOrderProposer=MinimumEESiteOrderProposer(),
    order_degeneracy_check::Bool=false
)

    spectra_map = all_spectrums(psi, sites)
    order = propose_siteordering(proposer, spectra_map, sites; degeneracy_check=order_degeneracy_check)

    return order
end

"""
    all_spectrums(psi::MPS, sites::Vector{Int})

Return the vector containing spectrums at all bonds for all possible orders of the sites.
    The # of orders is given by the permutation of the number of sites.

    Inputs:
        - psi: MPS
        - sites: Vector{Int}

    Outputs:
        - spectra_map: Dict{Vector{Int}, Vector{Spectrum}}

"""
function all_spectrums(psi::MPS, sites::Vector{Int})
    φ          = contract(map(i -> psi[i], sites))
    site_inds  = map(i -> siteind(psi, i), sites)
    first_link = setdiff(inds(psi[sites[1]]),
                vcat(inds(psi[sites[2]]), siteind(psi, sites[1])))

    spectra_map = Dict{Vector{Int}, Vector{Spectrum}}()
    cache = Dict{Tuple{Int,Vector{Int},Tuple{UInt64}},
                 Tuple{ITensor,Spectrum,Index}}()

    function _rec(curφ, link,
                  rem_sites::Vector{Int}, rem_inds::Vector{<:Index},
                  path::Vector{Int}, specs::Vector{Spectrum})

        if length(rem_sites) == 0
            spectra_map[vcat(path, rem_sites)] = copy(specs)
            return
        end

        for (k, s_left) in pairs(rem_sites)
            rest_sites = deleteat!(copy(rem_sites), k)
            rest_inds  = deleteat!(copy(rem_inds),  k)

            key = (s_left, sort(rest_sites), (id(link[1]),))

            if !haskey(cache, key)
                Linds            = vcat(link, rem_inds[k])
                U,S,V,spec,u,_   = svd(curφ, Linds; cutoff = 0.0)
                cache[key]       = (S*V, spec, u)
            end

            Vblk, spec_here, new_link = cache[key]

            _rec(Vblk, [new_link],
                 rest_sites, rest_inds,
                 push!(copy(path), s_left),
                 push!(copy(specs), spec_here))
        end
    end

    _rec(φ, first_link, copy(sites), site_inds, Int[], Spectrum[])
    return spectra_map
end