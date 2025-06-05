function siteordering(
    H::MPO,
    psi::MPS,
    sweeps::Sweeps;
    kwargs...
)
    N = length(siteinds(psi))
    whole_order = collect(1:N)
    prev_whole_order = copy(whole_order)
    while true
        H, psi, whole_order = siteordering(H, psi, whole_order)
        energy, psi = dmrg(H, psi, sweeps; kwargs...)
        if whole_order == prev_whole_order
            break
        end
        prev_whole_order = copy(whole_order)
    end
    return psi, H, whole_order
end

function siteordering(
    H::MPO,
    psi::MPS,
    whole_order::Vector{Int};
    nswapsites = 3, # that should be larger than 2.
    maxsweeps=100, # Number of sweeps
    observer=DMRGObserver(;energy_tol=1e-10, minsweeps=10),
    proposer::AbstractSiteOrderProposer=MinimumEESiteOrderProposer(),
)
    N = length(siteinds(psi))
    prev_whole_order = copy(whole_order)

    for ni in 1:maxsweeps
        for (b, ha) in sweepnext(N)
            ha == 1 && b > N - nswapsites + 1 && continue
            ni < maxsweeps && ha == 2 && b < nswapsites + 1 && continue
            ni == maxsweeps && ha == 2 && b < nswapsites && continue

            orthogonalize!(psi, b)

            ortho = ha == 1 ? "left" : "right"
            sites = ortho == "left" ? collect(b:b+nswapsites-1) : collect(b-nswapsites+1:b)

            order = swap_order(psi, sites, ortho; proposer=proposer)

            before_order = copy(whole_order)
            if ortho == "left"
                whole_order[b:b+nswapsites-1] = whole_order[order]
            else
                whole_order[b-nswapsites+1:b] = whole_order[order]
            end
            psi, H = swap_sites(psi, H, before_order, whole_order)
        end
        if prev_whole_order == whole_order
            break
        end
        prev_whole_order = copy(whole_order)
    end
    return H, psi, whole_order
end

function swap_sites(
    psi::MPS,
    H::MPO,
    sites::Vector{Int},
    order::Vector{Int},
)
    swaps = swap_sequence(sites, order)
    for swap in swaps
        orthogonalize!(psi, swap)
        psi = swapbondsites(psi, swap)
        H = swapbondsites(H, swap)
    end
    return psi, H
end


function swap_order(
    psi::MPS,
    sites::Vector{Int},
    ortho::String;
    proposer::AbstractSiteOrderProposer=MinimumEESiteOrderProposer(),
)
    spectra_map = all_spectrums(psi, sites)
    order = propose_siteordering(proposer, spectra_map, sites)
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
function all_spectrums(
    psi::MPS,
    sites::Vector{Int},
)
    function _recurse(
        cur_phi::ITensor,
        leftlinkind::Vector{Index{Int64}},
        rem_sites::Vector{Int},
        rem_inds::Vector{Index{Int64}},
        path_ids::Vector{Int},
        path_specs::Vector{Spectrum},
    )
        leftlinkind = typeof(leftlinkind) == Index ? [leftlinkind] : leftlinkind

        if length(rem_sites) == 1
            final_perm = copy(path_ids)
            push!(final_perm, rem_sites[1])
            spectra_map[final_perm] = copy(path_specs)
            return
        end

        for (k, new_left) in enumerate(rem_sites)
            rem2_sites = copy(rem_sites)
            deleteat!(rem2_sites, k)
            key = (new_left, sort(rem2_sites))

            if !haskey(cache, key)
                left_inds = vcat(rem_inds[k], leftlinkind)
                U, S, V, spec, u, v = svd(cur_phi, left_inds, cutoff=0.0)
                V = S*V
                cache[key] = (V, spec, u)
            end

            V_block, spec_here, v_block = cache[key]
            next_path_ids = copy(path_ids)
            push!(next_path_ids, new_left)
            next_path_specs = copy(path_specs)
            push!(next_path_specs, spec)

            rem2_inds = copy(rem_inds)
            rem2_sites  = copy(rem_sites)
            deleteat!(rem2_inds, k)
            deleteat!(rem2_sites,  k)

            leftlinkind = [v_block]

            _recurse(V_block, leftlinkind, rem2_sites, rem2_inds, next_path_ids, next_path_specs)
        end
    end

    phi = contract([psi[i] for i in sites])

    site_inds = [siteind(psi, i) for i in sites]

    leftlinkind = setdiff(inds(psi[sites[1]]), vcat(collect(inds(psi[sites[2]])), [siteind(psi, sites[1])]))

    spectra_map = Dict{Vector{Int}, Vector{Spectrum}}()

    cache = Dict{Tuple{Int, Vector{Int}}, Tuple{ITensor, Spectrum, Index{Int64}}}()

    _recurse(phi, leftlinkind, sites, site_inds, Vector{Int}(), Vector{Spectrum}())


    return spectra_map
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