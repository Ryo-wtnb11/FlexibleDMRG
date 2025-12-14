using Combinatorics
include("outputs.jl")


function sites_hamiltonian(
    sites::Vector{<:Index},
    interactions::Dict{<:Union{Tuple{Int,Int}, Int},
                        Vector{Tuple{Float64, Tuple{String, String}}}},
    whole_order::Vector{Int}
)
    particle_to_site = Dict{Int, Int}()
    for (site_idx, particle_id) in enumerate(whole_order)
        particle_to_site[particle_id] = site_idx
    end

    os = OpSum()
    
    for (particle_key, ham) in interactions
        if particle_key isa Tuple{Int, Int} # 2-body interaction
            pid_i, pid_j = particle_key

            site_i = particle_to_site[pid_i]
            site_j = particle_to_site[pid_j]

            for h in ham
                coefficient = first(h)
                op_pair = last(h)
                os += coefficient, op_pair[1], site_i, op_pair[2], site_j
            end

        elseif particle_key isa Int # on-site interaction
            pid_i = particle_key

            site_i = particle_to_site[pid_i]

            for h in ham
                coefficient = first(h)
                op_pair = last(h)

                os += coefficient, op_pair[1], site_i
            end
        end
    end

    H = MPO(os, sites)
    return H
end


function siteordering(
    psi::MPS,
    sweeps::Sweeps,
    interactions::Dict{<:Union{Tuple{Int,Int}, Int},
                        Vector{Tuple{Float64, Tuple{String, String}}}};
    sweeps_of_siteordering=nothing, # Sweeps for DMRG in site ordering. If nothing, then it will be set to the last sweep of the sweeps.
    nswapsites::Int=3, # that should be around 3 - 4.
    nsets::Int=100, # Number of siteordring sweeps
    order_degeneracy_check::Bool=true, # If true, then we consider degeneracy of the site ordering when it decide to swap sites order.
    convergence_count::Int=5, # Number of counts where energy is the same with that of previous sweeps sequentially to detect convergence of the energy.
    energy_tol::Float64=1e-10, # Tolerance for energy convergence.
    mps_cutoff::Float64=1e-15, # MPS truncation cutoff .
    mpo_cutoff::Float64=1e-15, # MPO truncation cutoff for emerging of the bond dimension introduced by swaps.
    proposer::AbstractSiteOrderProposer=MinimumEESiteOrderProposer(), # Site ordering proposer
    verbose::Int=1, # If true, then it will print the energy and site ordering at each sweep.
    dmrg_kwargs=(degeneracy_tol=1e-10,), # DMRG kwargs,
    save_outputs::Bool = true,
)
    verbose > 0 && begin
        println("Site ordering settings:")
        println("  nswapsites: $nswapsites")
        println("  nsets: $nsets")
        println("  order_degeneracy_check: $order_degeneracy_check")
        println("  convergence_count: $convergence_count")
        println("  energy_tol: $energy_tol")
    end

    observer=DMRGObserver(energy_tol=energy_tol)

    N = length(siteinds(psi))
    whole_order = collect(1:N) # Initial site ordering.
    prev_whole_order = copy(whole_order) # Previous site ordering.

    nswapsites == 1 && error("nswapsites must be greater than 1.")

    H = sites_hamiltonian(siteinds(psi), interactions, whole_order)
    energy, psi = dmrg(H, psi, sweeps; observer=observer, dmrg_kwargs...) # Initial DMRG to get the initial energy and MPS.

    maxdim = maxlinkdim(psi)

    verbose > 0 && println("Initial energy: $energy")
    prev_energy = energy

    n_converged = 0

    if sweeps_of_siteordering === nothing
        l = length(sweeps)
        sweeps_of_siteordering = sweeps
    end

    for ni in 1:nsets
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

            psi = swap_mps(psi, sites, order; cutoff=mps_cutoff, maxdim=maxdim)

        end
        
        # TODO: It is easy to generalize sites_hamiltonian to the case of higher-order interactions.
        H = sites_hamiltonian(siteinds(psi), interactions, whole_order)
        
        # we have not implemented sites_hamiltonian for the Hamiltonian with more than 2-body interactions.
        # H = swap_sites(H, prev_whole_order, whole_order; cutoff=mpo_cutoff) # we can use this any MPO.
        # NOTE: swap_sites is quit slow, and it may collapse the original Hamiltonian via applying many swap operators.

        energy, psi = dmrg(H, psi, sweeps_of_siteordering; observer=observer, dmrg_kwargs...)

        ΔE = abs(1.0 - energy/prev_energy)
        energy_converged = (ΔE < energy_tol)
        order_converged = (prev_whole_order == whole_order)

        if energy_converged && order_converged
            n_converged += 1
            if n_converged >= convergence_count
                break
            end
        else
            n_converged = 0
        end

        if ni == nsets
            break
        end

        prev_energy = energy
        prev_whole_order = copy(whole_order)
    end

    return energy, psi, whole_order
end

function swap_mps(
    psi::MPS,
    sites::Vector{Int},
    order::Vector{Int};
    cutoff::Float64=1e-13,
    maxdim::Int=maxlinkdim(psi),
)
    Ns = length(sites)
    
    φ = reduce(*, psi[sites])
    
    site_inds = [siteind(psi, s) for s in order]

    first_link = (sites[1] > 1) ? commonind(psi[sites[1]], psi[sites[1]-1]) : nothing
    L_link = first_link

    # 4. Sequential SVD
    for k in 1:(Ns - 1)
        L_inds = Index[site_inds[k]]

        if !isnothing(L_link)
            push!(L_inds, L_link)
        end

        U, S, V = svd(φ, L_inds; cutoff=cutoff, maxdim=maxdim, lefttags="Link,u")
        
        psi[sites[k]] = U

        φ = S * V
        L_link = commonind(U, S)
    end

    # Last site
    psi[sites[end]] = φ
    orthogonalize!(psi, sites[end])
    
    return psi
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
function all_spectrums_legacy(psi::MPS, sites::Vector{Int})
    φ = contract(map(i -> psi[i], sites))
    site_inds = map(i -> siteind(psi, i), sites)
    first_link = setdiff(
        inds(psi[sites[1]]),
        vcat(inds(psi[sites[2]]), siteind(psi, sites[1]))
    )

    spectra_map = Dict{Vector{Int}, Vector{Spectrum}}()
    cache = Dict{Tuple{Int,Vector{Int},Tuple{UInt64}},
                 Tuple{ITensor,Spectrum,Index}}()

    function _rec(
        curφ,
        link,
        rem_sites::Vector{Int},
        rem_inds::Vector{<:Index},
        path::Vector{Int},
        specs::Vector{Spectrum}
    )
        if length(rem_sites) == 0
            spectra_map[vcat(path, rem_sites)] = copy(specs)
            return
        end

        for (k, s_left) in pairs(rem_sites)
            rest_sites = deleteat!(copy(rem_sites), k)
            rest_inds  = deleteat!(copy(rem_inds),  k)

            key = (s_left, sort(rest_sites), (id(link[1]),))

            if !haskey(cache, key)
                Linds = vcat(link, rem_inds[k])
                U, S, V, spec, u, _ = svd(curφ, Linds; cutoff = 0.0)
                cache[key] = (S * V, spec, u)
            end

            Vblk, spec_here, new_link = cache[key]

            _rec(
                Vblk,
                [new_link],
                rest_sites,
                rest_inds,
                push!(copy(path), s_left),
                push!(copy(specs), spec_here)
            )
        end
    end

    _rec(φ, first_link, copy(sites), site_inds, Int[], Spectrum[])
    return spectra_map
end

function build_spectrum_lookup(psi::MPS, sites::Vector{Int})
    Ns = length(sites)
    
    # 1. Contract the entire region
    φ = contract([psi[s] for s in sites])
    
    site_inds_map = Dict(s => siteind(psi, s) for s in sites)
    
    # 2. Identify the left boundary indices (L_boundary)
    # Identify the boundary indices of the entire tensor φ that are connected to the first site (sites[1])
    all_phys_inds = [siteind(psi, s) for s in sites]
    cluster_boundary_inds = uniqueinds(inds(φ), all_phys_inds)
    L_boundary = commoninds(cluster_boundary_inds, inds(psi[sites[1]]))

    lookup_table = Dict{Vector{Int}, Spectrum}()

    # 3. Loop over subsets
    for i in 1:Ns
        for subset in combinations(sites, i)
            key = sort(subset)
            # Left side of SVD: [physical indices of subset] + [left boundary link]
            Linds = [site_inds_map[s] for s in subset]
            append!(Linds, L_boundary)

            # Execute SVD
            _, S, _ = svd(φ, Linds; cutoff=0.0)
            
            # [Fix] ITensors.svd returns a Spectrum with "squared singular values (eigenvalues)" by specification.
            # Since the legacy code uses it as-is, we need to square them here as well for storage.
            dims_S = minimum(dims(S))
            eigenvalues = Float64[S[k, k]^2 for k in 1:dims_S]
            
            # Create Spectrum with truncation error (truncerr) set to 0.0
            lookup_table[key] = Spectrum(eigenvalues, 0.0)
        end
    end

    return lookup_table
end

"""
    evaluate_all_permutations_using_lookup(sites::Vector{Int}, lookup_table::Dict)
"""
function evaluate_all_permutations_using_lookup(sites::Vector{Int}, lookup_table::Dict)
    Ns = length(sites)
    results = Dict{Vector{Int}, Vector{Spectrum}}()

    for p in permutations(sites)
        spectra_sequence = Vector{Spectrum}()
        sizehint!(spectra_sequence, Ns)
        
        current_left_set = Int[]
        for k in 1:Ns
            push!(current_left_set, p[k])
            key = sort(current_left_set)
            push!(spectra_sequence, lookup_table[key])
        end
        results[p] = spectra_sequence
    end
    return results
end

# optimized version of all_spectrums
function all_spectrums(psi::MPS, sites::Vector{Int})
    table = build_spectrum_lookup(psi, sites)
    return evaluate_all_permutations_using_lookup(sites, table)
end