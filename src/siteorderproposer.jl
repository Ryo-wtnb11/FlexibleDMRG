abstract type AbstractSiteOrderProposer end

struct MinimumEESiteOrderProposer <: AbstractSiteOrderProposer end

struct MinimumErrorSiteOrderProposer <: AbstractSiteOrderProposer end

function propose_siteordering(
    proposer::MinimumEESiteOrderProposer,
    spectrums::Dict{Vector{Int}, Vector{Spectrum}},
    order::Vector{Int};
    degeneracy_check::Bool=false,
)
    entropy = Vector{Float64}(undef, length(keys(spectrums)))
    orders = collect(keys(spectrums))
    for (i, order) in enumerate(orders)
        specs = spectrums[order]
        entropy[i] = 0.0
        for spec in specs
            entropy[i] += sum(-spec.eigs .* log.(spec.eigs))
        end
    end
    new_order = orders[argmin(entropy)]
    if degeneracy_check
        entropy[findfirst(==(order), orders)] ≈ minimum(entropy) && return order
    end
    return new_order
end

function propose_siteordering(
    proposer::MinimumErrorSiteOrderProposer,
    spectrums::Dict{Vector{Int}, Vector{Spectrum}},
    order::Vector{Int};
    degeneracy_check::Bool=false,
)
    error = Vector{Float64}(undef, length(keys(spectrums)))
    orders = collect(keys(spectrums))
    for (i, order) in enumerate(orders)
        specs = spectrums[order]
        error[i] = 0.0
        for spec in specs
            error[i] += 1.0 - sum(spec.eigs)
        end
    end
    new_order = orders[argmin(error)]
    if degeneracy_check
        error[findfirst(==(order), orders)] ≈ minimum(error) && return order
    end
    return new_order
end