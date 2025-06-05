abstract type AbstractSiteOrderProposer end

struct MinimumEESiteOrderProposer <: AbstractSiteOrderProposer end

struct MinimumErrorSiteOrderProposer <: AbstractSiteOrderProposer end

function propose_siteordering(proposer::MinimumEESiteOrderProposer, spectrums::Dict{Vector{Int}, Vector{Spectrum}}, order::Vector{Int})
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
    entropy[findfirst(==(order), orders)] ≈ minimum(entropy) && return order
    return new_order
end

function propose_siteordering(proposer::MinimumErrorSiteOrderProposer, spectrums::Dict{Vector{Int}, Vector{Spectrum}}, order::Vector{Int})
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
    error[findfirst(==(order), orders)] ≈ minimum(error) && return order
    return new_order
end