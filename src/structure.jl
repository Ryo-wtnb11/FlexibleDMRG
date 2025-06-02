function structure_candidate(M::MPS, b::Int, phi::ITensor)
    U, S, V, spec, u, v = svd(phi, [inds(M[b])], cutoff=0.0)
    return
end