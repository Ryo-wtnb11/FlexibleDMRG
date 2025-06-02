"""
    replacebond!(M::MPS, b::Int, phi::ITensor; kwargs...)

Factorize the ITensor `phi` and replace the ITensors
`b` and `b+1` of MPS `M` with the factors. Choose
the orthogonality with `ortho="left"/"right"`.
"""
function replacebond!(
    M::MPS,
    b::Int,
    phi::ITensor;
    normalize=nothing,
    swapsites=nothing,
    ortho=nothing,
    # Decomposition kwargs
    which_decomp=nothing,
    mindim=nothing,
    maxdim=nothing,
    cutoff=nothing,
    eigen_perturbation=nothing,
    # svd kwargs
    svd_alg=nothing,
    use_absolute_cutoff=nothing,
    use_relative_cutoff=nothing,
    min_blockdim=nothing,
    )
    normalize = NDTensors.replace_nothing(normalize, false)
    swapsites = NDTensors.replace_nothing(swapsites, false)
    ortho = NDTensors.replace_nothing(ortho, "left")

    indsMb = inds(M[b])
    if swapsites
        sb = siteind(M, b)
        sbp1 = siteind(M, b + 1)
        indsMb = replaceind(indsMb, sb, sbp1)
    end

    L, R, spec = factorize(
        phi,
        indsMb;
        mindim,
        maxdim,
        cutoff,
        ortho,
        which_decomp,
        eigen_perturbation,
        svd_alg,
        tags=tags(linkind(M, b)),
        use_absolute_cutoff,
        use_relative_cutoff,
        min_blockdim,
    )
    M[b] = L
    M[b + 1] = R
    if ortho == "left"
        leftlim(M) == b - 1 && setleftlim!(M, leftlim(M) + 1)
        rightlim(M) == b + 1 && setrightlim!(M, rightlim(M) + 1)
        normalize && (M[b + 1] ./= norm(M[b + 1]))
    elseif ortho == "right"
        leftlim(M) == b && setleftlim!(M, leftlim(M) - 1)
        rightlim(M) == b + 2 && setrightlim!(M, rightlim(M) - 1)
        normalize && (M[b] ./= norm(M[b]))
    else
        error(
        "In replacebond!, got ortho = $ortho, only currently supports `left` and `right`."
        )
    end
    return spec
    end
