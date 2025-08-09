using ITensors, ITensorMPS
using FlexibleDMRG: siteordering, MinimumEESiteOrderProposer
include("../src/hamiltonian.jl")
include("../src/parameters.jl")



function main()
    base_path = ARGS[1]
    parameters_path = joinpath(base_path, "parameters.toml")

    params = read_parameters(parameters_path)
    params = validate_parameters(params)

    output_path = joinpath(base_path, "outputs")
    mkpath(output_path)

    sweeps = create_sweeps_from_params(params)
    dmrg_kwargs = create_dmrg_kwargs_from_params(params)
    kwargs = create_siteordering_kwargs_from_params(params)

    N = params["system"]["system_size"]
    sitetype = params["system"]["site_type"]
    sites = siteinds(sitetype, N)

    # Heisenberg model
    H = create_xxz_hamiltonian_from_file(sites, joinpath(base_path, params["hamiltonian"]["coupling_file"]))

    # Create an initial random matrix product state
    psi0 = random_mps(sites)
    # Run the DMRG algorithm, returning energy

    output_keys = ["maxtruncerr", "maxlinkdim", "linkEE", "siteEE", "Sz"]
    # (dominant eigenvalue) and optimized MPS
    energy, psi, whole_order = siteordering(H, psi0, sweeps; kwargs..., dmrg_kwargs=dmrg_kwargs, output_path=output_path, output_keys=output_keys)
    
    return 0
end

main()