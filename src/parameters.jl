using TOML
using FlexibleDMRG: AbstractSiteOrderProposer, MinimumEESiteOrderProposer

"""
    read_parameters(filename)

Read parameters from a TOML file.

# Arguments
- `filename`: Path to the parameter file

# Returns
- `Dict`: Dictionary containing all parameters

# Example TOML file structure:
```toml
[system]
spin_size = 0.5
system_size = 10

[hamiltonian]
coupling_file = "couplings.dat"

[dmrg]
max_bond_dim = 100
energy_tolerance = 1e-8
convergence_sets = 3
max_sweeps = 10
max_sets = 20

[output]
filename = "results"
folder = "output"
```
"""
function read_parameters(filename)
    if !isfile(filename)
        error("Parameter file not found: $filename")
    end
    
    params = TOML.parsefile(filename)
    
    # Validate required sections
    required_sections = ["system", "hamiltonian", "dmrg"]
    for section in required_sections
        if !haskey(params, section)
            error("Missing required section: $section")
        end
    end
    
    return params
end

"""
    validate_parameters(params)

Validate parameter values and set defaults where appropriate.

# Arguments
- `params`: Parameter dictionary from read_parameters

# Returns
- `Dict`: Validated parameters with defaults applied
"""
function validate_parameters(params)
    validated = deepcopy(params)

    # System parameters
    if !haskey(validated["system"], "site_type")
        error("Missing site_type in system section")
    end
    if !haskey(validated["system"], "system_size")
        error("Missing system_size in system section")
    end

    # Hamiltonian parameters
    if !haskey(validated["hamiltonian"], "coupling_file")
        error("Missing coupling_file in hamiltonian section")
    end

    # DMRG parameters with defaults
    dmrg = validated["dmrg"]
    if !haskey(dmrg, "nsweeps")
        dmrg["nsweeps"] = 10
        warn("nsweeps not set in parameters.toml, using default value of 1")
    end
    if !haskey(dmrg, "maxdim")
        dmrg["maxdim"] = 50
        warn("maxdim not set in parameters.toml, using default value of 50")
    end
    if !haskey(dmrg, "cutoff")
        dmrg["cutoff"] = 1e-15
        warn("cutoff not set in parameters.toml, using default value of 1e-15")
    end
    if !haskey(dmrg, "noise")
        dmrg["noise"] = 0.0
        warn("noise not set in parameters.toml, using default value of 0.0")
    end
    if !haskey(dmrg, "nsets")
        dmrg["nsets"] = 20
        warn("nsets not set in parameters.toml, using default value of 10")
    end
    if !haskey(dmrg, "convergence_count")
        dmrg["convergence_count"] = 5
        warn("convergence_count not set in parameters.toml, using default value of 5")
    end
    if !haskey(dmrg, "energy_tol")
        dmrg["energy_tol"] = 1e-10
        warn("energy_tol not set in parameters.toml, using default value of 1e-10")
    end
    if !haskey(dmrg, "degeneracy_tol")
        dmrg["degeneracy_tol"] = 1e-10
        warn("degeneracy_tol not set in parameters.toml, using default value of 1e-10")
    end
    if !haskey(dmrg, "nswapsites")
        dmrg["nswapsites"] = 2
        warn("nswapsites not set in parameters.toml, using default value of 2")
    end

    return validated
end

function create_sweeps_from_params(params)
    dmrg = params["dmrg"]

    # Create a single sweep with the specified parameters
    sweeps = Sweeps(
        dmrg["nsweeps"],
        [
            "maxdim" "cutoff" "noise";
            dmrg["maxdim"] dmrg["cutoff"] dmrg["noise"]
        ]
    )

    return sweeps
end

function create_dmrg_kwargs_from_params(params)
    dmrg = params["dmrg"]
    return (
        degeneracy_tol = dmrg["degeneracy_tol"],
    )
end

function create_siteordering_kwargs_from_params(params)
    dmrg = params["dmrg"]
    return (
        nswapsites = dmrg["nswapsites"],
        nsets = dmrg["nsets"],
        convergence_count = dmrg["convergence_count"],
        energy_tol = dmrg["energy_tol"],
        dmrg_kwargs = create_dmrg_kwargs_from_params(params),
    )
end