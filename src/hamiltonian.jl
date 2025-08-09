function create_xxz_hamiltonian_from_file(sites, filename)
    coupling = read_xxz_data(filename)
    return create_xxz_hamiltonian(sites, coupling)
end

function create_xxz_hamiltonian(sites, coupling::Dict{Tuple{Int,Int}, Tuple{Float64, Float64}})
    
    os = OpSum()

    for (i, j) in keys(coupling)
        J_val, Δ_val = coupling[(i, j)]
        # Add XXZ terms
        if J_val != 0.0
            os += J_val, "Sx", i, "Sx", j
            os += J_val, "Sy", i, "Sy", j
            if Δ_val != 0.0
                os += J_val * Δ_val, "Sz", i, "Sz", j
            end
        end
    end

    return MPO(os, sites)
end

"""
    read_xxz_data(filename)

Read XXZ coupling data from a dat file.

# Arguments
- `filename`: Path to the data file

# File format
The data file should contain lines with the format:
```
i j J_{i,j} Δ_{i,j}
```
where i and j are site indices (1-based), and J_{i,j} and Δ_{i,j} are the coupling parameters.

# Returns
- `Dict{Tuple{Int,Int}, Tuple{Float64, Float64}}`: coupling dictionary

# Example file content:
```
1 2 1.0 0.5
2 3 1.0 0.5
3 4 1.0 0.5
4 5 1.0 0.5
5 6 1.0 0.5
```
"""
function read_xxz_data(filename)
    coupling = Dict{Tuple{Int,Int}, Tuple{Float64, Float64}}()

    open(filename, "r") do file
        for (line_num, line) in enumerate(eachline(file))
            # Skip empty lines and comments
            line = strip(line)
            if isempty(line) || startswith(line, "#")
                continue
            end

            # Parse the line
            try
                parts = split(line)
                if length(parts) != 4
                    error("Line $line_num: Expected 4 values (i j J_{i,j} Δ_{i,j}), got $(length(parts))")
                end

                i = parse(Int, parts[1])
                j = parse(Int, parts[2])
                J_val = parse(Float64, parts[3])
                Δ_val = parse(Float64, parts[4])
                # Validate site indices
                if i <= 0 || j <= 0
                    error("Line $line_num: Site indices must be positive integers")
                end

                ij = (i, j)
                if !haskey(coupling, ij)
                    coupling[ij] = (J_val, Δ_val)
                else
                    error("Line $line_num: Duplicate coupling found for site pair ($i, $j)")
                end

            catch e
                error("Line $line_num: Failed to parse line '$line': $e")
            end
        end
    end

    return coupling
end

