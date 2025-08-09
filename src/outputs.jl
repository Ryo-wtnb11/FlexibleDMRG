using JSON3
using DataFrames, CSV

function write_outputs(energy, whole_order, outputs, path::String)
    # ==== 1) core.json ====
    # Keep it simple: write compact JSON (one line, whole_order inline)
    core = (
        energy = energy,
        maxlinkdim = get(outputs, "maxlinkdim", nothing),
        maxtruncerr = get(outputs, "maxtruncerr", nothing),
    )
    open(joinpath(path, "core.json"), "w") do io
        JSON3.pretty(io, core)
        write(io, '\n')
    end

    # ==== 2) whole_order.json ====
    open(joinpath(path, "order.dat"), "w") do io
        write(io, join(whole_order, " "))
    end

    # ==== 2) sites.csv ====
    N = length(outputs["siteEE"])
    Sx_col = haskey(outputs, "Sx") ? outputs["Sx"] : fill(missing, N)
    Sy_col = haskey(outputs, "Sy") ? outputs["Sy"] : fill(missing, N)
    Sz_col = haskey(outputs, "Sz") ? outputs["Sz"] : fill(missing, N)

    sites_df = DataFrame(
        "site" => 1:N,
        "siteEE" => outputs["siteEE"],
        "Sx" => Sx_col,
        "Sy" => Sy_col,
        "Sz" => Sz_col,
    )
    CSV.write(joinpath(path, "sites.csv"), sites_df)

    # ==== 3) links.csv ====
    linkEE = haskey(outputs, "linkEE") ? outputs["linkEE"] : fill(missing, N-1)
    links_df = DataFrame(
        "site_i" => 1:N-1,
        "site_j" => 2:N,
        "linkEE" => linkEE,
    )
    CSV.write(joinpath(path, "links.csv"), links_df)
end