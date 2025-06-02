module FlexibleDMRG

include("imports.jl")

include("mps.jl")
export replacebond!

include("dmrg.jl")
export dmrg

end