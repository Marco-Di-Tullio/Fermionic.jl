module Fermionic
using SparseArrays
using LinearAlgebra

include("extra_file.jl")
include("base_functions.jl")
include("operators.jl")
include("states.jl")

export my_f, my_g

end # module
