module Fermionic
using SparseArrays
using LinearAlgebra

include("extra_file.jl")
include("base_functions.jl")
include("operators.jl")
include("states.jl")

export my_f, my_g
export  Op, dim, basis, cm, cdm, cdcm, cmcd, cmcm, cdcd
export State, State_sparse, st, ope, rhosp, eigensp, ssp
end # module
