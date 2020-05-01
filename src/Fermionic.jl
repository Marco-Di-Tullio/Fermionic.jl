module Fermionic
using SparseArrays
using LinearAlgebra

include("base_functions.jl")
include("operators.jl")
include("states.jl")
include("logic_gates.jl")

export  Op, dim, basis, cm, cdm, cdcm, cmcd, cmcm, cdcd, vacuum
export State, State_sparse, st, ope, rhosp, eigensp, ssp
export ucnot, hadamard, not, swap, phase
end # module
