module Fermionic
using SparseArrays
using LinearAlgebra

include("base_functions.jl")
include("operators.jl")
include("states.jl")
include("states_fixed.jl")
include("correlations.jl")
include("logic_gates.jl")
include("operators_fixed.jl")
include("mixed.jl")

export Op, dim, basis, cm, cdm, cdcm, cmcd, cmcm, cdcd, vacuum
export State, State_sparse, State_complex, State_sparse_complex, st, ope, rhosp, rhoqsp
export eigensp, ssp, eigenqsp, sqsp, majorization_sp, majorization_qsp, n_avg, rhom, rhomnd
export sigma_x, sigma_y, sigma_z, phase, hadamard, ucnot, swap
export fixed, basis_m, fixed_state, cdc
export State_fixed, State_sparse_fixed, State_complex_fixed, State_sparse_complex_fixed, nume
export rhosp_mixed, eigensp_mixed
end # module
