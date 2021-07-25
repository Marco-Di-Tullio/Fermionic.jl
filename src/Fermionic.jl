module Fermionic
using SparseArrays
using LinearAlgebra
using LazyStack

include("base_functions.jl")
include("operators.jl")
include("operators_fixed.jl")
include("states.jl")
include("states_fixed.jl")
include("correlations.jl")
include("logic_gates.jl")
include("mixed.jl")

export Op, dim, basis, a, ad, ada, aad, aa, adad, vacuum
export fixed, basis_m,  Op_fixed
export State, st, ope, typ, rhosp, rhoqsp, non_zero
export State_fixed, nume, fixed_state, unfixed_state, Op_semifixed
export eigensp, ssp, eigenqsp, sqsp, majorization_sp, majorization_qsp, n_avg, coef, rhom, rhomd, rhom2, trp
export sigma_x, sigma_y, sigma_z, phase, hadamard, ucnot, swap
export rhosp_mixed, eigensp_mixed
end # module
