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

export Op, dim, basis, a, ad, ada, aad, aa, adad, vacuum
export State, st, ope, typ, rhosp, rhoqsp, non_zero
export eigensp, ssp, eigenqsp, sqsp, majorization_sp, majorization_qsp, n_avg, rhom, rhomd, trp
export sigma_x, sigma_y, sigma_z, phase, hadamard, ucnot, swap
export fixed, basis_m, fixed_state, unfixed_state, cdc, ccd, Op_fixed
export State_fixed, nume
export rhosp_mixed, eigensp_mixed
end # module
