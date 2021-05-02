#states constructors

struct State_fixed{T<:AbstractVector}
    st::T
    ope::Op
    nume::Int64
end

st(s::State_fixed) = s.st

ope(s::State_fixed) = s.ope

nume(s::State_fixed) = s.nume

typ(s::State_fixed) = eltype(s.st)
