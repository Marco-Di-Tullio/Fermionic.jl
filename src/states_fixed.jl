#states constructors

struct State_fixed{T<:AbstractVector}
    st::T
    ope::Op_fixed
end
#
# struct State_fixed{T<:AbstractVector}
#     st::T
#     ope::Op
#     nume::Int64
# end

st(s::State_fixed) = s.st

ope(s::State_fixed) = s.ope

nume(s::State_fixed) = nume(s.ope)

typ(s::State_fixed) = eltype(s.st)

function rhosp(sta::State_fixed)
    precis = 15
    n = dim(ope(sta))
    num = nume(sta)
    rhospv = spzeros(typ(sta),n,n)
    estate = st(sta)
    o = ope(sta)
    for i in 1:n
        for j in 1:n
            rhospv[i,j] = round(estate'*ada(o, i, j)*estate, digits = precis)
        end
    end
    return rhospv
end

# This functions takes stat, an array
# in the full 2^n space, into the reduced
# space.
function fixed_state(stat::AbstractVector, nume::Int64)
    n = Int(log(2,length(stat)))
    posn = ffilter(n, nume)
    reduced = stat[setdiff(1:end,posn)]
    return reduced
end

function fixed_state(sta::State, nume::Int64)
    stat = st(sta)
    n = Int(log(2,length(stat)))
    posn = ffilter(n, nume)
    reduced = stat[setdiff(1:end,posn)]
    return State_fixed(reduced,Op_fixed(n,nume))
end

# This functions takes stat, an array
# in the reduced space, into the full 2^n
# space. Index is the index from basis_m
# and d the dimension
function unfixed_state(stat::AbstractVector, d::Int64, m::Int64)
    index = basis_m(d,m)[2]
    l = size(stat)[1]
    state2n = zeros(2^d)
    for (i,v) in enumerate(index)
        state2n[Int(v)]=stat[i]
    end
    return state2n
end

function unfixed_state(sta::State_fixed)
    d = dim(ope(sta))
    index = basis_m(d,nume(sta))[2]
    stat = st(sta)
    l = size(stat)[1]
    state2n = zeros(2^d)
    for (i,v) in enumerate(index)
        state2n[Int(v)]=stat[i]
    end
    return State(state2n,Op(d))
end
