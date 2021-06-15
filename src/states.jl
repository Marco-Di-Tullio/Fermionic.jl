#states constructors
struct State{T<:AbstractVector}
    st::T
    ope::Op
end

st(s::State) = s.st

ope(s::State) = s.ope

typ(s::State) = eltype(s.st)

function rhosp(sta::State)
    precis = 15
    n = dim(ope(sta))
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

function kqsp(sta::State)
    precis = 15
    n = dim(ope(sta))
    k = spzeros(typ(sta),n,n)
    estate = st(sta)
    o = ope(sta)
    for i in 1:n
        for j in 1:n
            k[i,j] = round(estate'*aa(o, j, i)*estate, digits = precis)
        end
    end
    return k
end

rhoqsp(s::State) = [rhosp(s) kqsp(s); kqsp(s)' I-rhosp(s)']

function non_zero(c, prec=15)
    a = similar(c, Int)
    cof = similar(c, Float64)
    count = 1
    @inbounds for i in eachindex(c)
        a[count] = i
        cof[count] = round(c[i], digits = prec)
        count += (round(c[i], digits = prec) != zero(eltype(c)))
    end
    return resize!(a, count-1), resize!(cof, count-1)
end


#states constructors for both states_fixed and states
#disadvantage is that working with unmutables is less efficient
#also I believe it would be better to have both structs
# mutable struct State{T<:AbstractVector}
#     st::T
#     ope::Op
#     nume::Union{Int64,Missing}
# end
#
# function State(st,ope;nume=missing)
#     return State(st,ope,nume)
# end

# non_zero inputs a vector c and a precision (number of digits)
# and outputs two vectors: the first one indicating the indexes
# of non zero elements (up to precision prec) and the second
# the coefficient in each of these non zero indexes
