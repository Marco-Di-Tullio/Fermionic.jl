#states constructors

struct State_fixed{T<:AbstractVector}
    st::T
    ope::Op
    nume::Int64
end

struct State_sparse_fixed{T}
    st::SparseVector{T,Int64}
    ope::Op
    nume::Int64
end

st(s::State_fixed) = s.st

ope(s::State_fixed) = s.ope

nume(s::State_fixed) = s.nume

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
            rhospv[i,j] = round(estate'*fixed(cdcm(o, i, j), num)*estate, digits = precis)
        end
    end
    return rhospv
end
