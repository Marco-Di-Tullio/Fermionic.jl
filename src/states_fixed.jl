#states constructors

struct State_fixed{T}
    st::Array{T,1}
    ope::Op
    nume::Int64
end

struct State_sparse_fixed{T}
    st::SparseVector{T,Int64}
    ope::Op
    nume::Int64
end

st(s::State_fixed) = s.st
st(s::State_sparse_fixed) = s.st

ope(s::State_fixed) = s.ope
ope(s::State_sparse_fixed) = s.ope

nume(s::State_fixed) = s.nume
nume(s::State_sparse_fixed) = s.nume

typ(s::State_fixed) = eltype(s.st)
typ(s::State_sparse_fixed) = eltype(s.st)

function rhosp(sta::State_fixed)
    precis = 15
    n = dim(ope(sta))
    num = nume(sta)
    rhospv = zeros(typ(sta),n,n)
    estate = st(sta)
    o = ope(sta)
    for i in 1:n
        for j in 1:n
            rhospv[i,j] = round(estate'*fixed(cdcm(o, i, j), num)*estate, digits = precis)
        end
    end
    return rhospv
end

function rhosp(sta::State_sparse_fixed)
    precis = 15
    n = dim(ope(sta))
    num = nume(sta)
    rhosps = spzeros(typ(sta),n,n)
    estate = st(sta)
    o = ope(sta)
    for i in 1:n
        for j in 1:n
            rhosps[i,j] = round(estate'*fixed(cdcm(o, i, j), num)*estate, digits = precis)
        end
    end
    return rhosps
end
