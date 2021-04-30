#states constructors

struct State{T}
    st::Array{T,1}
    ope::Op
end

struct State_sparse{T}
    st::SparseVector{T,Int64}
    ope::Op
end

st(s::State) = s.st
st(s::State_sparse) = s.st

ope(s::State) = s.ope
ope(s::State_sparse) = s.ope

typ(s::State) = eltype(s.st)
typ(s::State_sparse) = eltype(s.st)

function rhosp(sta::State)
    precis = 15
    n = dim(ope(sta))
    rhospv = zeros(typ(sta),n,n)
    estate = st(sta)
    o = ope(sta)
    for i in 1:n
        for j in 1:n
            rhospv[i,j] = round(estate'*cdcm(o, i, j)*estate, digits = precis)
        end
    end
    return rhospv
end

function rhosp(sta::State_sparse)
    precis = 15
    n = dim(ope(sta))
    rhosps = spzeros(typ(sta),n,n)
    estate = st(sta)
    o = ope(sta)
    for i in 1:n
        for j in 1:n
            rhosps[i,j] = round(estate'*cdcm(o, i, j)*estate, digits = precis)
        end
    end
    return rhosps
end

function kqsp(sta::State)
    precis = 15
    n = dim(ope(sta))
    k = zeros(typ(sta),n,n)
    estate = st(sta)
    o = ope(sta)
    for i in 1:n
        for j in 1:n
            k[i,j] = round(estate'*cmcm(o, j, i)*estate, digits = precis)
        end
    end
    return k
end

function kqsp(sta::State_sparse)
    precis = 15
    n = dim(ope(sta))
    k = spzeros(typ(sta),n,n)
    estate = st(sta)
    o = ope(sta)
    for i in 1:n
        for j in 1:n
            k[i,j] = round(estate'*cmcm(o, j, i)*estate, digits = precis)
        end
    end
    return k
end

rhoqsp(s::State) = [rhosp(s) kqsp(s); kqsp(s)' I-rhosp(s)']
rhoqsp(s::State_sparse) = [rhosp(s) kqsp(s); kqsp(s)' I-rhosp(s)']
