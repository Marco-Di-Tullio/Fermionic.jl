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
            rhospv[i,j] = round(estate'*cdcm(o, i, j)*estate, digits = precis)
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
            k[i,j] = round(estate'*cmcm(o, j, i)*estate, digits = precis)
        end
    end
    return k
end

rhoqsp(s::State) = [rhosp(s) kqsp(s); kqsp(s)' I-rhosp(s)']
