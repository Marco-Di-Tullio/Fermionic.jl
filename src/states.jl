#states constructors
mutable struct State{T<:AbstractVector}
    st::T
    ope::Op
    nume::Union{Int64,Missing}
end

function State(st,ope;nume=missing)
    return State(st,ope,nume)
end

st(s::State) = s.st

ope(s::State) = s.ope

typ(s::State) = eltype(s.st)

nume(s::State) = s.nume

function rhosp(sta::State)
    precis = 15
    n = dim(ope(sta))
    num = nume(sta)
    rhospv = spzeros(typ(sta),n,n)
    estate = st(sta)
    o = ope(sta)
    if isequal(nume(sta),missing)
        for i in 1:n
            for j in 1:n
                rhospv[i,j] = round(estate'*cdcm(o, i, j)*estate, digits = precis)
            end
        end
    else
        if length(st(sta)) == binomial(n,num)
            for i in 1:n
                for j in 1:n
                    rhospv[i,j] = round(estate'*fixed(cdcm(o, i, j), num)*estate, digits = precis)
                end
            end
        else
            for i in 1:n
                for j in 1:n
                    rhospv[i,j] = round(estate'*cdcm(o, i, j)*estate, digits = precis)
                end
            end
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
