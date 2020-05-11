#states constructors

struct State
    st::Array{Float64,1}
    ope::Op
end;

struct State_sparse
    st::SparseVector{Float64,Int64}
    ope::Op
end;

struct State_complex
    st::Array{Complex{Float64},1}
    ope::Op
end;

struct State_sparse_complex
    st::SparseVector{Complex{Float64},Int64}
    ope::Op
end;

st(s::State) = s.st
st(s::State_sparse) = s.st
st(s::State_complex) = s.st
st(s::State_sparse_complex) = s.st

ope(s::State) = s.ope
ope(s::State_sparse) = s.ope
ope(s::State_complex) = s.ope
ope(s::State_sparse_complex) = s.ope


function rhosp(sta::State)
    precis = 15
    n = dim(ope(sta))
    rhospv = zeros(n,n)
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
    rhosps = spzeros(n,n)
    estate = st(sta)
    o = ope(sta)
    for i in 1:n
        for j in 1:n
            rhosps[i,j] = round(estate'*cdcm(o, i, j)*estate, digits = precis)
        end
    end
    return rhosps
end

function rhosp(sta::State_complex)
    precis = 15
    n = dim(ope(sta))
    rhospvc = zeros(Complex,n,n)
    estate = st(sta)
    o = ope(sta)
    for i in 1:n
        for j in 1:n
            rhospvc[i,j] = round(estate'*cdcm(o, i, j)*estate, digits = precis)
        end
    end
    return rhospvc
end

function rhosp(sta::State_sparse_complex)
    precis = 15
    n = dim(ope(sta))
    rhospsc = spzeros(Complex{Float64},n,n)
    estate = st(sta)
    o = ope(sta)
    for i in 1:n
        for j in 1:n
            rhospsc[i,j] = round(estate'*cdcm(o, i, j)*estate, digits = precis)
        end
    end
    return rhospsc
end

function kqsp(sta::State)
    precis = 15
    n = dim(ope(sta))
    k = zeros(n,n)
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
    k = spzeros(n,n)
    estate = st(sta)
    o = ope(sta)
    for i in 1:n
        for j in 1:n
            k[i,j] = round(estate'*cmcm(o, j, i)*estate, digits = precis)
        end
    end
    return k
end

function kqsp(sta::State_complex)
    precis = 15
    n = dim(ope(sta))
    k = zeros(Complex{Float64},n,n)
    estate = st(sta)
    o = ope(sta)
    for i in 1:n
        for j in 1:n
            k[i,j] = round(estate'*cmcm(o, j, i)*estate, digits = precis)
        end
    end
    return k
end

function kqsp(sta::State_sparse_complex)
    precis = 15
    n = dim(ope(sta))
    k = spzeros(Complex{Float64},n,n)
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
rhoqsp(s::State_complex) = [rhosp(s) kqsp(s); kqsp(s)' I-rhosp(s)']
rhoqsp(s::State_sparse) = [rhosp(s) kqsp(s); kqsp(s)' I-rhosp(s)']
rhoqsp(s::State_sparse_complex) = [rhosp(s) kqsp(s); kqsp(s)' I-rhosp(s)']
