#states constructors

struct State_fixed
    st::Array{Float64,1}
    ope::Op
    nume::Int64
end;

struct State_sparse_fixed
    st::SparseVector{Float64,Int64}
    ope::Op
    nume::Int64
end;

struct State_complex_fixed
    st::Array{Complex{Float64},1}
    ope::Op
    nume::Int64
end;

struct State_sparse_complex_fixed
    st::SparseVector{Complex{Float64},Int64}
    ope::Op
    nume::Int64
end;

st(s::State_fixed) = s.st
st(s::State_sparse_fixed) = s.st
st(s::State_complex_fixed) = s.st
st(s::State_sparse_complex_fixed) = s.st

ope(s::State_fixed) = s.ope
ope(s::State_sparse_fixed) = s.ope
ope(s::State_complex_fixed) = s.ope
ope(s::State_sparse_complex_fixed) = s.ope

nume(s::State_fixed) = s.nume
nume(s::State_sparse_fixed) = s.nume
nume(s::State_complex_fixed) = s.nume
nume(s::State_sparse_complex_fixed) = s.nume

function rhosp(sta::State_fixed)
    precis = 15
    n = dim(ope(sta))
    num = nume(sta)
    rhospv = zeros(n,n)
    estate = st(sta)
    o = ope(sta)
    for i in 1:n
        for j in 1:n
            rhospv[i,j] = round(estate'*fixed(cdcm(o, i, j), n, num)*estate, digits = precis)
        end
    end
    return rhospv
end

function rhosp(sta::State_sparse_fixed)
    precis = 15
    n = dim(ope(sta))
    num = nume(sta)
    rhosps = spzeros(n,n)
    estate = st(sta)
    o = ope(sta)
    for i in 1:n
        for j in 1:n
            rhosps[i,j] = round(estate'*fixed(cdcm(o, i, j), n, num)*estate, digits = precis)
        end
    end
    return rhosps
end

function rhosp(sta::State_complex_fixed)
    precis = 15
    n = dim(ope(sta))
    num = nume(sta)
    rhospvc = zeros(Complex,n,n)
    estate = st(sta)
    o = ope(sta)
    for i in 1:n
        for j in 1:n
            rhospvc[i,j] = round(estate'*fixed(cdcm(o, i, j), n, num)*estate, digits = precis)
        end
    end
    return rhospvc
end

function rhosp(sta::State_sparse_complex_fixed)
    precis = 15
    n = dim(ope(sta))
    num = nume(sta)
    rhospsc = spzeros(Complex{Float64},n,n)
    estate = st(sta)
    o = ope(sta)
    for i in 1:n
        for j in 1:n
            rhospsc[i,j] = round(estate'*fixed(cdcm(o, i, j), n, num)*estate, digits = precis)
        end
    end
    return rhospsc
end
