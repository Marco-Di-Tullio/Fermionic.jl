#states constructors

struct State
    st::Array{Float64,1}
    ope::Op
end

struct State_sparse
    st::SparseVector{Float64,Int64}
    ope::Op
end

st(s::State) = s.st
st(s::State_sparse) = s.st

ope(s::State) = s.ope
ope(s::State_sparse) = s.ope


function rhosp(sta::State)
    precis = 15
    n = dim(ope(sta))
    rhospv = zeros(n,n)
    for i in 1:n
        for j in 1:n
            rhospv[i,j] = round(st(sta)'*cdcm(ope(sta), i, j)*st(sta), digits = precis)
        end
    end
    return rhospv
end

function rhosp(sta::State_sparse)
    precis = 15
    n = dim(ope(sta))
    rhosps = spzeros(n,n)
    for i in 1:n
        for j in 1:n
            rhosps[i,j] = round(st(sta)'*cdcm(ope(sta), i, j)*st(sta), digits = precis)
        end
    end
    return rhosps
end

prec = 15 #precision
eigensp(s::State) = [round(eigvals(rhosp(s))[i], digits = prec) for i in 1:dim(ope(s))]
#For some reason, I can not compute eigenvalues
#directly from sparse matrices
eigensp(s::State_sparse) = [round(eigvals(Matrix(rhosp(s)))[i], digits = prec) for i in 1:dim(ope(s))]

function ssp(sta::State)
    eigen = eigensp(sta)
    lene = length(eigen)
    s = 0
    for i in 1:lene
        if eigen[i] != 0 && eigen[i] != 1
            s = s - (eigen[i]*log(2,eigen[i]) + (1 - eigen[i])*log(2,1-eigen[i]))
        end
    end
    return s/lene
end

function ssp(sta::State_sparse)
    eigen = eigensp(sta)
    lene = length(eigen)
    s = 0
    for i in 1:lene
        if eigen[i] != 0 && eigen[i] != 1
            s = s - (eigen[i]*log(2, eigen[i]) + (1 - eigen[i])*log(2, 1-eigen[i]))
        end
    end
    return s/lene
end
