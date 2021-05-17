# Op_fixed initializes all the c^dagger c operators
# in the fixed particle number subspace
struct Op_fixed
    dim::Int
    nume::Int
    cdctot::SparseMatrixCSC{Float64,Int64}
    le::Int
    basis::SparseMatrixCSC{Float64,Int64}
    Op_fixed(dim,nume) = new(dim, nume, operators_fixed(dim,nume)...)
end

dim(o::Op_fixed) = o.dim
nume(o::Op_fixed) = o.nume
cdctot(o::Op_fixed) = o.cdctot
le(o::Op_fixed) = o.le
basis(o::Op_fixed) = o.basis

cdc(o::Op_fixed,i,j) = cdctot(o)[(1+(i-1)*le(o)):i*le(o),(1+(j-1)*le(o)):j*le(o)]
ccd(o::Op_fixed,i,j) = cdc(o,j,i)

# Here we create the full matrix with all
# fixed particle operators definitions
function operators_fixed(n::Int,m::Int)
    l1 = binomial(n,m)
    l2 = n*l1
    b,ind = basis_m(n,m)
    a=spzeros(l2,l2)
    for i in 1:n
        for j in i:n
            if i==j
                a[(1+(i-1)*l1):i*l1,(1+(j-1)*l1):j*l1] = 1/2*cdc(b,ind,i,j)
            else
                a[(1+(i-1)*l1):i*l1,(1+(j-1)*l1):j*l1] = cdc(b,ind,i,j)
            end
        end
    end
    op_general = a+a'
    return op_general, l1, b
end


#operator c^dagger c for a system with fixed particle number.
# Much faster than using fiexed over the whole operators

function cdc(n::Int64, m::Int64, i::Int64, j::Int64)
    base, ind = basis_m(n,m)
    indice = ind .-1
    l = binomial(n,m)
    op = spzeros(l,l)
    if i==j
        for k in 1:l
            if base[k,:][j] == 1
                op[k,k] = 1
            end
        end
    else
        for k in 1:l
            if base[k,:][j] == 1 && base[k,:][i] == 0
                ill = Int(indice[k]-2^(n-j)+2^(n-i))
                ill2 = myfind(indice, ill)[1]
                sign = (-1)^(sum(base[k,1:(j-1)])+sum(base[ill2,1:(i-1)]))
                op[ill2,k] = sign
            end
        end
    end
    return op
end

# When using many operators, it is easier to predefine the basis outside
# the function
function cdc(base::SparseArrays.SparseMatrixCSC{Float64,Int64}, ind::SparseArrays.SparseVector{Float64,Int64}, i::Int64, j::Int64)
    l = size(base)[1]
    n = size(base)[2]
    op = spzeros(l,l)
    indice = ind .-1
    if i==j
        for k in 1:l
            if base[k,:][j] == 1
                op[k,k] = 1
            end
        end
    else
        for k in 1:l
            if base[k,:][j] == 1 && base[k,:][i] == 0
                ill = Int(indice[k]-2^(n-j)+2^(n-i))
                ill2 = myfind(indice, ill)[1]
                sign = (-1)^(sum(base[k,1:(j-1)])+sum(base[ill2,1:(i-1)]))
                op[ill2,k] = sign
            end
        end
    end
    return op
end

# the c c^dagger operator, is the
function ccd(n::Int64, m::Int64, i::Int64, j::Int64)
    return cdc(n,m,j,i)
end

function ccd(base::SparseArrays.SparseMatrixCSC{Float64,Int64}, indice::SparseArrays.SparseVector{Float64,Int64}, i::Int64, j::Int64)
    return cdc(base,indice,j,i)
end


#This functions trimms the operators matrix
#for fixed particle number subspace
function fixed(o, nume::Int64)
    n = Int(log(2,size(o)[1]))
    posn = ffilter(n, nume)
    reduced = o[setdiff(1:end,posn), setdiff(1:end,posn)]
    return reduced
end

function fixed_state(stat::AbstractVector, nume::Int64)
    n = Int(log(2,length(stat)))
    posn = ffilter(n, nume)
    reduced = stat[setdiff(1:end,posn)]
    return reduced
end

function fixed_state(sta::State, nume::Int64)
    stat = st(sta)
    o = ope(sta)
    n = Int(log(2,length(stat)))
    posn = ffilter(n, nume)
    reduced = stat[setdiff(1:end,posn)]
    return State_fixed(reduced,o,nume)
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

function ffilter(n::Int64, m::Int64)
    _,_,baso,_ = integer_digits(n)
    posn = zeros(Int64, Int(2^n-binomial(n,m)))
    countern = 1
    for i in 1:2^n
        if sum(baso[i,:]) != m
            posn[countern] = Int(i)
            countern = countern + 1
        end
    end
    return posn
end

function myfind(c,j)
    a = similar(c, Int)
    count = 1
    @inbounds for i in eachindex(c)
        a[count] = i
        count += (c[i] == j)
    end
    return a[1:(count-1)]
end
