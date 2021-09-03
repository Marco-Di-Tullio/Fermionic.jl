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

ada(o::Op_fixed,i,j) = cdctot(o)[(1+(i-1)*le(o)):i*le(o),(1+(j-1)*le(o)):j*le(o)]
aad(o::Op_fixed,i,j) = ada(o,j,i)

# Here we create the full matrix with all
# fixed particle operators definitions
function operators_fixed(n::Int,m::Int)
    l1 = binomial(n,m)
    l2 = n*l1
    b,ind = basis_m(n,m)
    z = spzeros(l2,l2)
    for i in 1:n
        for j in i:n
            if i==j
                z[(1+(i-1)*l1):i*l1,(1+(j-1)*l1):j*l1] = 1/2*ada(b,ind,i,j)
            else
                z[(1+(i-1)*l1):i*l1,(1+(j-1)*l1):j*l1] = ada(b,ind,i,j)
            end
        end
    end
    op_general = z+z'
    return op_general, l1, b
end

#operator c^dagger c for a system with fixed particle number.
# Much faster than using fiexed over the whole operators

function ada(n::Int64, m::Int64, i::Int64, j::Int64)
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
function ada(base::SparseArrays.SparseMatrixCSC{Float64,Int64}, ind::SparseArrays.SparseVector{Float64,Int64}, i::Int64, j::Int64)
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
function aad(n::Int64, m::Int64, i::Int64, j::Int64)
    return -ada(n,m,j,i)
end

function aad(base::SparseArrays.SparseMatrixCSC{Float64,Int64}, indice::SparseArrays.SparseVector{Float64,Int64}, i::Int64, j::Int64)
    return -ada(base,indice,j,i)
end

# This is the creation operator adagger
# it serves as ancilla creation or in the same
# sp space. m is the number of particles in The
# starting state
function ad(n::Int64, m::Int64, i::Int64)
    base1, ind1 = basis_m(n,m)
    l2 = binomial(n,m)
    # if I am in the original sp space
    if i<=n
        base2, ind2 = basis_m(n,m+1)
        l1 = binomial(n,m+1)
        op = spzeros(l1,l2)
        for k in 1:l2
            if base1[k,:][i] == 0
                j = myfind(ind2,(ind1[k]+2^(n-i)))[1]
                sign = (-1)^(sum(base1[k,1:(i-1)]))
                op[j,k] = sign
            end
        end
    # else, create an ancilla
    else
        base2, ind2 = basis_m(i,m+1)
        l1 = binomial(i,m+1)
        op = spzeros(l1,l2)
        for k in 1:l2
            x = join(Int.(append!(Array(base1[k,:]),1)))
            x = parse(Int, x; base=2) + 1
            j = myfind(ind2,x)[1]
            sign = (-1)^m
            op[j,k] = sign
        end
    end
    return op
end

# This is the destruction operator.
# m is the number of particles in The
# starting state. It can also be built from ad
# computing the transpose
function a(n::Int64, m::Int64, i::Int64)
    if i <= n
        base1, ind1 = basis_m(n,m)
        base2, ind2 = basis_m(n,m-1)
        l1 = binomial(n,m-1)
        l2 = binomial(n,m)
        op = spzeros(l1,l2)
        for k in 1:l2
            if base1[k,:][i] != 0
                j = myfind(ind2,(ind1[k]-2^(n-i)))[1]
                sign = (-1)^(sum(base1[k,1:(i-1)]))
                op[j,k] = sign
            end
        end
        return op
    else
        return spzeros(binomial(n,m))
    end
end

# This serves for applying several creations op
# at once. Modes will be a vector with the chosen
# modes
function ad(n::Int64, m::Int64, modes::Array{Int64,1})
    l = length(modes)
    op  = ad(n,m,modes[1])
    for k in 2:l
        op = ad(n,m+k-1,modes[k])*op
    end
    return op
end

# This serves for applying several destruction op
# at once. Modes will be a vector with the chosen
# modes
function a(n::Int64, m::Int64, modes::Array{Int64,1})
    l = length(modes)
    op  = a(n,m,modes[1])
    for k in 2:l
        op = a(n,m-k+1,modes[k])*op
    end
    return op
end

# The following function initializes the creation
# operators for each mode given a fixed dimension d
# and the starting number of particles m
function adf(d::Int, m::Int)
    adtot = [ad(d,m,i) for i in 1:d]
    return adtot
end

# The following function initializes the destruction
# operators for each mode given a fixed dimension d
# and the starting number of particles m
function af(d::Int, m::Int)
    atot = [a(d,m,i) for i in 1:d]
    return atot
end

#This functions trimms the operators matrix
#for fixed particle number subspace
function fixed(o, nume::Int64)
    n = Int(log(2,size(o)[1]))
    posn = ffilter(n, nume)
    reduced = o[setdiff(1:end,posn), setdiff(1:end,posn)]
    return reduced
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


## What follows is the op_semifixed struct
# that I previously used for initalizing all
# creation and destruction op in fixed
# particle number subspaces, switching
# from one to the other.
# I prefer now initalizing only adf and af
# as it is clearer and spanns most cases

# You should add LazyStack to the dependancies:
#LazyStack = "1fad7336-0346-5a1a-a56f-a06ba010965b"
# in Fermionic.jl:
# using LazyStack

# # The following struct initialized the creation
# # and destruction operators, given a diemension
# struct Op_semifixed
#     d::Int
#     adtot::Array{Array{Float64,3},1}
#     Op_semifixed(d) = new(d, opnfixed(d))
# end
#
# ad(o::Op_semifixed, m::Int, i::Int) = sparse(o.adtot[m][:,:,i])
# # a(o::Op_semifixed, m::Int, i::Int) = sparse(o.adtot[m-1][:,:,i]')
#
# function a(o::Op_semifixed, m::Int, i::Int)
#     if i <= o.d
#         return sparse(o.adtot[m-1][:,:,i]')
#     else
#         return spzeros(binomial(o.d,m))
#     end
# end

# # This function is the core of Op_semifixed
# # initialization. It uses rstack from the package
# # Lazystack. It looses the sparse structure in doing so
# # A pending task is to improve this mechanism
# # We only initialize ad, and compute a as the transpose
# # Changing this function would allow deleting
# # a package from the dependancies of the package
# function opnfixed(d::Int)
#     r = [rstack([ad(d,j,i) for i in 1:d];fill=0.) for i in 1:d for j in 1:(d-1)]
#     return r
# end

# # This serves for applying several creations op
# # at once. Modes will be a vector with the chosen
# # modes
# function ad(o::Op_semifixed, m::Int64, modes::Array{Int64,1})
#     l = length(modes)
#     op  = ad(o,m,modes[1])
#     for k in 2:l
#         op = ad(o,m+k-1,modes[k])*op
#     end
#     return op
# end

# # This serves for applying several destruction op
# # at once. Modes will be a vector with the chosen
# # modes
# function a(o::Op_semifixed, m::Int64, modes::Array{Int64,1})
#     l = length(modes)
#     op  = a(o,m,modes[1])
#     for k in 2:l
#         op = a(o,m-k+1,modes[k])*op
#     end
#     return op
# end

# The documentation:
# - **Op_semifixed(d::Int)**: initializes the fermionic operators $a_i$ and $a_i^\dagger$ corresponding to dimension **d** in the fixed particle subspace. Once initizialized, the following functions are made avaiable.
#     - **a(o::Op_semifixed, n::Int, i::Int)**: returns the fermionic destruction operator in mode i, i.e. $a_i$, corresponding to the previously defined **o** in the fixed particle subspace with $n$ starting particles.
#     - **a(o::Op_semifixed,  n::Int, mod::Array{Int,1})**: returns the consecutive fermionic destruction operator of modes in vecor **mod**, i.e. $a_{mod_m}...a_{mod_1}$, corresponding to the previously defined **o** in the fixed particle subspace with $n$ starting particles.
#     - **ad(o::Op_semifixed,  n::Int, i::Int)**: returns the fermionic creation operator in mode **i**, i.e. $a^\dagger_i$, corresponding to the previously defined **o** in the fixed particle subspace with $n$ starting particles.
#      - **ad(o::Op_semifixed,  n::Int, mod::Array{Int,1})**: returns the consecutive fermionic creation operator of modes in vecor **mod**, i.e. $a^\dagger_{mod_m}...a^\dagger_{mod_1}$, corresponding to the previously defined **o** in the fixed particle subspace with $n$ starting particles.
