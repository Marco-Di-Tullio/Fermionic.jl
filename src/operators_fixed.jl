#This functions trimms the operators matrix
#for fixed particle number subspace
function fixed(o, nume)
    n = Int(log(2,size(o)[1]))
    posn = ffilter(n, nume)
    reduced = o[setdiff(1:end,posn), setdiff(1:end,posn)]
    return reduced
end

function fixed_state(stat, nume)
    n = Int(log(2,length(stat)))
    posn = ffilter(n, nume)
    reduced = stat[setdiff(1:end,posn)]
    return reduced
end


function ffilter(n, m)
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

function myfind(c,j)
    a = similar(c, Int)
    count = 1
    @inbounds for i in eachindex(c)
        a[count] = i
        count += (c[i] == j)
    end
    return a[1:(count-1)]
end
