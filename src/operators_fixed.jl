#This functions trimms the operators matrix
#for fixed particle number subspace
function fixed(o, n, nume)
    posn = filter(n, nume)
    reduced = o[setdiff(1:end,posn), setdiff(1:end,posn)]
    return reduced
end

function filter(n, m)
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


function basis_m(o::Op, m::Int)
    d = dim(o)
    basm = spzeros(binomial(d,m),d)
    counter = 1
    baso = basis(o)
    for i in 1:2^d
        if sum(baso[i,:]) == m
            basm[counter,:] = baso[i,:]
            counter = counter + 1
        end
    end
    return basm
end
