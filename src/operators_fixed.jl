#This functions trimms the operators matrix
#for fixed particle number subspace
function fixed(o, nume)
    n = Int(log(2,size(o)[1]))
    posn = filter(n, nume)
    reduced = o[setdiff(1:end,posn), setdiff(1:end,posn)]
    return reduced
end

function fixed_state(stat, nume)
    n = Int(log(2,length(stat)))
    posn = filter(n, nume)
    reduced = stat[setdiff(1:end,posn)]
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
