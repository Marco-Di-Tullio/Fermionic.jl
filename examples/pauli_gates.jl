using Fermionic
using SparseArrays
o = Op(3)

sigma_x = cdcm(o,1,2) + cdcm(o,2,1)

sigma_x2 = cdcd(o,1,2) + cmcm(o,1,2)

sigma_x3 = cdcm(o,1,2) + cdcm(o,2,1) + cdcd(o,1,2) + cmcm(o,1,2)

#esta es la negación de 2 modos
#Quiero hacer la negación de 1 modo
o = Op(3)

function not(o, mode) #c is control, t is target

    base = basis(o)
    d = dim(o)
    l = 2^d

    if mode > d
        throw(ArgumentError("Your modes must be within your dimensions!"))
    end

    row = spzeros(l)
    col = spzeros(l)
    data = spzeros(l)
    for k in 1:l
        row[k] = k
        if base[k,mode] != 0
            col[k] = k - 2^(d-mode)
            data[k] = 1.0
        else
            col[k] = k + 2^(d-mode)
            data[k] = 1.0
        end
    end
    not = sparse(row, col, data)
    return not
end

#Swap:
function swap(o, mode1, mode2) #c is control, t is target

    base = basis(o)
    d = dim(o)
    l = 2^d

    if mode1 > d || mode2 > d
        throw(ArgumentError("Your modes must be within your dimensions!"))
    end

    if mode1 == mode2
        throw(ArgumentError("Control must be different from target"))
    end

    row = spzeros(l)
    col = spzeros(l)
    data = spzeros(l)
    for k in 1:l
        row[k] = k
        if base[k, mode1] != 0
            if base[k, mode2] != 0
                col[k] = k
                data[k] = 1.0
            else
                col[k] = k - 2^(d-mode1) + 2^(d-mode2)
                data[k] = 1.0
            end
        else
            if base[k, mode2] != 0
                col[k] = k + 2^(d-mode1) - 2^(d-mode2)
                data[k] = 1.0
            else
                col[k] = k
                data[k] = 1.0
            end
        end
    end
    swap = sparse(row, col, data)
    return swap
end
