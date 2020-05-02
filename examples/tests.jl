using Fermionic
using SparseArrays

o = Op(2)
mode = 1
function sigma_y(o, mode) #c is control, t is target

    base = basis(o)
    d = dim(o)
    l = 2^d

    if mode > d
        throw(ArgumentError("Your mode must be within your dimensions!"))
    end

    row = spzeros(l)
    col = spzeros(l)
    data = spzeros(Complex{Float64},l)
    for k in 1:l
        row[k] = k
        if base[k,mode] != 0
            col[k] = k - 2^(d-mode)
            data[k] = 1.0*im
        else
            col[k] = k + 2^(d-mode)
            data[k] = -1.0*im
        end
    end
    sigma_y = sparse(row, col, data)
    return sigma_y
end

function sigma_z(o, mode) #c is control, t is target

    base = basis(o)
    d = dim(o)
    l = 2^d

    if mode > d
        throw(ArgumentError("Your mode must be within your dimensions!"))
    end

    row = spzeros(l)
    col = spzeros(l)
    data = spzeros(l)
    for k in 1:l
        row[k] = k
        if base[k,mode] != 0
            col[k] = k
            data[k] = -1.0
        else
            col[k] = k
            data[k] = 1.0
        end
    end
    sigma_z = sparse(row, col, data)
    return sigma_z
end
