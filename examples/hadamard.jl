using Fermionic
using SparseArrays


#U control not
function ucnot(o::Op, c::Int64, t::Int64) #c is control, t is target

    #we first check that control is different from target:
    if c == t
        throw(ArgumentError("Control must be different from target"))
    end

    base = basis(o)
    d = dim(o)
    l = 2^d
    row = spzeros(l)
    col = spzeros(l)
    data = spzeros(l)
    for k in 1:l
        row[k] = k
        if base[k,c] != 0
            if base[k,t] == 1.0
                col[k] = k - 2^(d-t)
                data[k] = 1.0
            else
                col[k] = k + 2^(d-t)
                data[k] = 1.0
            end
        else
            col[k] = k
            data[k] = 1.0
        end
    end
    mat = sparse(row, col, data)
    return mat
end


function hadamard(o, mode1, mode2)
    #we first check that the input modes are different
    if mode1 == mode2
        throw(ArgumentError("Modes must be differents"))
    end

    i = min(mode1, mode2)
    j = max(mode1, mode2)
    base = basis(o)
    d = dim(o)
    l = 2^d
    row = []
    col = []
    data = []
    for k in 1:l
        if base[k,i] != 0
            if base[k,j] != 0
                append!(row, k)
                append!(col, k)
                append!(data, 1.0)
            else
                append!(row, [k, k])
                append!(col, [k, k - 2^(d-i) + 2^(d-j)])
                append!(data, [1/sqrt(2), 1/sqrt(2)])
            end
        elseif base[k,j] != 0
            append!(row, [k, k])
            append!(col, [k, k + 2^(d-i) - 2^(d-j)])
            append!(data, [1/sqrt(2), 1/sqrt(2)])
        else
            append!(row, k)
            append!(col, k)
            append!(data, 1.0)
        end
    end
    mat = sparse(row,col,data)
    return mat
end
