
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
