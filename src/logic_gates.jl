
#U control not
function ucnot(o::Op, control::Int64, target::Int64) #c is control, t is target

    c = control
    t = target
    base = basis(o)
    d = dim(o)

    #we first check that control is different from target:
    if c == t
        throw(ArgumentError("Control must be different from target"))
    end

    if c > d || t > d
        throw(ArgumentError("Your modes must be within your dimensions!"))
    end

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
    u_cnot = sparse(row, col, data)
    return u_cnot
end

#Hadamard
function hadamard(o::Op, mode1::Int64, mode2::Int64)
    #we first check that the input modes are different
    if mode1 == mode2
        throw(ArgumentError("Modes must be differents"))
    end

    i = min(mode1, mode2)
    j = max(mode1, mode2)
    base = basis(o)
    d = dim(o)
    l = 2^d

    if mode1 > d || mode2 > d
        throw(ArgumentError("Your modes must be within your dimensions!"))
    end

    row = spzeros(l+2^(d-1))
    col = spzeros(l+2^(d-1))
    data = spzeros(l+2^(d-1))
    counter = 1
    for k in 1:l
        if base[k,i] != 0
            if base[k,j] != 0
                row[counter] = k
                col[counter] = k
                data[counter] = 1.0
                counter = counter + 1
            else
                row[counter] = k
                col[counter] = k
                data[counter] = 1.0/sqrt(2)
                counter = counter + 1
                row[counter] = k
                col[counter] = k - 2^(d-i) + 2^(d-j)
                data[counter] = 1.0/sqrt(2)
                counter = counter + 1
            end

        elseif base[k,j] != 0
            row[counter] = k
            col[counter] = k
            data[counter] = 1.0/sqrt(2)
            counter = counter + 1
            row[counter] = k
            col[counter] = k + 2^(d-i) - 2^(d-j)
            data[counter] = 1.0/sqrt(2)
            counter = counter + 1
        else
            row[counter] = k
            col[counter] = k
            data[counter] = 1.0
            counter = counter + 1
        end
    end
    had = sparse(row,col,data)
    return had
end
