#core.jl

function operators(n)
    rowb, colb, base, vacio = integer_digits(n)
    l = 2^n
    lb = n*2^(n-1)
    row = spzeros(lb+1)
    col = spzeros(lb+1)
    data = spzeros(lb+1)
    row[lb+1] = l
    col[lb+1] = 1
    data[lb+1] = 0
    for i in 1:lb
        j = floor(Int, rowb[i])-2^(n-floor(Int, colb[i]))
        sign = (-1)^(sum(base[floor(Int, rowb[i]),1:floor(Int, colb[i])])+1)
        row[i] = j
        col[i] = floor(Int, rowb[i])+(floor(Int, colb[i])-1)*l
        data[i] = sign
    end
    cm_tot = sparse(row, col, data)
    cd_tot = sparse(cm_tot')
    return cm_tot, cd_tot, l, base
end

#This function outputs elements for building
#the sparse matrix of basis
function integer_digits(n)
    rowb = spzeros(n*2^(n-1))
    colb = spzeros(n*2^(n-1))
    data = spzeros(n*2^(n-1))
    counter = 1
    for i in 0:(2^n-1)
        binary_base = bitstring(i)[65-n:64]
        bin_vector = splitter(binary_base)
        for j in 1:n
            if bin_vector[j] == 1
                rowb[counter] = (i+1)::Int
                colb[counter] = j::Int
                data[counter] = 1
                counter = counter + 1
            end
        end
    end
    base = sparse(rowb, colb, data)
    vacio = spzeros(2^n)
    vacio[1] = 1
    return rowb, colb, base, vacio
end


function splitter(x)
    vectorized = []
    for i in x
        str_to_int = parse(Int, i)
        append!(vectorized, str_to_int)
    end
    return vectorized
end
