prec = 15 #precision

##Single Particle
eigensp(s::State) = [round(eigvals(rhosp(s))[i], digits = prec) for i in 1:dim(ope(s))]
eigensp(s::State_complex) = [round(eigvals(rhosp(s))[i], digits = prec) for i in 1:dim(ope(s))]
#For some reason, I can not compute eigenvalues
#directly from sparse matrices
eigensp(s::State_sparse) = [round(eigvals(Matrix(rhosp(s)))[i], digits = prec) for i in 1:dim(ope(s))]
eigensp(s::State_sparse_complex) = [round(eigvals(Matrix(rhosp(s)))[i], digits = prec) for i in 1:dim(ope(s))]

function ssp(sta::State)
    eigen = eigensp(sta)
    lene = length(eigen)
    s = 0
    for i in 1:lene
        if eigen[i] != 0 && eigen[i] != 1
            s = s - (eigen[i]*log(2,eigen[i]) + (1 - eigen[i])*log(2,1-eigen[i]))
        end
    end
    return s/lene
end

function ssp(sta::State_complex)
    eigen = eigensp(sta)
    lene = length(eigen)
    s = 0
    for i in 1:lene
        if eigen[i] != 0 && eigen[i] != 1
            s = s - (eigen[i]*log(2,eigen[i]) + (1 - eigen[i])*log(2,1-eigen[i]))
        end
    end
    return s/lene
end

function ssp(sta::State_sparse)
    eigen = eigensp(sta)
    lene = length(eigen)
    s = 0
    for i in 1:lene
        if eigen[i] != 0 && eigen[i] != 1
            s = s - (eigen[i]*log(2, eigen[i]) + (1 - eigen[i])*log(2, 1-eigen[i]))
        end
    end
    return s/lene
end

function ssp(sta::State_sparse_complex)
    eigen = eigensp(sta)
    lene = length(eigen)
    s = 0
    for i in 1:lene
        if eigen[i] != 0 && eigen[i] != 1
            s = s - (eigen[i]*log(2, eigen[i]) + (1 - eigen[i])*log(2, 1-eigen[i]))
        end
    end
    return s/lene
end

## Quasi Particles
eigenqsp(s::State) = [round(eigvals(rhoqsp(s))[i], digits = prec) for i in 1:(2*dim(ope(s)))]
eigenqsp(s::State_complex) = [round(eigvals(rhoqsp(s))[i], digits = prec) for i in 1:(2*dim(ope(s)))]
#For some reason, I can not compute eigenvalues
#directly from sparse matrices
eigenqsp(s::State_sparse) = [round(eigvals(Matrix(rhoqsp(s)))[i], digits = prec) for i in 1:(2*dim(ope(s)))]
eigenqsp(s::State_sparse_complex) = [round(eigvals(Matrix(rhoqsp(s)))[i], digits = prec) for i in 1:(2*dim(ope(s)))]


function sqsp(sta::State)
    eigen = eigenqsp(sta)
    lene = length(eigen)
    s = 0
    for i in 1:lene
        if eigen[i] != 0 && eigen[i] != 1
            s = s - (eigen[i]*log(2,eigen[i]) + (1 - eigen[i])*log(2,1-eigen[i]))
        end
    end
    return s/lene
end

function sqsp(sta::State_complex)
    eigen = eigenqsp(sta)
    lene = length(eigen)
    s = 0
    for i in 1:lene
        if eigen[i] != 0 && eigen[i] != 1
            s = s - (eigen[i]*log(2,eigen[i]) + (1 - eigen[i])*log(2,1-eigen[i]))
        end
    end
    return s/lene
end

function sqsp(sta::State_sparse)
    eigen = eigenqsp(sta)
    lene = length(eigen)
    s = 0
    for i in 1:lene
        if eigen[i] != 0 && eigen[i] != 1
            s = s - (eigen[i]*log(2, eigen[i]) + (1 - eigen[i])*log(2, 1-eigen[i]))
        end
    end
    return s/lene
end

function sqsp(sta::State_sparse_complex)
    eigen = eigenqsp(sta)
    lene = length(eigen)
    s = 0
    for i in 1:lene
        if eigen[i] != 0 && eigen[i] != 1
            s = s - (eigen[i]*log(2, eigen[i]) + (1 - eigen[i])*log(2, 1-eigen[i]))
        end
    end
    return s/lene
end

#This macro serves for printing as a string the variable's name
macro Name(arg)
   string(arg)
end

#This function outputs 1 if s1 majorizes s2,
#-1 if s2 majorizes s1,
#0 if there is no majorization relation
function majorization_sp(s1, s2)
    e1 = sort(eigensp(s1), rev = true)
    e2 = sort(eigensp(s2), rev = true)
    count1 = 0
    count2 = 0
    d = length(e1)
    for j in 1:(d-1) #el último valor a veces hace problemas
        if sum([round(e1[i]-e2[i], digits = 10) for i in 1:j]) < 0
            count1 = count1 + 1
        elseif sum([round(e2[i]-e1[i], digits = 10) for i in 1:j]) < 0
            count2 = count2 + 1
        end
    end
    if count1 == (d-1)
        println(@Name(s2)," sp majorizes ", @Name(s1))
        return -1
    elseif count2 == (d-1)
        println(@Name(s1)," sp majorizes ", @Name(s2))
        return 1
    else
        println("There is no sp majorization relation")
        return 0
    end
end


function majorization_qsp(s1, s2)
    e1 = sort(eigenqsp(s1), rev = true)
    e2 = sort(eigenqsp(s2), rev = true)
    count1 = 0
    count2 = 0
    d = length(e1)
    for j in 1:(d-1) #el último valor a veces hace problemas
        if sum([round(e1[i]-e2[i], digits = 10) for i in 1:j]) < 0
            count1 = count1 + 1
        elseif sum([round(e2[i]-e1[i], digits = 10) for i in 1:j]) < 0
            count2 = count2 + 1
        end
    end
    if count1 == (d-1)
        println(@Name(s2)," qsp majorizes ", @Name(s1))
        return -1
    elseif count2 == (d-1)
        println(@Name(s1)," qsp majorizes ", @Name(s2))
        return 1
    else
        println("There is no qsp majorization relation")
        return 0
    end
end
