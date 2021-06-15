#prec = 15 #precision

##Single Particle
eigensp(s::State, prec=15) = sort([round(eigvals(Matrix(rhosp(s)))[i], digits = prec) for i in 1:dim(ope(s))], rev=true)
eigensp(s::State_fixed, prec=15) = sort([round(eigvals(Matrix(rhosp(s)))[i], digits = prec) for i in 1:dim(ope(s))], rev=true)
# I can not compute eigenvalues
#directly from sparse matrices

function ssp(sta::State, prec=15)
    eigen = eigensp(sta, prec)
    lene = length(eigen)
    s = 0
    for i in 1:lene
        if eigen[i] != 0 && eigen[i] != 1
            s = s - (eigen[i]*log(2,eigen[i]) + (1 - eigen[i])*log(2,1-eigen[i]))
        end
    end
    return s/lene
end

function ssp(sta::State_fixed, prec=15)
    eigen = eigensp(sta, prec)
    lene = length(eigen)
    s = 0
    for i in 1:lene
        if eigen[i] != 0 && eigen[i] != 1
            s = s - (eigen[i]*log(2,eigen[i]) + (1 - eigen[i])*log(2,1-eigen[i]))
        end
    end
    return s/lene
end
## Quasi Particles
eigenqsp(s::State, prec=15) = sort([round(eigvals(Matrix(rhoqsp(s)))[i], digits = prec) for i in 1:(2*dim(ope(s)))], rev=true)
# I can not compute eigenvalues
#directly from sparse matrices

function sqsp(sta::State, prec=15)
    eigen = eigenqsp(sta, prec)
    lene = length(eigen)
    s = 0
    for i in 1:lene
        if eigen[i] != 0 && eigen[i] != 1
            s = s - (eigen[i]*log(2,eigen[i]) + (1 - eigen[i])*log(2,1-eigen[i]))
        end
    end
    return s/lene
end

#This function outputs 1 if s1 majorizes s2,
#-1 if s2 majorizes s1,
#0 if there is no majorization relation
function majorization_sp(s1, s2)
    e1 = eigensp(s1)
    e2 = eigensp(s2)
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
        #println(@Name(s2)," sp majorizes ", @Name(s1))
        return -1
    elseif count2 == (d-1)
        #println(@Name(s1)," sp majorizes ", @Name(s2))
        return 1
    else
        #println("There is no sp majorization relation")
        return 0
    end
end

function majorization_qsp(s1, s2)
    e1 = eigenqsp(s1)
    e2 = eigenqsp(s2)
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
        #println(@Name(s2)," qsp majorizes ", @Name(s1))
        return -1
    elseif count2 == (d-1)
        #println(@Name(s1)," qsp majorizes ", @Name(s2))
        return 1
    else
        #println("There is no qsp majorization relation")
        return 0
    end
end

function n_avg(s::State)
    navg = Real(tr(rhosp(s)))
    return navg
end

####################################################
#The following are the rhom matrices without diagonalization
# i.e. in the original basis
function rhom(s::State, m::Int64)
    o = ope(s)
    num = Int(round(n_avg(s),digits = 8))
    d = dim(o) #tengo que pasar a matrix lamentablemente
    rhomd = spzeros(typ(s),binomial(d, m), binomial(d, m))
    bas,_ = basis_m(d, m)
    estat = st(s)
    and = sparse([non_diag_ops(o, d, num, m, bas, j) for j in 1:binomial(d, m)])
    for i in 1:binomial(d, m)
        for j in 1:binomial(d,m)
            rhomd[i,j] = round(estat'*and[j]*and[i]'*estat, digits = 14)
        end
    end
    return rhomd
end

# método directo (muy lento)
# function rhom2(s::State)
#     d = dim(ope(s))
#     o = ope(s)
#     bas,_ = basis_m(d,2)
#     l = Int(binomial(d,2))
#     rhom = spzeros(l,l);
#     stat = st(s)
#     for i in 1:l
#         indi = indx(bas[i,:])
#         for j in 1:l
#             indj = indx(bas[j,:])
#                 rhom[i,j]=(stat'*ad(o,Int(indi[2]))*ad(o,Int(indi[1]))*a(o,Int(indj[1]))*a(o,Int(indj[2]))*stat)
#         end
#     end
#     return rhom
# end

#rhom(State(1/sqrt(2)*(ad(o,1)*ad(o,2)*ad(o,3)*ad(o,4)+ad(o,3)*ad(o,4)*ad(o,5)*ad(o,6))*vacuum(o),o),2)

# It works, but it is not efficient,
# we need to optimieze ir for state_fixed
function rhom(s::State_fixed, o::Op, m::Int64)
    d = dim(o)
    #o = Op(d)
    num = nume(s)
    rhomd = spzeros(typ(s),binomial(d, m), binomial(d, m))
    bas,_ = basis_m(d, m)
    estat = st(s)
    and = sparse([non_diag_ops(o, d, num, m, bas, j) for j in 1:binomial(d, m)])
    for i in 1:binomial(d, m)
        for j in 1:binomial(d,m)
            rhomd[i,j] = round(estat'*fixed(and[j]*and[i]',num)*estat, digits = 14)
        end
    end
    return rhomd
end

#main function
function rhomd(s::State, m::Int64)
    o = ope(s)
    num = Int(round(n_avg(s),digits = 8))
    d = dim(o) #tengo que pasar a matrix lamentablemente
    rhomd = spzeros(typ(s),binomial(d, m), binomial(d, m))
    bas, _ = basis_m(d, m)
    estat = st(s)
    ty = typ(s)
    u, _ = svd(cof_mat(o, d, num, m, bas, estat, ty))
    and = sparse([non_diag_ops(o, d, num, m, bas, j) for j in 1:binomial(d, m)])
    adi = sparse([diag_ops(o, d, num, m, bas, u, and, j) for j in 1:binomial(d, m)])
    #adi[1] es evaluar la funcion diagonal en i=1
    for i in 1:binomial(d, m)
        rhomd[i,i] = round(estat'*adi[i]*adi[i]'*estat, digits = 14)
    end
    return rhomd
end

# It works, but it is not efficient,
# we need to optimieze ir for state_fixed
function rhomd(s::State_fixed, o::Op, m::Int64)
    d = dim(o)
    #o = Op(d)
    num = nume(s)
    rhomd = spzeros(typ(s),binomial(d, m), binomial(d, m))
    bas_tot, _ = basis_m(d, num)
    bas, _ = basis_m(d, m)
    estat = st(s)
    ty = typ(s)
    u, _ = svd(cof_mat_fixed(o, d, num, m, bas, bas_tot, estat, ty))
    and = sparse([non_diag_ops(o, d, num, m, bas, i) for i in 1:binomial(d, m)])
    adi = sparse([diag_ops(o, d, num, m, bas, u, and, i) for i in 1:binomial(d, m)])
    #ad[1] es evaluar la funcion diagonal en i=1
    for i in 1:binomial(d, m)
        rhomd[i,i] = round(estat'*fixed(adi[i]*adi[i]',num)*estat, digits = 14)
    end
    return rhomd
end

#A more efficient 2 body DM in fixed space
function rhom2(s::State_fixed, prec=15)
    m = 2
    o = ope(s)
    num = nume(s)
    d = dim(o)
    rhomnd = spzeros(typ(s),binomial(d, m), binomial(d, m))
    bas,_ = basis_m(d,m)
    estat = st(s)
    for i in 1:binomial(d, m)
        indi = indx(bas[i,:])
        for j in 1:binomial(d,m)
            indj = indx(bas[j,:])
            #rhomd[i,j] = round(estat'*fixed(and[j]*and[i]',num)*estat, digits = 14)
            if indi[2] != indj[1] && indi[1] != indj[2]
                rhomnd[i,j] = round(estat'*ada(o,indj[1],indi[1])*ada(o,indj[2],indi[2])*estat, digits = prec)
            else
                rhomnd[i,j] = round(-estat'*ada(o,indj[1],indi[2])*ada(o,indj[2],indi[1])*estat, digits = prec)
            end
        end
    end
    return rhomnd
end

#the following are m-body matrices: work on progress
#for example, for a vector [0,1,0,0,1]->[2,5]
function indx(arr)
    l = length(arr)
    ind = spzeros(Int,0)
    for i=1:l
        if arr[i] != 0
            ind = sparse([ind; i])
        end
    end
    return ind
end

function cof_mat(o, d, num, m, bas, estat, ty)
    cmat = zeros(ty, binomial(d,m),binomial(d,num-m))
    am = sparse([bas[i,:] for i in 1:binomial(d, m)])
    bas_nume_m, _ = basis_m(d, num - m)
    amn = sparse([bas_nume_m[i,:] for i in 1:binomial(d, num-m)])
    for i in 1:binomial(d, m)
        ai = am[i]
        indi = indx(ai)
        for j in 1:binomial(d, num-m)
            aj = amn[j]
            if ai'*aj != 0
                cmat[i,j] = 0
            else
                permut = 0
                for k in 1:length(indi)
                    permut = permut + sum(aj[1:floor(Int,indi[k]-1)])
                end
                indij = indx(ai+aj)
                elem = Int(sum([2^(d-i) for i in indij]) + 1)
                cmat[i,j] = estat[elem] * (-1)^permut
            end
        end
    end
    return cmat
end

function cof_mat_fixed(o, d, num, m, bas, bas_tot, estat, ty)
    cmat = zeros(ty,binomial(d,m),binomial(d,num-m))
    am = sparse([bas[i,:] for i in 1:binomial(d, m)])
    bas_nume_m, _ = basis_m(d, num - m)
    amn = sparse([bas_nume_m[i,:] for i in 1:binomial(d, num-m)])
    for i in 1:binomial(d, m)
        ai = am[i]
        indi = indx(ai)
        for j in 1:binomial(d, num-m)
            aj = amn[j]
            if ai'*aj != 0
                cmat[i,j] = 0
            else
                permut = 0
                for k in 1:length(indi)
                    permut = permut + sum(aj[1:floor(Int,indi[k]-1)])
                end
                indij = indx(ai+aj)
                elem = 1
                while indx(bas_tot[elem,:]) != indij
                    elem = elem + 1
                end
                cmat[i,j] = estat[elem] * (-1)^permut
            end
        end
    end
    return cmat
end

function non_diag_ops(o, d, nume, m, bas, i)
    vec = indx(bas[i, :])
    #le es nume-1 me parece
    mel = 1
    for j in 1:m
        mel = mel*ad(o, Int(vec[j]))
    end
    return mel
end

function diag_ops(o, d, nume, m, bas, u, and, i)
    adi = spzeros(2^d, 2^d)
    for j in 1:binomial(d,m)
        if round(u[j,i], digits=14) != 0
            adi = adi + u[j, i]*and[j]
        end
    end
    return adi
end

#= Work in Progress: non diagonal operators with fixed states
function rhomnd(s::Union{State_fixed,State_sparse_fixed}, m::Int64)
    o = ope(s)
    num = nume(s)
    d = dim(o) #tengo que pasar a matrix lamentablemente
    rhomd = spzeros(binomial(d, m), binomial(d, m))
    bas, _ = basis_m(d, m)
    estat = st(s)
    and = sparse([non_diag_ops(o, d, num, m, bas, j) for j in 1:binomial(d, m)])
    for i in 1:binomial(d, m)
        for j in 1:binomial(d,m)
            rhomd[i,j] = round(estat'*and[j]*and[i]'*estat, digits = 14)
        end
    end
    return rhomd
end

#main function
function rhomnd(s::Union{State_complex_fixed,State_sparse_complex_fixed}, m::Int64)
    o = ope(s)
    num = nume(s)
    d = dim(o) #tengo que pasar a matrix lamentablemente
    rhomd = spzeros(Complex{Float64},binomial(d, m), binomial(d, m))
    bas, _ = basis_m(d, m)
    estat = st(s)
    and = sparse([non_diag_ops(o, d, num, m, bas, j) for j in 1:binomial(d, m)])
    for i in 1:binomial(d, m)
        for j in 1:binomial(d,m)
            rhomd[i,j] = round(estat'*and[j]*and[i]'*estat, digits = 14)
        end
    end
    return rhomd
end
=#

#Partial trace in a given basis of modes
function trp(state::State, modos::Array{Int64,1})
    d = dim(ope(state))
    bas = basis(ope(state))
    sta = st(state)
    lm = length(modos)
    zvr = zeros(typ(state),2^lm,2^(d-lm))
    full = [i for i in 1:d]
    lista = sort(modos)
    listb = filter(x->x ∉ lista,full)
    for k in 1:2^d
        indicea = parse(Int,join([Int(bas[k,i]) for i in lista]), base=2)+1
        indiceb = parse(Int,join([Int(bas[k,i]) for i in listb]), base=2)+1
        signi = 0
        for i in lista
            listac = filter(x->x <= i-1,listb)
            signl = 0
            for l in listac
                signl = signl+bas[k,l]
            end
            signi = signi + signl*bas[k,i]
        end
        sign = (-1)^(signi)
        zvr[indicea,indiceb] = sta[k]*sign
    end
    rhoa = zvr*zvr'
    #rhob=zvr'*zvr;
    return rhoa
end

#= Work in Progress: partial trace for fixed particle number
function trp(state::Union{State_fixed,State_sparse_fixed},modos::Array{Int64,1})
    d = dim(ope(state))
    bas = basis(ope(state))
    sta = st(state)
    lm = length(modos)
    zvr = zeros(Complex{Float64},2^lm,2^(d-lm))
    full = [i for i in 1:d]
    lista = sort(modos)
    listb = filter(x->x ∉ lista,full)
    for k in 1:2^d
        indicea=parse(Int,join([Int(bas[k,i]) for i in lista]), base=2)+1
        indiceb=parse(Int,join([Int(bas[k,i]) for i in listb]), base=2)+1
        signi = 0
        for i in lista
            listac = filter(x->x <= i-1,listb)
            signl = 0
            for l in listac
                signl = signl+bas[k,l]
            end
            signi = signi + signl*bas[k,i]
        end
        sign = (-1)^(signi)
        zvr[indicea,indiceb]=sta[k]*sign
    end
    rhoa=zvr*conj(zvr')
    #rhob=zvr'*zvr;
    return rhoa
end

function trp(state::Union{State_complex_fixed,State_sparse_complex_fixed},modos::Array{Int64,1})
    d = dim(ope(state))
    bas = basis(ope(state))
    sta = st(state)
    lm = length(modos)
    zvr = zeros(Complex{Float64},2^lm,2^(d-lm))
    full = [i for i in 1:d]
    lista = sort(modos)
    listb = filter(x->x ∉ lista,full)
    for k in 1:2^d
        indicea=parse(Int,join([Int(bas[k,i]) for i in lista]), base=2)+1
        indiceb=parse(Int,join([Int(bas[k,i]) for i in listb]), base=2)+1
        signi = 0
        for i in lista
            listac = filter(x->x <= i-1,listb)
            signl = 0
            for l in listac
                signl = signl+bas[k,l]
            end
            signi = signi + signl*bas[k,i]
        end
        sign = (-1)^(signi)
        zvr[indicea,indiceb]=sta[k]*sign
    end
    rhoa=zvr*zvr'
    #rhob=zvr'*zvr;
    return rhoa
end
=#

#= This code is faster for many iterations (majorization)
count = 0
@time for i in 1:repeticiones
    #definimos estado premedida
    state_ran = spzeros(Complex{Float64},l)
    for i in 1:l
        if sum(basis(o)[i,:]) == nume
            state_ran[i] = 2*rand(Complex{Float64},1)[1]-1-im
        end
    end
    state_ran = state_ran/sqrt(state_ran'*state_ran)
    state_ran = State_sparse_complex(state_ran, o);
    #e = sort(eigvals(Matrix(rhoqsp(state_ran))), rev=true)

    #ahora medimos
    state_ran_medido = operator*st(state_ran);
    state_ran_medido = state_ran_medido/sqrt(state_ran_medido'*state_ran_medido)
    state_ran_medido = State_sparse_complex(state_ran_medido, o)
    #e_medido = sort(eigvals(Matrix(rhoqsp(state_ran_medido))), rev=true)

    #y chequeamos mayorización
    global count = count + majorization_qsp(state_ran_medido, state_ran)
end

if count == repeticiones
    println("Measured state qsp majorizes unmeasured state")
elseif count == -repeticiones
    println("Unmeasured state qsp majorizes measured state")
else
    println("There is no qsp majorization")
end

=#

#This macro serves for printing as a string the variable's name
#=
macro Name(arg)
   string(arg)
end
=#
