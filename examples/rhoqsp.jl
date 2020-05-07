using Fermionic
using SparseArrays
using LinearAlgebra


function kqsp(sta)
    precis = 15
    n = dim(ope(sta))
    k = spzeros(n,n)
    estate = st(sta)
    o = ope(sta)
    for i in 1:n
        for j in 1:n
            k[i,j] = round(estate'*cmcm(o, j, i)*estate, digits = precis)
        end
    end
    return k
end

function kqspc(sta)
    precis = 15
    n = dim(ope(sta))
    k = spzeros(Complex{Float64},n,n)
    estate = st(sta)
    o = ope(sta)
    for i in 1:n
        for j in 1:n
            k[i,j] = round(estate'*cmcm(o, j, i)*estate, digits = precis)
        end
    end
    return k
end

kcqsp2 = -conj(kqspc(sta))


#no es necesaria esta funcion
function kcqsp(sta)
    precis = 15
    n = dim(ope(sta))
    k = spzeros(n,n)
    estate = st(sta)
    o = ope(sta)
    for i in 1:n
        for j in 1:n
            k[i,j] = -round(estate'*cdcd(o, j, i)*estate, digits = precis)
        end
    end
    return k
end


#no es necesaria esta funcion
function kcqspc(sta)
    precis = 15
    n = dim(ope(sta))
    k = spzeros(Complex{Float64},n,n)
    estate = st(sta)
    o = ope(sta)
    for i in 1:n
        for j in 1:n
            k[i,j] = conj(round(estate'*cdcd(o, j, i)*estate, digits = precis))
        end
    end
    return k
end


function rhoqsp(sta)
    rho = [rhosp(sta) kqsp(sta); conj(kqsp(sta))' I-rhosp(sta)]
    return rho
end


op4 = Op(4)
cd1 = cdm(op4,1);
cd2 = cdm(op4,2);
cd3 = cdm(op4,3);
cd4 = cdm(op4,4);
vac = vacuum(op4)


#initzialize a slater determinant
state_ds0 = cd1*cd2*vac;
state_ds = State_sparse(state_ds0, op4) #or just State if you want to work with arrays

#initialize a state with no fixed fermionic number
state_qsp1 = (cd1*cd2+cd3*cd4+cd1*cd2*cd3*cd4)*vac + vac
state_qsp1 = state_qsp1/sqrt(state_qsp1'*state_qsp1)
state_qsp1 = State_sparse(state_qsp1, op4)

#initialize a state with no fixed fermionic number
state_qsp2 = (cd1*cd2+im*cd1*cd4+cd1*cd2*cd3*cd4)*vac
state_qsp2 = state_qsp2/sqrt(state_qsp2'*state_qsp2)
state_qsp2 = State_sparse_complex(state_qsp2, op4)

#initzialize a maximally entangled state
state_ent = (cd1*cd2 + cd3*cd4)*vac/sqrt(2)
state_ent = State_sparse(state_ent, op4); #or just State if you want to work with arrays


#initialization of a random state of 2 particles
#sparse real:
using SparseArrays

l= length(vac)
nume = 2 #number of particles
state_ran0 = spzeros(l)

for i in 1:l
    if sum(basis(op4)[i,:]) == nume
        state_ran0[i] = 2*rand(Float64,1)[1]-1
    end
end

state_ran = State_sparse(state_ran0, op4)

#random complex sparse
using SparseArrays

l= length(vac)
nume = 2 #number of particles
state_ran_com = spzeros(Complex{Float64},l)

for i in 1:l
    if sum(basis(op4)[i,:]) == nume
        state_ran_com[i] = 2*rand(Complex{Float64},1)[1]-1
    end
end
state_ran_com = state_ran_com/sqrt(state_ran_com'*state_ran_com)
state_ran_com = State_sparse_complex(state_ran_com, op4)



#random real array state
nume = 2
l = length(vac)

random_state = [0.0 for i in 1:l]

for i in 1:l
    if sum(basis(op4)[i,:]) == nume
        random_state[i] = 2*rand(Float64,1)[1]-1
    end
end
random_state = random_state/sqrt(random_state'*random_state)
random_state = State(random_state,op4)


#Random complex array state
nume = 2
l = length(vac)

random_state_com = zeros(Complex{Float64},l)

for i in 1:l
    if sum(basis(op4)[i,:]) == nume
        random_state_com[i] = 2*rand(Complex{Float64},1)[1]-1
    end
end
random_state_com = random_state_com/sqrt(random_state_com'*random_state_com)
random_state_com = State_complex(random_state_com,op4)
