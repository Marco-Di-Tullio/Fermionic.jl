using Fermionic

o = Op(6)
vac = vacuum(o)

cd1 = cdm(o,1)
cd2 = cdm(o,2)
cd3 = cdm(o,3)
cd4 = cdm(o,4)
cd5 = cdm(o,5)
cd6 = cdm(o,6)

p1 = sqrt(0.9)
p2 = sqrt(1-p1^2)
a = 1/sqrt(2)
b = sqrt(1-a^2)

state0 = (a*cd5 + b*cd6)*(p1*cd1*cd3+p2*cd2*cd4)*vac

state1 = ucnot(o, 5, 1)*ucnot(o, 5, 2)*state0

state2 = hadamard(o, 5, 6)*state1

state3 = ucnot(o, 5, 3)*ucnot(o, 5, 4)*state2
