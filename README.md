# Welcome to Fermionic.jl

[![Build Status](https://travis-ci.com/Marco-Di-Tullio/Fermionic.jl.svg?branch=master)](https://travis-ci.com/Marco-Di-Tullio/Fermionic.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/Marco-Di-Tullio/Fermionic.jl?svg=true)](https://ci.appveyor.com/project/Marco-Di-Tullio/Fermionic-jl)
[![Codecov](https://codecov.io/gh/Marco-Di-Tullio/Fermionic.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/Marco-Di-Tullio/Fermionic.jl)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/Marco-Di-Tullio/Fermionic.jl/master)


Fermionic is a Julia toolkit for implementing fermionic simulations and exploring its quantum information properties.

Everything relating to fermions can be expressed in terms of annihilation and creation operators. This package numerically constructs fermionic operators and can therefore be used for building any fermionic operator and state. It also defines the main Quantum Gates.

The corresponding fermionic operators can be constructed both in the full Fock space or in the fixed particle number subspace, which is far more efficient and _the secret weapon_ of this library. Then you can define states in the corresponding base and calculate several properties. 

Many interesting quantities can be obtained from states in Fermionic, such as one body matrices entropy, partially traced systems, m-bodies density matrices, one body entropies, majorization relations, average particle number and more.

You can also apply any unitary operation with this package and use some of the most common logical gates (Pauli matrices, phase shift, Hadamard, Ucnot, SWAP) to perform fermionic quantum computation.


## Installation

For installing this package, you must first access the pkg REPL (by typing ']' in your command line) and then execute

```add Fermionic```

The pkg manager will automatically download the package. Then you can initialize it by typing

```using Fermionic```

Alternatively, you can install the package from an editor/Jupyter notebook by typing

```import Pkg```

```Pkg.add("Fermionic")```

```using Fermionic```

## Getting started

For instructions on how to use this package, you can read the tutorials located in the folder 'examples\'.  For a quick preview, check out the following snippet for a simple example of Fermionic in action:


```julia
using Fermionic

# We initialize the fermionic operators in dimension d
d = 4
o = Op(d)

# We build the state by applying creation operators on the vacuum
vac = vacuum(o)
st = ad(o,1)*ad(o,2)*ad(o,3)*vac
```
Output
```julia
16-element SparseArrays.SparseVector{Float64,Int64} with 1 stored entry:
  [15]  =  1.0
```
as we are working on the Fock Space, were

```julia
Matrix(basis(o))[15,:]
```

is equal to

```julia
4-element Array{Float64,1}:
 1.0
 1.0
 1.0
 0.0
 ```
 
A more interesting exaple: Simulating a superconductor in the fixed particle subspace.

```julia
using Fermionic

#=
We will construct the superconducting Hamiltonian, obtain its fundamental state and some property
=#

# This system will have 4 levels with a double degeneracy
d = 8
nume = Int(d/2)

o = Op_fixed(d,nume)

# We select the energy of each level and the coupling between them
e0 = 1.0
g = 5.0

epsilon = [e0*(i-d/4-1/2) for i in 1:d/2]
epsilon = sort([epsilon; epsilon])

# The non interacting Hamiltonian
h0 = sum([epsilon[i]*(ada(o,i,i) + ada(o,i+1,i+1)) for i in 1:2:(Int(d)-1)]) 

# The interacting Hamiltonian
hi = sum([sum([if i==j spzeros(binomial(d,nume), binomial(d,nume)) else -(ada(o,j,i+1)*ada(o,j+1,i)) end
                    for i in 1:2:(Int(d)-1)]) for j in 1:2:(Int(d)-1)]) 

# The full Hamiltonian
h = h0 - g*hi

# We compute the fundamental state of the full Hamiltonian
fundamental = eigvecs(Matrix(h))[:,1]

# We normalize the state
fundamental = fundamental/sqrt(fundamental'*fundamental)

# We initialize a state with fundamental
fund = State_fixed(fundamental,o)

# We can now access some properties, for instance the eigenvalues of the one body matrix
eigensp(fund)
```

Output

```julia
8-element Array{Float64,1}:
 0.598475566126938
 0.598475566126938
 0.534545046798091
 0.534545046798091
 0.465454953201909
 0.465454953201909
 0.401524433873062
 0.401524433873062
```

This can also be done for varying couplings or even symbolically (non fixed g and e), as you can see in 'examples\6. Solving a superconducting system (Example)'. 

![](/images/quantuminfo.png)
