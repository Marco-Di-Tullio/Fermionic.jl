# Welcome to Fermionic.jl

[![Build Status](https://travis-ci.com/Marco-Di-Tullio/Fermionic.jl.svg?branch=master)](https://travis-ci.com/Marco-Di-Tullio/Fermionic.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/Marco-Di-Tullio/Fermionic.jl?svg=true)](https://ci.appveyor.com/project/Marco-Di-Tullio/Fermionic-jl)
[![Codecov](https://codecov.io/gh/Marco-Di-Tullio/Fermionic.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/Marco-Di-Tullio/Fermionic.jl)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/Marco-Di-Tullio/Fermionic.jl/master)


Fermionic is a Julia toolkit for implementing fermionic simulations and exploring its quantum information properties.

Everything relating to fermions can be expressed in terms of annilhation and creation operators. This package numerically constructs fermionic operators and can therefore be used for building any fermionic operator and state. It also defines the main Quantum Gates.

The corresponding fermionic operators can be constructed both in the full Fock space or in the fixed particle number subspace. Then you can define states in the corresponding base and calculate several properties. 

Many interesting quantities can be obtained from states in Fermionic, such as one body matrices entropy, partially traced systems, m-bodies density matrices, one body entropies, majorization relations, average particle number and more.

You can also perform any unitary operation with this package and use some of the most commons logical gates (Pauli matrices, phase shift, Hadamard, Ucnot, SWAP) to perform fermionic quantum computation.

## Installation

For installing this package, you must first acces the pkg REPL (by typing ']' in your command line) and then execute

```add Fermionic```

The pkg manager will automatically download the package. Then you can initialize it by typing

```using Fermionic```

Alternatively, you can install the package from an editor/Jupyter notebook by typing

```import Pkg```

```Pkg.add("Fermionic")```

```using Fermionic```

For instructions on how to use this package, you can read the tutorials located in the folder 'examples\'

![](/images/quantuminfo.png)
