# Welcome to the Fermionic Package (Julia)

[![Build Status](https://travis-ci.com/Marco-Di-Tullio/Fermionic.jl.svg?branch=master)](https://travis-ci.com/Marco-Di-Tullio/Fermionic.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/Marco-Di-Tullio/Fermionic.jl?svg=true)](https://ci.appveyor.com/project/Marco-Di-Tullio/Fermionic-jl)
[![Codecov](https://codecov.io/gh/Marco-Di-Tullio/Fermionic.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/Marco-Di-Tullio/Fermionic.jl)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/Marco-Di-Tullio/Fermionic.jl/master)

Everything relating to fermions can be expressed in terms of annilhation and creation operators. This package numerically defines fermionic operators and can therefore be used for building any fermionic operator and state. It also defines the main Quantum Gates.

The only input the program needs is the dimension, and the corresponding one body fermionic operators will be automatically defined as sparse matrices. Then you can define states in the base of possible states of the specified dimension, and calculate some of its properties.

You can perform any unitary operation with this package. You can also use some of the most commons logical gates (Pauli matrices, phase shift, Hadamard, Ucnot, SWAP).

For installing this package, you must first acces the pkg REPL (by typing ']' in your command line) and then execute

```add Fermionic```

The pkg manager will automatically download the package. Then you can initialize it by typing

```using Fermionic```

Alternatively, you can install the package from an editor/Jupyter notebook by typing

```import Pkg```
```Pkg.add("Fermionic")```

For instructions on how to use this package, you can read the tutorials located in the folder 'examples\'

![](/images/quantuminfo.png)
