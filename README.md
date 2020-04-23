# Welcome to the Fermionic Package

[![Build Status](https://travis-ci.com/Marco-Di-Tullio/Fermionic.jl.svg?branch=master)](https://travis-ci.com/Marco-Di-Tullio/Fermionic.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/Marco-Di-Tullio/Fermionic.jl?svg=true)](https://ci.appveyor.com/project/Marco-Di-Tullio/Fermionic-jl)
[![Codecov](https://codecov.io/gh/Marco-Di-Tullio/Fermionic.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/Marco-Di-Tullio/Fermionic.jl)


This Julia package numerically defines fermionic operators through a sparse matrix representation. The only input the program needs is the dimension, and the corresponding one body fermionic operators will be automatically defined. Then you can define states in the base of possible states of the specified dimension, and calculate some of its properties. 

For installing this package, you must first acces the pkg REPL (by typing ']' in your command line) and then execute

```add Fermionic```

The pkg manager will automatically download the package. Then you can initialize it by typing 

```using Fermionic```

For instructions on how to use this package, you can read the 'Fermionic tutorial' located in the folder 'examples\'

![](/images/quantuminfo.png)
