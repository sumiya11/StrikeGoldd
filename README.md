# StrikeGoldd

**StrikeGoldd** is sandbox Julia package for analyzing structural identifiability of dynamical models described by ordinary differential equations (ODEs). 

**StrikeGoldd** is (to an extent) a direct port of the MATBLAB package [STRIKE-GOLDD](https://github.com/afvillaverde/strike-goldd). The package [STRIKE-GOLDD](https://github.com/afvillaverde/strike-goldd) is licensed under the GNU GPL Version 3. 

## Installation

In the Julia REPL, execute

```julia
import Pkg; Pkg.add(url="")
```

## Usage Example

This package provides the function `find_some_symmetries`. It can be used like this:

```julia
ode = @ODEmodel(
    ...
)

symmetries = find_some_symmetries(ode)
```

This will print the following:

```julia

```
