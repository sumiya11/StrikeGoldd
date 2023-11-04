# StrikeGoldd

**StrikeGoldd** is sandbox Julia package for analyzing structural identifiability of dynamical models described by ordinary differential equations (ODEs). 

**StrikeGoldd** is (to an extent) a direct port of the MATBLAB package [STRIKE-GOLDD](https://github.com/afvillaverde/strike-goldd). The package [STRIKE-GOLDD](https://github.com/afvillaverde/strike-goldd) is licensed under the GNU GPL Version 3. 

## Installation

To install **StrikeGoldd**, execute this command in the Julia REPL

```julia
import Pkg; Pkg.add(url="https://github.com/sumiya11/StrikeGoldd")
```

## Usage Example

**StrikeGoldd** provides the function `find_some_symmetries`, which can be applied to polynomial ODE models. For example:

```julia
using StrikeGoldd   # load the package

ode = @ODEmodel(    # create the ODE
    x1'(t) = x1(t) + p1 + p2,
    y(t) = x1(t)
)

symmetries = find_some_symmetries(ode)
```

Running the above will print, among other things, the following:

```julia

```
