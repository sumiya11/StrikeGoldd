# StrikeGoldd

**StrikeGoldd** is sandbox Julia package for analyzing structural identifiability of dynamical models described by ordinary differential equations (ODEs). 

**StrikeGoldd** is (to an extent) a direct port of the MATBLAB package [STRIKE-GOLDD](https://github.com/afvillaverde/strike-goldd). The package [STRIKE-GOLDD](https://github.com/afvillaverde/strike-goldd) is licensed under the GNU GPL Version 3. 

## Installation

To install **StrikeGoldd**, execute this command in the Julia REPL

```julia
import Pkg; Pkg.add(url="https://github.com/sumiya11/StrikeGoldd")
```

## Usage Example

**StrikeGoldd** provides the function `find_some_symmetries`, which can be applied to polynomial ODE models. For example, you can run in Julia:

```julia
using StrikeGoldd   # load the package

ode = @ODEmodel(    # create the ODE
    x1'(t) = x1(t) + p1 + p2,
    y(t) = x1(t)
)

symmetries = find_some_symmetries(ode)
```

Running this code will print, among other things, the following:

```julia
[ Info: Start searching for symmetries
[ Info: Creating the infinitesimal generator
┌ Info: ODE has 1 states and 2 parameters.
└ Creating univariate infinitesimals of degree 2.
┌ Info: Introducing 9 unknown coefficients:
└ [r10,r11,r12,r20,r21,r22,r30,r31,r32].
[ Info: Extending the infinitesimal generator
[ Info: Constructing the linear system to solve
[ Info: Constructed the linear system with 11 equations
[ Info: Solving the system
[ Info: The dimension of the right-kernel is 1.
┌ Info: Found 1 infinitesimal generators:
└ (1)   (0)∂/∂x1 + (-1)∂/∂p2 + (1)∂/∂p1
[ Info: Integrating infinitesimal generator (1)
[ Info: Found 1 symmetries
[ Info: The search for symmetries concluded in 5.571 seconds
1-element Vector{Any}:
 Dict{Any, Any}(p1 => p1 + ε, p2 => p2 - ε, x1 => x1)
```
