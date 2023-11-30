# StrikeGoldd

**StrikeGoldd** is a sandbox Julia package for analyzing structural identifiability of dynamical models described by ordinary differential equations (ODEs). 

**StrikeGoldd** is (to an extent) a direct port of the MATBLAB package [STRIKE-GOLDD](https://github.com/afvillaverde/strike-goldd). The package [STRIKE-GOLDD](https://github.com/afvillaverde/strike-goldd) is licensed under the GNU GPL Version 3. 

## Installation

To install **StrikeGoldd**, execute this command in the Julia REPL

```julia
import Pkg; Pkg.add(url="https://github.com/sumiya11/StrikeGoldd")
```

## Usage Example

**StrikeGoldd** provides two functions: `find_some_identifiable_functions` and `find_some_symmetries`.

#### `find_some_identifiable_functions`

**StrikeGoldd** provides the function `find_some_identifiable_functions`, which can be applied to polynomial ODE models to find some at least locally identifiable functions of parameters. For example, you can run in Julia:

```julia
using StrikeGoldd     # load the package

ode = @ODEmodel(      # create the ODE
    x1'(t) = x1*x2 + a*b,
    x2'(t) = x2 + a*c,
    y(t) = x1
)

funcs = find_some_identifiable_functions(ode)
```

Running this code will print:

```julia
[ Info: Searching for some (at least locally) identifiable functions
┌ Info: Constructing the observability matrix with respect to parameters:
└ Nemo.QQMPolyRingElem[a, b, c]
[ Info: The observability matrix of size (3, 3) has rank 2
┌ Info: Searching for identifiable functions of degree 2 of the form:
│     a^2*r5 + a*b*r6 + a*c*r7 + a*r2 + b^2*r8 + b*c*r9 + b*r3 + c^2*r10 + c*r4 + r1        
└ where Nemo.QQMPolyRingElem[r1, r2, r3, r4, r5, r6, r7, r8, r9, r10] are the unknown coefficients.
┌ Info: Constructed the following system in Nemo.QQMPolyRingElem[r1, r2, r3, r4, r5, r6, r7, r8, r9, r10]:
└       -2*a^2*r5 - a*r2 + 2*b^2*r8 + 2*b*c*r9 + b*r3 + 2*c^2*r10 + c*r4 = 0
┌ Info: Transformed into a linear system:
│       r4 = 0
│       -2*r5 = 0
│       2*r8 = 0
│       -r2 = 0
│       r3 = 0
│       2*r9 = 0
└       2*r10 = 0
┌ Info: There are 3 independent solutions
│   solutions =
│    3-element Vector{Dict{Nemo.QQMPolyRingElem, Nemo.QQFieldElem}}:
│     Dict(r9 => 0, r10 => 0, r1 => 1, r4 => 0, r3 => 0, r2 => 0, r5 => 0, r6 => 0, r7 => 0, r8 => 0…)
│     Dict(r9 => 0, r10 => 0, r1 => 0, r4 => 0, r3 => 0, r2 => 0, r5 => 0, r6 => 1, r7 => 0, r8 => 0…)
└     Dict(r9 => 0, r10 => 0, r1 => 0, r4 => 0, r3 => 0, r2 => 0, r5 => 0, r6 => 0, r7 => 1, r8 => 0…)
[ Info: Constructing the identifiable functions, normalizing the coefficients to 1
3-element Vector{Any}:
 1
 a*b
 a*c
```

#### `find_some_symmetries`:

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
