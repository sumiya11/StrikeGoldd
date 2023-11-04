module StrikeGoldd

using AbstractAlgebra
using Nemo
using MacroTools

include("utils.jl")
include("ODE.jl")
include("stuff.jl")

"""
    find_some_symmetries(ode)

Tries to find some Lie symmetries of `ode`.

## Input

- `ode`: an object created by the `@ODEmodel` macro, the input ODE
    system

## Output

Returns a (possibly empty) array of variable transformations, each variable
transformation encodes a symmetry.

## Example

Running the following code

```jldoctest
using StrikeGoldd

ode = @ODEmodel(
    x1'(t) = x1(t) + p1 + p2,
    y(t) = x1(t)
)

find_some_symmetries(ode)
```

will print on the console

```jldoctest
1-element Vector{Any}:
 Dict{Any, Any}(p1 => p1 + ε, p2 => p2 - ε, x1 => x1)
```
"""
function find_some_symmetries(ode)
    ###
    # (0)
    # The input system of differential equations
    #   x1' = f1(x)
    #   ...
    #   xn' = fn(x)
    #   y  = g(x)
    # where
    # - x is a vector of state variables and time-independent parameters
    # - g are the output functions (given as polynomials in x)
    # - f1,...,fn is the vector field (given as polynomials in x)
    @info "Start searching for symmetries"
    time_start = time()

    ###
    # (1)
    # Create infinitesimal generator of form 
    #   X = η1 ∂/∂x1 + ... + ηn ∂/∂xn
    # We take infinitesimals η1, ..., ηn from a specific family of (polynomial)
    # functions parametrized by unknown parameters. The goal is to find such
    # values for parameters that a certain invariance condition is satisfied, or
    # prove that there are none
    @info "Creating the infinitesimal generator"
    infinitesimal_generator = create_infinitesimal_generator(ode)

    ###
    # (2)
    # Produce the extended infinitesimal generator X' of form
    #   X' = X + x1' ∂η1/∂x1 + ... + xn' ∂ηn/∂xn
    @info "Extending the infinitesimal generator"
    infinitesimal_extension = extend_infinitesimal_generator(infinitesimal_generator)

    ###
    # (3)
    # Find solutions of the system induced by the invariance criterion
    #   X' (x1' - f1(x)) = 0
    #   ...
    #   X' (xn' - fn(x)) = 0
    #   X  (y - g(x))    = 0
    @info "Constructing the linear system to solve"
    system = invariance_criterion_system(ode, infinitesimal_generator, infinitesimal_extension)

    ###
    # (4)
    # Solve the system from (3).
    # Each solution of the system gives rise to one infinitesimal generator X,
    # or, equivalently, to one Lie symmetry of the input ode
    @info "Solving the system"
    kernel_dimension, solutions = solve_system(system)

    if iszero(kernel_dimension)
        @info "No symmetries found."
        return nothing
    end

    ###
    # (5)
    # Integrate the infinitesimal generators by solving the IVP for x*:
    #   dx1*/dε = η1
    #   ...
    #   dxn*/dε = ηn
    #   xi* = xi when ε = 0
    # We handle only some simple types of integrands η, and abort the
    # integration otherwise. As a solution, we obtain the sought group of
    # transformations x1*,...,xn*
    @info "Integrating infinitesimals"
    symmetries = integrate_infinitesimals(ode, infinitesimal_generator, solutions)

    @info "Found $(length(symmetries)) symmetries"
    @info "The search for symmetries concluded in $(round(time() - time_start, digits=3)) seconds"
    return symmetries
end

export @ODEmodel
export find_some_symmetries

end
