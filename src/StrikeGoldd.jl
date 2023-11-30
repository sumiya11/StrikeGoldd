module StrikeGoldd

using AbstractAlgebra
using Nemo
using MacroTools
using Combinatorics

include("utils.jl")
include("ODE.jl")
include("stuff.jl")
include("combinations.jl")

#=
The algorithm:
    (1) Construct the observability matrix
    (2) Compute the kernel of the observability matrix
    (3) Generate an ansatz for an identifiable function
    (4) Solve the linear equations on the ansatz: the partial derivatives of the
        ansatz must be orthogonal to the kernel.

Loosely follows the algorithm from 
    An Algorithm for Generating Locally Identifiable Reparameterisations of
    Unidentifiable Systems
    https://doi.org/10.1016/S1474-6670(17)35488-5
=#
function find_some_identifiable_functions(ode)
    # NOTE: The efficiency of kernel finding can probably be improved by
    # applying evaluation and interpolation.
    # TODO: Doing some experiments would be nice.
    # TODO: Comparing this to `find_identifiable_functions` from
    # StructuralIdentifiability.jl would be nice.
    # NOTE: We can perhaps generalize the algorithm to handle the ODE states.
    # NOTE: Instead of searching for vectors orthogonal to the kernel of the
    # matrix, we can instead find a basis of its rowspace, which would give us
    # the partial derivatives of the identifiable functions.

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
    # - f1,...,fn is the vector field (given as polynomials in x).
    @info "Searching for some (at least locally) identifiable functions"

    # For now, assume that there is only one output
    @assert length(ode.y_equations) == 1

    ###
    # (1)
    # Construct the observability matrix of the form
    #   | ∂g/∂p1        ∂g/∂p2    ...      ∂g/∂pm  |
    #   | ∂Lg/∂p1      ∂Lg/∂p2    ...     ∂Lg/∂pm  |
    #                           ...
    #   | ∂L^r g/∂p1   ∂L^r g/∂p2 ...   ∂L^r g/∂pm |
    # where 
    #   - p1,...,pm are variables of interest. Usually, either parameters of the
    #     ODE, or parameters with the states adjoined.
    #   - L is the Lie derivative with respect to f1,...,fn
    #   - r is the maximal required order, r <= m
    obs_matrix, ode_vars = observability_matrix_symbolic(ode, with_respect_to=:parameters)

    ###
    # (2)
    # Find the kernel of the observability matrix over a fraction field that
    # includes parameters. 
    # The resulting kernel is given by its basis, the array of vectors
    #   v1,...,vq
    # such that for all i, j, we have
    #   vi * ∇ L^j g = 0
    kernel_dimension, kernel_vectors = matrix_kernel_symbolic(obs_matrix)

    if iszero(kernel_dimension)
        @info "The observability matrix is full rank."
        return nothing
    end

    ###
    # (3)
    # Generate an ansatz together with a system of equations on the unknown
    # coefficients of the ansatz. We shall search for identifiable functions as
    # the solutions to the equality
    #   vi * (∂A/∂p1 ∂A/∂p2 ... ∂A/∂pm) = 0,
    # where vi are the elements of the kernel, and A is the ansatz.
    ode_vars, ansatz_expr, derivatives, unknown_coeffs = generate_ansatz_polynomial(ode, ode_vars, degree=2)

    ###
    # (4)
    # Solve the linear system from (3) for the coefficients of the ansatz. Note
    # that the solution space can be infinite!
    solutions_dimension, solutions = solve_linear_system(ode, ode_vars, derivatives, kernel_vectors, unknown_coeffs)

    if iszero(solutions_dimension)
        @info "No identifiable functions found."
        return nothing
    end

    ###
    # (5)
    # Consrtuct the answer
    identifiable_functions = construct_identifiable_functions(ode, ode_vars, ansatz_expr, solutions, unknown_coeffs)

    return identifiable_functions
end

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
function find_some_symmetries(ode; type=:univariate, degree=2)
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
    infinitesimal_generator = create_infinitesimal_generator(ode, type=type, degree=degree)

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
    symmetries = integrate_infinitesimals(ode, infinitesimal_generator, solutions)

    @info "Found $(length(symmetries)) symmetries"
    @info "The search for symmetries concluded in $(round(time() - time_start, digits=3)) seconds"
    return symmetries
end

export @ODEmodel
export find_some_symmetries

end
