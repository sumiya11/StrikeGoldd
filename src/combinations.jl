
"""
    lie_derivative(f, ode)

Input:
- `f` - rational function in states, parameters (not inputs) of `ode`
- `ode' - an ODE model

Output:
- Lie derivative of `f` with respect to `ode` 
"""
function lie_derivative(
    f,
    ode,
)
    res = zero(parent(f))
    for (x, eq) in ode.x_equations
        res += derivative(f, x) * eq
    end
    return res
end

function observability_matrix_symbolic(ode; with_respect_to=:parameters)
    ode_vars = if with_respect_to === :parameters
        ode.parameters
    else
        vcat(ode.x_vars, ode.parameters)
    end
    @info "Constructing the observability matrix with respect to $with_respect_to:\n$ode_vars"
    func = first(values(ode.y_equations))
    @debug "" func typeof(func)
    func_derivatives = [func]
    for i in 1:length(ode_vars)-1
        new_derivative = lie_derivative(last(func_derivatives), ode)
        push!(func_derivatives, new_derivative)
    end
    @debug "" ode_vars func_derivatives
    n = length(func_derivatives)
    R = parent(first(ode_vars))
    S = MatrixSpace(FractionField(R), n, n)
    J = zero(S)
    @debug "" n R
    for (i, f) in enumerate(func_derivatives)
        for (j, x) in enumerate(ode_vars)
            J[i, j] = derivative(f, x)
        end
    end
    @debug "" J
    return J, ode_vars
end

function matrix_kernel_symbolic(matrix)
    kerdim, ker = Nemo.kernel(matrix)
    solutions = [ker[:, i] for i in 1:kerdim]
    solutions = map(s -> s .* reduce(lcm, map(denominator, s)), solutions)
    @debug "" kerdim solutions
    @info "The observability matrix of size $(size(matrix)) has rank $(size(matrix, 1) - kerdim)"
    return kerdim, solutions
end

function generate_ansatz_polynomial(ode, ode_vars; degree=2)
    ring = ode.poly_ring
    all_variables = unique(vcat(ode.x_vars, ode_vars))
    n = length(all_variables)
    degree = 2
    new_vars = begin
        all_monoms = get_all_monoms_up_to_total_degree(ring, ode_vars, degree)
        ["r$j" for j in 1:length(all_monoms)]
    end

    @debug "" new_vars

    # ∂φ/∂a = r2 + 2 r4 a + r5 b
    # ∂φ/∂b = r3 + 2 r6 b + r5 a

    # -a (r2 + 2 r4 a + r5 b) + b (r3 + 2 r6 b + r5 a) = 0
    # a (-r2) + b (r3) + a^2 (-2 r4) + a b (-r5 + r5) + b^2 (2 r6) = 0
    # r2 = r3 = r4 = r6 = 0
    # r1 = Any
    # r2 = Any

    extended_ring, ext_vars = PolynomialRing(base_ring(ring), vcat(map(string, all_variables), new_vars))
    delta = length(all_variables) - length(ode_vars)
    ode_vars, unknown_coeffs = ext_vars[1+delta:n], ext_vars[n+1:end]
    all_monoms = get_all_monoms_up_to_total_degree(extended_ring, ode_vars, degree)
    @info "" ode_vars unknown_coeffs

    @debug "" all_monoms
    all_monoms = map(m -> parent_ring_change(m, extended_ring), all_monoms)
    derivatives = Dict()
    count = length(all_monoms)
    orig = zero(extended_ring)
    orig = sum([unknown_coeffs[j] * all_monoms[j] for j in 1:count])
    for (i, variable) in enumerate(ode_vars)
        derivatives[variable] = derivative(orig, variable)
    end
    @debug "" orig
    @info """
    Searching for identifiable functions of degree $degree of the form:
        $orig
    where $unknown_coeffs are the unknown coefficients.
    """

    return ode_vars, orig, derivatives, unknown_coeffs
end

function solve_linear_system(ode, ode_vars, derivatives, kernel_vectors, unknown_coeffs)
    #=
    ∂φ/∂x1 v1 + ... + ∂φ/∂xn vq = 0
    for v in kernel vectors

    v = (0, -1, 1)
    φ = c2 a + c3 b + c4 a*b

    0 - (c2 + c4*b) + (c3 + c4*a) = 0
    c3 - c2 = 0
    c4 = 0
    =#

    extended_ring = parent(first(ode_vars))

    @debug "" derivatives ode_vars unknown_coeffs

    system_eqs = []
    for v in kernel_vectors
        s = zero(extended_ring)
        for (i, variable) in enumerate(ode_vars)
            s += derivatives[variable] * parent_ring_change(v[i], extended_ring)
        end
        push!(system_eqs, s)
    end

    system_eqs = map(s -> (@assert isone(denominator(s)); numerator(s)), system_eqs)
    @debug "" system_eqs

    @debug """
    Constructed the following system in $unknown_coeffs:
    $(join(map(s -> "\t"*string(s) * " = 0", system_eqs), "\n"))"""

    system = []
    polyvars = collect(keys(derivatives))
    for eq in system_eqs
        s = collect(values(extract_coefficients(eq, polyvars)))
        s = map(k -> parent_ring_change(k, extended_ring), s)
        append!(system, s)
    end

    # TODO: A bug here !!
    @info """
    Transformed into a linear system:
    $(join(map(s -> "\t"*string(s) * " = 0", system), "\n"))"""
    @debug "" system

    kernel_dimension, solutions = solve_system(system, unknown_coeffs)
    @debug "" solutions

    @debug "There are $kernel_dimension independent solutions" solutions

    kernel_dimension, solutions
end

function construct_identifiable_functions(ode, ode_vars, orig, solutions, unknown_coeffs)
    @info "Constructing the identifiable functions, normalizing the coefficients to $(1)"
    extended_ring = parent(first(ode_vars))
    identifiable_functions = []
    ring = extended_ring
    polyvars = Nemo.gens(ring)
    n = length(ode.x_vars) + length(ode.parameters)
    for solution in solutions
        solution = Dict(parent_ring_change(k, ring) => v for (k, v) in solution)
        point = copy(polyvars)
        for var in polyvars[1:n]
            for (i, var__) in enumerate(polyvars)
                if var__ in unknown_coeffs
                    point[i] = var__
                end
                if haskey(solution, var__)
                    point[i] = ring(solution[var__])
                end
            end
        end
        func = evaluate(orig, point)
        if iszero(func)
            continue
        end
        push!(identifiable_functions, func)
    end
    identifiable_functions
end
