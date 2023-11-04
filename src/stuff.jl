
# Returns a map from variable xi to an infinitesimal ηi 
function create_infinitesimal_generator(ode; degree=2)
    ring = ode.poly_ring
    x_vars, parameters = ode.x_vars, ode.parameters
    all_variables = vcat(x_vars, parameters)
    n = length(all_variables)
    new_vars = ["r$i$j" for i in 1:n for j in 0:degree]
    extended_ring, ext_vars = PolynomialRing(base_ring(ring), vcat(map(string, all_variables), new_vars, map(string, ode.u_vars)))
    ode_vars, unknown_coeffs = ext_vars[1:n], ext_vars[n+1:end-length(ode.u_vars)]
    @info """
    ODE has $(length(x_vars)) states and $(length(parameters)) parameters.
    Creating univariate infinitesimals of degree $degree.
    """
    @info """
    Introducing $(length(unknown_coeffs)) unknown coefficients:
    [$(join(map(string, unknown_coeffs), ','))]."""
    # Maps xi to ηi
    infinitesimals = Dict()
    for (i, variable) in enumerate(ode_vars)
        infinitesimals[variable] = sum([unknown_coeffs[(i-1)*(degree+1)+j+1] * variable^j for j in 0:1:degree])
    end
    @debug """
    Inifinitesimals:
    $infinitesimals"""
    return infinitesimals
end

# Returns a map from variable xi to ∂ηi/∂xi
function extend_infinitesimal_generator(infinitesimals)
    extension = empty(infinitesimals)
    for (variable, infinitesimal) in infinitesimals
        extension[variable] = derivative(infinitesimal, variable)
    end
    @debug """
    Extension:
    $extension"""
    return extension
end

# Returns a system to be solved
function invariance_criterion_system(ode, infinitesimals, extension)
    x_equations, y_equations = ode.x_equations, ode.y_equations
    ext_ring = parent(first(keys(infinitesimals)))
    x_equations = Dict(parent_ring_change(k, ext_ring) => parent_ring_change(v, ext_ring) for (k, v) in x_equations)
    y_equations = Dict(k => parent_ring_change(v, ext_ring) for (k, v) in y_equations)
    X_y = []
    for (_, eq) in y_equations
        poly = sum([infinitesimals[var] * derivative(eq, var) for var in keys(infinitesimals)])
        push!(X_y, poly)
    end
    dX_x = []
    for (_var, eq) in x_equations
        poly = sum([infinitesimals[var] * derivative(eq, var) for var in keys(infinitesimals)])
        ext_poly = sum([x_equations[_var] * derivative(infinitesimals[_var], var) for var in keys(infinitesimals)])
        push!(dX_x, ext_poly - poly)
    end
    @debug """
    X (y - g(x)) = 0 is equivalent to
    $X_y = 0

    X' (x' - f(x)) = 0 is equivalent to
    $dX_x = 0"""
    system = []
    polyvars = collect(keys(infinitesimals))
    for eq in X_y
        append!(system, collect(values(extract_coefficients(eq, polyvars))))
    end
    for eq in dX_x
        append!(system, collect(values(extract_coefficients(eq, polyvars))))
    end
    @debug """
    Constructed the system
    $(join(map(x -> string(x) * " = 0", system), "\n"))"""
    @info """
    Constructed the linear system with $(length(system)) equations"""
    return system
end

function solve_system(system)
    monoms = reduce(union, map(f -> collect(monomials(f)), system))
    sort!(monoms, rev=true)
    n = length(monoms)
    m = length(system)
    @debug "Existing monomials in unknowns are:\n$monoms"
    system_numeric = map(f -> map(monom -> coeff(f, monom), monoms), system)
    system_numeric = reduce(hcat, system_numeric)
    S = Nemo.MatrixSpace(base_ring(parent(monoms[1])), m, n)
    system_nemo = S(Matrix(transpose(system_numeric)))
    @debug """
    Finding the kernel of
    $(join(map(x -> string(x), [system_nemo[i, :] for i in 1:m]), "\n"))"""
    kerdim, ker = kernel(system_nemo)
    @info "The dimension of the right-kernel is $kerdim."
    @debug "The kernel is\n$ker"
    solutions = [Dict(monoms .=> ker[:, i]) for i in 1:kerdim]
    @debug "System solutions are\n$solutions"
    kerdim, solutions
end

function integrate_infinitesimals(ode, infinitesimals, solutions)
    infinitesimal_generators = []
    ring = parent(first(keys(infinitesimals)))
    polyvars = Nemo.gens(ring)
    n = length(ode.x_vars) + length(ode.parameters)
    for solution in solutions
        solution = Dict(parent_ring_change(k, ring) => v for (k, v) in solution)
        inft = []
        for var in polyvars[1:n]
            inft_ = infinitesimals[var]
            point = copy(polyvars)
            for (i, var__) in enumerate(polyvars)
                if haskey(solution, var__)
                    point[i] = ring(solution[var__])
                end
            end
            inft_eval = evaluate(inft_, point)
            push!(inft, inft_eval)
        end
        push!(infinitesimal_generators, inft)
    end
    symbs = map(i -> "∂/∂$(polyvars[i])", 1:n)
    s = map(gen -> join(map(x -> "(" * string(x) * ")", gen) .* symbs, " + "), infinitesimal_generators)
    for i in 1:length(s)
        s[i] = "($i)\t" * s[i]
    end
    @info """
    Found $(length(infinitesimal_generators)) infinitesimal generators:
    $(join(s, "\n"))"""
    newring, newvars = PolynomialRing(base_ring(ring), vcat(map(string, polyvars), "ε"))
    newvars, epsilon = newvars[1:end-1], newvars[end]
    symmetries = []
    for (idx, gen) in enumerate(infinitesimal_generators)
        @info "Integrating infinitesimal generator ($(idx))"
        transform = Dict()
        skip = false
        for i in 1:n
            # If the infinitesimal is zero, then there is no change
            if iszero(gen[i])
                int = parent_ring_change(polyvars[i], newring)
                transform[polyvars[i]] = int
                continue
            end
            cfs = extract_coefficients(gen[i], [polyvars[i]])
            # If more than 1 coefficient
            if length(cfs) > 1
                skip = true
                @info "Cannot integrate generator ($idx): too large"
                break
            end
            cf = first(values(cfs))
            if degree(gen[i], polyvars[i]) == 0
                @debug "Integrated ($idx), infinitesimal $i: translation"
                int = parent_ring_change(polyvars[i], newring)
                int = int + epsilon * parent_ring_change(cf, newring)
                transform[polyvars[i]] = int
            elseif degree(gen[i], polyvars[i]) == 1
                # @debug "Integrated ($idx), infinitesimal $i: scaling"
                # int = parent_ring_change(polyvars[i], newring)
                # int = int * epsilon * parent_ring_change(cf, newring)
                # transform[polyvars[i]] = int
                @info "Cannot integrate ($idx): exponentiation is not supported :("
                skip = true
                break
            else
                skip = true
                @info "Cannot integrate ($idx): too large"
                break
            end
        end
        if skip
            continue
        end
        push!(symmetries, transform)
    end
    return symmetries
end
