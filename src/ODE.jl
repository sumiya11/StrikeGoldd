# This file is adapted from StructuralIdentifiability.jl.
# The license is MIT.
# Copyright (c) 2020, R. Dong, C. Goodbrake, H. Harrington, G. Pogudin

# P is the type of polynomials in the rhs of the ODE system
"""
The main structure that represents input ODE system.

Stores information about states (`x_vars`), outputs (`y_vars`), inputs (`u_vars`), parameters (`parameters`) and the equations.

This structure is constructed via `@ODEmodel` macro.
"""
struct ODE{P}
    poly_ring::MPolyRing
    x_vars::Array{P,1}
    y_vars::Array{P,1}
    u_vars::Array{P,1}
    parameters::Array{P,1}
    x_equations::Dict{P,<:Union{P,Generic.Frac{P}}}
    y_equations::Dict{P,<:Union{P,Generic.Frac{P}}}

    function ODE{P}(
        x_vars::Array{P,1},
        y_vars::Array{P,1},
        x_eqs::Dict{P,<:Union{P,Generic.Frac{P}}},
        y_eqs::Dict{P,<:Union{P,Generic.Frac{P}}},
        inputs::Array{P,1},
    ) where {P<:MPolyElem{<:FieldElem}}
        # Initialize ODE
        # x_eqs is a dictionary x_i => f_i(x, u, params)
        # y_eqs is a dictionary y_i => g_i(x, u, params)
        if isempty(y_eqs)
            @info "Could not find output variables in the model."
        end
        poly_ring = parent(first(vcat(y_vars, x_vars)))
        u_vars = inputs
        parameters = filter(
            v -> (!(v in x_vars) && !(v in u_vars) && !(v in y_vars)),
            gens(poly_ring),
        )
        new{P}(poly_ring, x_vars, y_vars, u_vars, parameters, x_eqs, y_eqs)
    end

    function ODE{P}(
        x_eqs::Dict{P,<:Union{P,Generic.Frac{P}}},
        y_eqs::Dict{P,<:Union{P,Generic.Frac{P}}},
        inputs::Array{P,1},
    ) where {P<:MPolyElem{<:FieldElem}}
        x_vars = collect(keys(x_eqs))
        y_vars = collect(keys(y_eqs))
        return ODE{P}(x_vars, y_vars, x_eqs, y_eqs, inputs)
    end
end

#------------------------------------------------------------------------------

function _extract_aux!(funcs, all_symb, eq, ders_ok=false)
    aux_symb = Set([:(+), :(-), :(=), :(*), :(^), :t, :(/), :(//)])
    MacroTools.postwalk(
        x -> begin
            if @capture(x, f_'(t))
                if !ders_ok
                    throw(
                        Base.ArgumentError(
                            "Derivative are not allowed in the right-hand side",
                        ),
                    )
                end
                push!(all_symb, f)
            elseif @capture(x, f_(t))
                push!(funcs, f)
            elseif (x isa Symbol) && !(x in aux_symb)
                push!(all_symb, x)
            end
            return x
        end,
        eq,
    )
end

"""
  For an expression of the form f'(t) or f(t) returns (f, true) and (f, false), resp
"""
function _get_var(expr)
    if @capture(expr, f_'(t))
        return (f, true)
    end
    if @capture(expr, f_(t))
        return (f, false)
    end
    error("cannot extract the single function name from $expr")
end

function macrohelper_extract_vars(equations::Array{Expr,1})
    funcs, all_symb = Set(), Set()
    x_vars, y_vars = Vector(), Vector()
    aux_symb = Set([:(+), :(-), :(=), :(*), :(^), :t, :(/), :(//)])
    for eq in equations
        if eq.head != :(=)
            _extract_aux!(funcs, all_symb, eq)
        else
            lhs, rhs = eq.args[1:2]
            _extract_aux!(funcs, all_symb, lhs, true)
            _extract_aux!(funcs, all_symb, rhs)
            (v, is_state) = _get_var(lhs)
            if is_state
                push!(x_vars, v)
            else
                push!(y_vars, v)
            end
        end
    end
    u_vars = setdiff(funcs, vcat(x_vars, y_vars))
    all_symb = collect(all_symb)
    return x_vars, y_vars, collect(u_vars), collect(all_symb)
end

function macrohelper_extract_vars(equations::Array{Symbol,1})
    return macrohelper_extract_vars(map(Expr, equations))
end

#------------------------------------------------------------------------------

function macrohelper_clean(ex::Expr)
    ex = MacroTools.postwalk(x -> @capture(x, f_'(t)) ? f : x, ex)
    ex = MacroTools.postwalk(x -> @capture(x, f_(t)) ? f : x, ex)
    ex = MacroTools.postwalk(x -> x == :(/) ? :(//) : x, ex)
    ex = MacroTools.postwalk(x -> x isa Float64 ? rationalize(x) : x, ex)
    return ex
end

#------------------------------------------------------------------------------

"""
    macro ODEmodel

Macro for creating an ODE from a list of equations.
It also injects all variables into the global scope.

## Example

Creating a simple `ODE`:

```jldoctest
using StrikeGoldd

ode = @ODEmodel(
    x1'(t) = a * x1(t) + u(t),
    x2'(t) = b * x2(t) + c*x1(t)*x2(t),
    y(t) = x1(t)
)
```

Here,
- `x1`, `x2` are state variables
- `y` is an output variable
- `u` is an input variable
- `a`, `b`, `c` are time-indepdendent parameters

"""
macro ODEmodel(ex::Expr...)
    equations = [ex...]
    x_vars, y_vars, u_vars, all_symb = macrohelper_extract_vars(equations)

    # creating the polynomial ring
    vars_list = :([$(all_symb...)])
    R = gensym()
    vars_aux = gensym()
    exp_ring = :(
        ($R, $vars_aux) = StrikeGoldd.Nemo.PolynomialRing(
            StrikeGoldd.Nemo.QQ,
            map(string, $all_symb),
        )
    )
    assignments = [:($(all_symb[i]) = $vars_aux[$i]) for i in 1:length(all_symb)]

    # setting x_vars and y_vars in the right order
    vx = gensym()
    vy = gensym()
    x_var_expr = :($vx = Vector{StrikeGoldd.Nemo.fmpq_mpoly}([$(x_vars...)]))
    y_var_expr = :($vy = Vector{StrikeGoldd.Nemo.fmpq_mpoly}([$(y_vars...)]))

    # preparing equations
    equations = map(macrohelper_clean, equations)
    x_dict = gensym()
    y_dict = gensym()
    x_dict_create_expr = :(
        $x_dict = Dict{
            StrikeGoldd.Nemo.fmpq_mpoly,
            Union{
                StrikeGoldd.Nemo.fmpq_mpoly,
                StrikeGoldd.AbstractAlgebra.Generic.Frac{
                    StrikeGoldd.Nemo.fmpq_mpoly,
                },
            },
        }()
    )
    y_dict_create_expr = :(
        $y_dict = Dict{
            StrikeGoldd.Nemo.fmpq_mpoly,
            Union{
                StrikeGoldd.Nemo.fmpq_mpoly,
                StrikeGoldd.AbstractAlgebra.Generic.Frac{
                    StrikeGoldd.Nemo.fmpq_mpoly,
                },
            },
        }()
    )
    eqs_expr = []
    for eq in equations
        if eq.head != :(=)
            throw("Problem with parsing at $eq")
        end
        lhs, rhs = eq.args[1:2]
        loc_all_symb = macrohelper_extract_vars([rhs])[4]
        to_insert = undef
        if lhs in x_vars
            to_insert = x_dict
        elseif lhs in y_vars
            to_insert = y_dict
        else
            throw("Unknown left-hand side $lhs")
        end

        uniqueness_check_expr = quote
            if haskey($to_insert, $lhs)
                throw(
                    DomainError(
                        $lhs,
                        "The variable occurs twice in the left-hand-side of the ODE system",
                    ),
                )
            end
        end
        push!(eqs_expr, uniqueness_check_expr)
        if isempty(loc_all_symb)
            push!(eqs_expr, :($to_insert[$lhs] = $R($rhs)))
        else
            push!(eqs_expr, :($to_insert[$lhs] = ($rhs)))
        end
    end

    params = setdiff(all_symb, union(x_vars, y_vars, u_vars))
    allnames = map(
        string,
        vcat(collect(x_vars), collect(params), collect(u_vars), collect(y_vars)),
    )
    for n in allnames
        if !Base.isidentifier(n)
            throw(
                ArgumentError(
                    "The names of the variables will be injected into the global scope, so their name must be allowed Julia names, $n is not",
                ),
            )
        end
    end
    logging_exprs = [
        :(@info "Summary of the model:"),
        :(@info "State variables: " * $(join(map(string, collect(x_vars)), ", "))),
        :(@info "Parameters: " * $(join(map(string, collect(params)), ", "))),
        :(@info "Inputs: " * $(join(map(string, collect(u_vars)), ", "))),
        :(@info "Outputs: " * $(join(map(string, collect(y_vars)), ", "))),
    ]

    # creating the ode object
    ode_expr = :(StrikeGoldd.ODE{StrikeGoldd.Nemo.fmpq_mpoly}(
        $vx,
        $vy,
        $x_dict,
        $y_dict,
        Array{StrikeGoldd.Nemo.fmpq_mpoly}([$(u_vars...)]),
    ))

    result = Expr(
        :block,
        logging_exprs...,
        exp_ring,
        assignments...,
        x_var_expr,
        y_var_expr,
        x_dict_create_expr,
        y_dict_create_expr,
        eqs_expr...,
        ode_expr,
    )
    return esc(result)
end

#------------------------------------------------------------------------------

function Base.show(io::IO, ode::ODE)
    varstr =
        Dict(x => var_to_str(x) * "(t)" for x in vcat(ode.x_vars, ode.u_vars, ode.y_vars))
    merge!(varstr, Dict(p => var_to_str(p) for p in ode.parameters))
    R_print, vars_print = Nemo.PolynomialRing(
        base_ring(ode.poly_ring),
        [varstr[v] for v in gens(ode.poly_ring)],
    )
    for x in ode.x_vars
        print(io, var_to_str(x) * "'(t) = ")
        print(io, evaluate(ode.x_equations[x], vars_print))
        print(io, "\n")
    end
    for y in ode.y_vars
        print(io, var_to_str(y) * "(t) = ")
        print(io, evaluate(ode.y_equations[y], vars_print))
        print(io, "\n")
    end
end
