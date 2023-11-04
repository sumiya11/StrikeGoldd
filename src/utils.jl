
# ------------------------------------------------------------------------------

function str_to_var(s::String, ring::MPolyRing)
    ind = findfirst(v -> (string(v) == s), symbols(ring))
    if ind === nothing
        throw(Base.KeyError("Variable $s is not found in ring $ring"))
    end
    return gens(ring)[ind]
end

# ------------------------------------------------------------------------------

function var_to_str(v::MPolyElem; xs=gens(parent(v)))
    ind = findfirst(vv -> vv == v, xs)
    return string(symbols(parent(v))[ind])
end

# ------------------------------------------------------------------------------

"""
    parent_ring_change(poly, new_ring)

Converts a polynomial to a different polynomial ring
Input
- `poly` - a polynomial to be converted
- `new_ring` - a polynomial ring such that every variable name appearing in poly appears among the generators

Output:
- a polynomial in `new_ring` “equal” to `poly`
"""
function parent_ring_change(
    poly::MPolyElem,
    new_ring::MPolyRing;
    matching=:byname,
    shift=0
)
    old_ring = parent(poly)
    # Construct a mapping for the variable indices.
    # Zero indicates no image of the old variable in the new ring  
    var_mapping = zeros(Int, max(nvars(old_ring), nvars(new_ring)))
    if matching === :byname
        old_symbols, new_symbols = symbols(old_ring), symbols(new_ring)
        for i in 1:length(old_symbols)
            u = old_symbols[i]
            found = findfirst(v -> (u === v), new_symbols)
            isnothing(found) && continue
            var_mapping[i] = found
        end
    elseif matching === :byindex
        var_mapping[1:(nvars(new_ring)-shift)] .= (1+shift):nvars(new_ring)
    else
        throw(Base.ArgumentError("Unknown matching type: $matching"))
    end
    # Hoist the compatibility check out of the loop
    for i in 1:nvars(old_ring)
        if degree(poly, i) > 0 && iszero(var_mapping[i])
            throw(
                Base.ArgumentError(
                    """
                    The polynomial $poly contains a variable $(gens(old_ring)[i]) not present in the new ring.
                    New ring variables are $(gens(new_ring)))""",
                ),
            )
        end
    end
    bring = base_ring(new_ring)
    exps = Vector{Vector{Int}}(undef, length(poly))
    coefs = map(c -> bring(c), coefficients(poly))
    @inbounds for i in 1:length(poly)
        evec = exponent_vector(poly, i)
        new_exp = zeros(Int, nvars(new_ring))
        for i in 1:length(evec)
            iszero(var_mapping[i]) && continue
            new_exp[var_mapping[i]] = evec[i]
        end
        exps[i] = new_exp
    end
    return new_ring(coefs, exps)
end

function parent_ring_change(
    f::Generic.Frac{<:MPolyElem},
    new_ring::MPolyRing;
    matching=:byname
)
    n, d = unpack_fraction(f)
    return parent_ring_change(n, new_ring; matching=matching) //
           parent_ring_change(d, new_ring; matching=matching)
end

function unpack_fraction(f::MPolyElem)
    return (f, one(parent(f)))
end

function unpack_fraction(f::Generic.Frac{<:MPolyElem})
    return (numerator(f), denominator(f))
end

# ------------------------------------------------------------------------------

"""
    extract_coefficients(poly, variables)

Input:
- `poly` - multivariate polynomial
- `variables` - a list of variables from the generators of the ring of p
Output:
- dictionary with keys being tuples of length `lenght(variables)` and values being polynomials in the variables other than those which are the coefficients at the corresponding monomials (in a smaller polynomial ring)
"""
function extract_coefficients(poly::P, variables) where {P<:MPolyElem}
    xs = gens(parent(poly))
    @assert all(in(xs), variables)
    cut_indices = map(v -> findfirst(x -> x == v, xs), variables)
    coeff_indices = setdiff(collect(1:length(xs)), cut_indices)
    coeff_vars = xs[coeff_indices]

    K = base_ring(parent(poly))
    new_ring, _ = Nemo.PolynomialRing(K, map(vv -> var_to_str(vv, xs=xs), coeff_vars))
    FieldType = elem_type(K)

    result = Dict{Vector{Int},Tuple{Vector{Vector{Int}},Vector{FieldType}}}()

    @inbounds for i in 1:length(poly)
        coef = coeff(poly, i)
        evec = exponent_vector(poly, i)
        var_slice = [evec[i] for i in cut_indices]
        if !haskey(result, var_slice)
            monom_vect, coef_vect = Vector{Vector{Int}}(), Vector{FieldType}()
            sizehint!(monom_vect, 8)
            sizehint!(coef_vect, 8)
            result[var_slice] = (monom_vect, coef_vect)
        end
        monom_vect, coef_vect = result[var_slice]
        new_monom = Vector{Int}(undef, length(coeff_vars))
        for i in 1:length(new_monom)
            new_monom[i] = evec[coeff_indices[i]]
        end
        push!(monom_vect, new_monom)
        push!(coef_vect, coef)
    end

    return Dict(k => new_ring(v[2], v[1]) for (k, v) in result)
end
