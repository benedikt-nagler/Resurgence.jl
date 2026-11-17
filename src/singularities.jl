# Darboux subtraction for a supplied singularity model; not a proof of analytic structure.
# Γ(x) for real x (poles at nonpositive integers excluded), via _gamma and the
# recurrence for x ≤ 0.
function _gamma_real(x::Real)
    x > 0 && return _gamma(x)
    (isinteger(x) || iszero(x)) &&
        throw(InvalidArgument("Γ has a pole at $x"))
    m = Int(ceil(-float(x))) + 1
    g = _gamma(float(x) + m)
    z = float(x) + m
    for _ in 1:m
        z -= 1
        g /= z
    end
    g
end

"""
    subtract_singularity(B::BorelSeries, zeta0; beta = 1, c = nothing,
                         order = nothing, method = :wynn, accel_order = nothing)
        -> NamedTuple{(:peeled, :c, :zeta0, :beta)}
    subtract_singularity(B::BorelSeries; kwargs...)

Darboux peeling: model the dominant singularity of the *reduced* Borel function as
``c\\,(1 - ζ/ζ_0)^{-β}`` (Taylor coefficients ``s_n = c\\,(β)_n ζ_0^{-n}/n!``),
subtract it coefficient-by-coefficient, and return the peeled
[`BorelSeries`](@ref) (same `B.beta`, constant term, and source variable) - so
`poles(pade(out.peeled))` or `large_order_fit(inverse_borel(out.peeled))` can
resolve the **subleading** singularity/action that the dominant one masked.

Here `beta` is the singularity exponent ``β`` - **not** the `B.beta` of the Borel
normalization; `beta = 1` is a simple pole. The amplitude `c` is estimated when not
given: for `beta = 1` from the Padé residue nearest `zeta0`
(``c = -\\mathrm{res}/ζ_0``, exact for a rational Borel function); for general
`beta` by [`accelerate`](@ref)ing the sequence ``b_n ζ_0^n\\, n!/(β)_n → c``
(kwargs `method`, `accel_order`). Without an explicit `zeta0` the
smallest-modulus Borel–Padé pole is used (`order` is passed to [`pade`](@ref));
throws [`NoSingularityFound`](@ref) when there is none.
"""
function subtract_singularity(B::BorelSeries, zeta0::Number;
                              beta::Real = 1,
                              c::Union{Nothing,Number} = nothing,
                              order::Union{Nothing,Integer} = nothing,
                              method::Symbol = :wynn,
                              accel_order::Union{Nothing,Integer} = nothing)
    iszero(zeta0) && throw(InvalidArgument("the singularity position must be ≠ 0"))
    beta > 0 || throw(InvalidArgument("the singularity exponent beta must be > 0"))
    b = B.series.coeffs
    N = length(b)
    camp = c
    if camp === nothing
        if beta == 1
            rs = residues(pade(B; order, reduce = true))
            isempty(rs) &&
                throw(NoSingularityFound("subtract_singularity: the Borel–Padé " *
                                         "approximant has no poles; pass c"))
            _, res = argmin(p -> abs(first(p) - zeta0), rs)
            camp = -res / zeta0
        else
            # c_n = b_n ζ₀ⁿ n!/(β)_n → c, built by the recurrence
            # f_{n+1} = f_n ζ₀ (n+1)/(β+n) on the coefficient-free factor
            f = one(float(b[1])) * one(float(zeta0))
            cn = Vector{typeof(float(b[1]) * f)}(undef, N)
            for n in 0:(N - 1)
                cn[n + 1] = b[n + 1] * f
                f *= zeta0 * (n + 1) / (beta + n)
            end
            camp = accelerate(cn; method, order = accel_order)
        end
    end
    # peel: t_n = c (β)_n ζ₀^{-n}/n!, with t_{n+1} = t_n (β+n)/((n+1) ζ₀)
    t = camp * one(b[1])
    peeled = Vector{typeof(b[1] - t)}(undef, N)
    for n in 0:(N - 1)
        peeled[n + 1] = b[n + 1] - t
        t *= (beta + n) / ((n + 1) * zeta0)
    end
    out = BorelSeries(FormalSeries(peeled, variable(B);
                                   power_offset = B.beta - 1),
                      B.beta, oftype(peeled[1], B.constant_term), B.source_var)
    (peeled = out, c = camp, zeta0 = zeta0, beta = beta)
end

function subtract_singularity(B::BorelSeries;
                              order::Union{Nothing,Integer} = nothing, kwargs...)
    ζs = poles(pade(B; order, reduce = true))
    isempty(ζs) &&
        throw(NoSingularityFound("subtract_singularity: the Borel–Padé approximant " *
                                 "has no poles; pass zeta0 explicitly"))
    subtract_singularity(B, ζs[1]; order, kwargs...)
end
