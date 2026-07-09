# Exact Borel transform. The transform maps ħ^{n+β} ↦ ζ^{n+β-1}/Γ(n+β); we normalize
# by rising factorials (β)_n = Γ(n+β)/Γ(β) instead of Γ(n+β) so the stored
# coefficients stay exact (rational β, any exact coefficient type). The lone 1/Γ(β)
# factor is applied numerically only inside the Laplace step (laplace.jl).

"""
    BorelSeries{T} <: AbstractSeries

The Borel transform of a [`FormalSeries`](@ref), as a distinct type so that a double
Borel transform is a type error. Fields:

- `series::FormalSeries{T}` - the ζ-plane series ``ζ^{β-1} Σ b_n ζ^n`` with the
  *reduced* coefficients ``b_n = a_n / (β)_n`` (rising factorial
  ``(β)_n = β(β+1)⋯(β+n-1)``). The true Borel transform is `1/Γ(β)` times this;
  that constant is applied inside [`laplace_sum`](@ref).
- `beta::Rational{Int}` - the leading ħ-power ``β`` of the source series.
- `constant_term::T` - the ``ħ^0`` coefficient of the source (a constant does not
  Borel-transform; it is split off here and re-added by summation).
- `source_var::Symbol` - the variable of the source series (restored by
  [`inverse_borel`](@ref)).
"""
struct BorelSeries{T} <: AbstractSeries
    series::FormalSeries{T}
    beta::Rational{Int}
    constant_term::T
    source_var::Symbol
end

n_terms(B::BorelSeries) = n_terms(B.series)
coefficients(B::BorelSeries) = coefficients(B.series)
variable(B::BorelSeries) = variable(B.series)
is_exact(B::BorelSeries) = is_exact(B.series)
Base.getindex(B::BorelSeries, n::Integer) = B.series[n]
Base.eltype(::Type{BorelSeries{T}}) where {T} = T

Base.:(==)(A::BorelSeries, B::BorelSeries) =
    A.series == B.series && A.beta == B.beta &&
    A.constant_term == B.constant_term && A.source_var == B.source_var

# (β)_n as exact Rational{BigInt}; poch[n+1] = (β)_n.
function _pochhammer_table(beta::Rational{Int}, n::Integer)
    poch = Vector{Rational{BigInt}}(undef, n)
    poch[1] = 1
    for k in 0:(n - 2)
        poch[k + 2] = poch[k + 1] * (Rational{BigInt}(beta) + k)
    end
    poch
end

"""
    borel(Φ::FormalSeries; beta = power_offset(Φ), var = :ζ) -> BorelSeries

Exact Borel transform of ``Φ = ħ^β Σ a_n ħ^n``:

```math
\\mathcal{B}[Φ](ζ) = \\sum_n \\frac{a_n}{Γ(n+β)} ζ^{n+β-1}
                   = \\frac{ζ^{β-1}}{Γ(β)} \\sum_n \\frac{a_n}{(β)_n} ζ^n .
```

Requires `beta ≥ 0`. For `beta == 0` the constant term ``a_0`` is split off into the
`constant_term` field and the rest is transformed with effective ``β = 1``; the
inverse and the Laplace sum restore it. The transform is exact in every coefficient
type: `inverse_borel(borel(Φ)) == Φ` holds coefficient-by-coefficient.
"""
function borel(Φ::FormalSeries{T};
               beta::Union{Rational{Int},Integer} = power_offset(Φ),
               var::Symbol = :ζ) where {T}
    β = Rational{Int}(beta)
    β ≥ 0 || throw(InvalidArgument("borel requires a leading power β ≥ 0, got $β"))
    if iszero(β)
        constant = Φ.coeffs[1]
        length(Φ.coeffs) ≥ 2 ||
            throw(InvalidArgument("a constant series has no Borel transform; need ≥ 2 terms"))
        a = Φ.coeffs[2:end]
        β = 1//1
    else
        constant = zero(Φ.coeffs[1])
        a = Φ.coeffs
    end
    poch = _pochhammer_table(β, length(a))
    b = [a[k] / poch[k] for k in eachindex(a)]
    # dividing by Rational{BigInt} pochhammers may promote the coefficient type
    BorelSeries(FormalSeries(b, var; power_offset = β - 1), β,
                oftype(b[1], constant), Φ.var)
end

borel(::BorelSeries; kwargs...) =
    throw(InvalidArgument("double Borel transform: the argument is already a BorelSeries"))

"""
    inverse_borel(B::BorelSeries) -> FormalSeries

Exact inverse of [`borel`](@ref): multiplies the reduced coefficients back by the
rising factorials and restores the source variable, power offset, and any split-off
constant term. `inverse_borel(borel(Φ)) == Φ` in every coefficient type.
"""
function inverse_borel(B::BorelSeries{T}) where {T}
    b = B.series.coeffs
    poch = _pochhammer_table(B.beta, length(b))
    a = [b[k] * poch[k] for k in eachindex(b)]
    if iszero(B.constant_term)
        FormalSeries(a, B.source_var; power_offset = B.beta)
    else
        # β was 0 at transform time: constant term split off, effective β = 1
        FormalSeries(vcat(B.constant_term, a), B.source_var; power_offset = 0//1)
    end
end
