# Generic truncated formal power series x^{power_offset} Σ a_n x^n with fully generic
# coefficient type T — Rational{BigInt}, Complex{BigFloat}, and AbstractAlgebra ring
# elements all work by duck typing (no AbstractAlgebra dependency).

"""
    AbstractSeries

Abstract supertype of all series-like objects in Resurgence.jl
([`FormalSeries`](@ref), [`BorelSeries`](@ref)).
"""
abstract type AbstractSeries end

"""
    FormalSeries{T} <: AbstractSeries

A truncated formal power series ``x^{β_0} \\sum_{n=0}^{N} a_n x^n`` in the variable
`var` (default `:ħ`), with coefficients `coeffs::Vector{T}` (so `coeffs[k]` is
``a_{k-1}``) and rational leading power `power_offset` ``= β_0``.

The coefficient type `T` is fully generic: exact rationals, `BigFloat`s, complex
numbers, and ring elements from computer-algebra systems all work.

# Examples
```julia
Φ = FormalSeries([1//1, 1, 2, 6], :ħ)            # Σ n! ħ^n, 4 terms
Φ = FormalSeries([1.0, 1.0], :x; power_offset = 1//2)   # √x (1 + x)
Φ[2]                                              # coefficient a₂ == 2
```
"""
struct FormalSeries{T} <: AbstractSeries
    coeffs::Vector{T}
    var::Symbol
    power_offset::Rational{Int}

    function FormalSeries{T}(coeffs::Vector{T}, var::Symbol,
                             power_offset::Rational{Int}) where {T}
        isempty(coeffs) && throw(InvalidArgument("a FormalSeries needs ≥ 1 coefficient"))
        new{T}(copy(coeffs), var, power_offset)
    end
end

function FormalSeries(coeffs::Vector{T}, var::Symbol = :ħ;
                      power_offset::Union{Rational{Int},Integer} = 0//1) where {T}
    FormalSeries{T}(coeffs, var, Rational{Int}(power_offset))
end

# -- accessors ---------------------------------------------------------------------

"""
    n_terms(Φ::AbstractSeries) -> Int

Number of stored coefficients (the truncation order is `n_terms(Φ) - 1`).
"""
n_terms(Φ::FormalSeries) = length(Φ.coeffs)

"""
    coefficients(Φ::AbstractSeries) -> Vector

The coefficient vector ``[a_0, a_1, …]`` (a copy; series are immutable).
"""
coefficients(Φ::FormalSeries) = copy(Φ.coeffs)

"""
    variable(Φ::AbstractSeries) -> Symbol

The variable the series is written in.
"""
variable(Φ::FormalSeries) = Φ.var

"""
    power_offset(Φ::FormalSeries) -> Rational{Int}

The rational leading power ``β_0`` in ``x^{β_0} Σ a_n x^n``.
"""
power_offset(Φ::FormalSeries) = Φ.power_offset

"""
    is_exact(Φ::AbstractSeries) -> Bool

`true` when the coefficient type supports exact arithmetic (integers, rationals, and
non-`AbstractFloat` number types such as computer-algebra ring elements).
"""
is_exact(Φ::FormalSeries) = _is_exact_type(eltype(Φ.coeffs))

_is_exact_type(::Type{<:Union{Integer,Rational}}) = true
_is_exact_type(::Type{Complex{S}}) where {S} = _is_exact_type(S)
_is_exact_type(::Type{<:AbstractFloat}) = false
_is_exact_type(::Type) = true   # duck-typed ring elements are presumed exact

Base.eltype(::Type{FormalSeries{T}}) where {T} = T

"""
    Φ[n]

Coefficient ``a_n`` of ``x^{β_0 + n}`` (0-based; `n ≥ n_terms(Φ)` returns zero).
"""
function Base.getindex(Φ::FormalSeries{T}, n::Integer) where {T}
    n < 0 && throw(InvalidArgument("coefficient index must be ≥ 0, got $n"))
    n < length(Φ.coeffs) ? Φ.coeffs[n + 1] : zero(Φ.coeffs[1])
end

Base.:(==)(Φ::FormalSeries, Ψ::FormalSeries) =
    Φ.var == Ψ.var && Φ.power_offset == Ψ.power_offset && Φ.coeffs == Ψ.coeffs

# -- arithmetic --------------------------------------------------------------------

function _check_compatible(Φ::FormalSeries, Ψ::FormalSeries)
    Φ.var == Ψ.var || throw(IncompatibleSeries(:var, Φ.var, Ψ.var))
    nothing
end

# Align two series to a common power offset (must differ by an integer for +/−).
function _aligned_coeffs(Φ::FormalSeries{T}, Ψ::FormalSeries{T}) where {T}
    shift = Φ.power_offset - Ψ.power_offset
    denominator(shift) == 1 ||
        throw(IncompatibleSeries(:power_offset, Φ.power_offset, Ψ.power_offset))
    s = Int(numerator(shift))
    # pad the series with the larger offset by |s| leading zeros
    if s ≥ 0
        offset = Ψ.power_offset
        a = vcat(fill(zero(Φ.coeffs[1]), s), Φ.coeffs)
        b = Ψ.coeffs
    else
        offset = Φ.power_offset
        a = Φ.coeffs
        b = vcat(fill(zero(Ψ.coeffs[1]), -s), Ψ.coeffs)
    end
    # a common truncation order: the shorter tail wins
    N = min(length(a), length(b))
    a[1:N], b[1:N], offset
end

function Base.:+(Φ::FormalSeries{T}, Ψ::FormalSeries{T}) where {T}
    _check_compatible(Φ, Ψ)
    a, b, offset = _aligned_coeffs(Φ, Ψ)
    FormalSeries(a .+ b, Φ.var; power_offset = offset)
end

function Base.:-(Φ::FormalSeries{T}, Ψ::FormalSeries{T}) where {T}
    _check_compatible(Φ, Ψ)
    a, b, offset = _aligned_coeffs(Φ, Ψ)
    FormalSeries(a .- b, Φ.var; power_offset = offset)
end

Base.:-(Φ::FormalSeries) = FormalSeries(-Φ.coeffs, Φ.var; power_offset = Φ.power_offset)

"""
    Φ * Ψ

Cauchy product, truncated to `min(n_terms(Φ), n_terms(Ψ))` terms; power offsets add.
"""
function Base.:*(Φ::FormalSeries{T}, Ψ::FormalSeries{T}) where {T}
    _check_compatible(Φ, Ψ)
    N = min(length(Φ.coeffs), length(Ψ.coeffs))
    c = [sum(Φ.coeffs[k + 1] * Ψ.coeffs[n - k + 1] for k in 0:n) for n in 0:(N - 1)]
    FormalSeries(c, Φ.var; power_offset = Φ.power_offset + Ψ.power_offset)
end

Base.:*(λ::Number, Φ::FormalSeries) =
    FormalSeries(λ .* Φ.coeffs, Φ.var; power_offset = Φ.power_offset)
Base.:*(Φ::FormalSeries, λ::Number) = λ * Φ

"""
    truncate(Φ::FormalSeries, n) -> FormalSeries

The series truncated to its first `n` coefficients.
"""
function Base.truncate(Φ::FormalSeries, n::Integer)
    1 ≤ n ≤ length(Φ.coeffs) ||
        throw(InvalidArgument("truncation length must be in 1:$(length(Φ.coeffs)), got $n"))
    FormalSeries(Φ.coeffs[1:n], Φ.var; power_offset = Φ.power_offset)
end

"""
    derivative(Φ::FormalSeries) -> FormalSeries

The termwise derivative ``∂_x``: ``x^{β_0} \\sum a_n x^n \\mapsto x^{β_0 - 1}
\\sum (β_0 + n) a_n x^n``. Exact in the coefficient type whenever the powers are
integral. The companion operator ``x^2 ∂_x`` — the one that governs transseries sector
equations, see [`resonant_solve`](@ref) — is `_shift_power(derivative(Φ), 2)`.
"""
function derivative(Φ::FormalSeries)
    β = Φ.power_offset
    c = [_power_factor(β + (k - 1)) * Φ.coeffs[k] for k in eachindex(Φ.coeffs)]
    FormalSeries(c, Φ.var; power_offset = β - 1)
end

# an integral power multiplies exactly in T; a fractional one has to go through Rational
_power_factor(p::Rational{Int}) = denominator(p) == 1 ? numerator(p) : p

# multiply by x^s (s an integer): shift the leading power, keep the coefficients
_shift_power(Φ::FormalSeries, s::Integer) =
    FormalSeries(Φ.coeffs, Φ.var; power_offset = Φ.power_offset + s)

# the coefficient of x^p, honoring power_offset (zero when the power grid misses p)
function _coeff_at(Φ::FormalSeries{T}, p::Rational{Int}) where {T}
    m = p - Φ.power_offset
    (denominator(m) == 1 && 0 ≤ numerator(m)) || return zero(Φ.coeffs[1])
    Φ[Int(numerator(m))]
end

_is_zero_series(Φ::FormalSeries) = all(iszero, Φ.coeffs)

"""
    evaluate(Φ::FormalSeries, x) -> Number

Evaluate the truncated series at `x` (Horner on the polynomial part, times
`x^power_offset`). This is *formal* evaluation of the truncation — for the resummed
value of a divergent series use [`borel_sum`](@ref).
"""
function evaluate(Φ::FormalSeries, x)
    p = foldr((c, acc) -> c + x * acc, Φ.coeffs; init = zero(Φ.coeffs[1]) * zero(x))
    iszero(Φ.power_offset) ? p : x^Φ.power_offset * p
end
