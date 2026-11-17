# One-parameter transseries F(σ, ħ) = Σ_n σ^n e^{-nA/ħ} Φ_n(ħ) with dense integer
# exponential weights: sectors[k] is Φ_{k-1}, and any ħ^{β_n} prefactor lives in that
# sector's power_offset. The parameter σ is *not* stored - a Transseries is the
# family {Φ_n}; σ enters only at summation time (transseries_sum), which is where
# the Stokes automorphism acts (alien.jl).

"""
    Transseries{T,A} <: AbstractSeries

A one-parameter transseries ``F(σ, ħ) = \\sum_{n≥0} σ^n e^{-nA/ħ}\\, Φ_n(ħ)`` with
action `action::A` and sectors `sectors::Vector{FormalSeries{T}}` (dense in the
integer exponential weight: `sectors[k]` is ``Φ_{k-1}``). The prefactor ``ħ^{β_n}``
of each sector lives in its `power_offset`. The transseries parameter ``σ`` is not
stored; it is supplied to [`transseries_sum`](@ref), and the Stokes automorphism
acts on it (see [`stokes_automorphism`](@ref)).

# Examples
```julia
F = Transseries(:euler, 20)                       # named course-canon transseries
F = Transseries(1//1, [Φ₀, Φ₁])                   # from explicit sectors
sector(F, 1)                                      # Φ₁
```
"""
struct Transseries{T,A<:Number} <: AbstractSeries
    action::A
    sectors::Vector{FormalSeries{T}}
    var::Symbol

    function Transseries{T,A}(action::A, sectors::Vector{FormalSeries{T}},
                              var::Symbol) where {T,A<:Number}
        iszero(action) &&
            throw(InvalidArgument("a Transseries needs a nonzero action, got A = 0"))
        isempty(sectors) &&
            throw(InvalidArgument("a Transseries needs ≥ 1 sector"))
        for Φ in sectors
            variable(Φ) == var || throw(IncompatibleSeries(:var, var, variable(Φ)))
        end
        new{T,A}(action, copy(sectors), var)
    end
end

function Transseries(action::A, sectors::Vector{FormalSeries{T}},
                     var::Symbol = isempty(sectors) ? :ħ : variable(sectors[1])) where {T,A<:Number}
    Transseries{T,A}(action, sectors, var)
end

# -- accessors ---------------------------------------------------------------------

"""
    action(F::Transseries) -> Number

The action ``A`` of the exponential weights ``e^{-nA/ħ}``. Its sign convention: `A`
is the location of the leading Borel-plane singularity of ``Φ_0`` (so the weights
``e^{-nA/ħ}`` are exponentially small in the directions where the perturbative
sector dominates).
"""
action(F::Transseries) = F.action

"""
    n_sectors(F::Transseries) -> Int

Number of stored sectors (the highest stored exponential weight is
`n_sectors(F) - 1`).
"""
n_sectors(F::Transseries) = length(F.sectors)

"""
    sector(F::Transseries, n::Integer) -> FormalSeries

The sector ``Φ_n`` of weight ``e^{-nA/ħ}`` (0-based; `n ≥ n_sectors(F)` returns the
zero series, mirroring `Φ[n]` on [`FormalSeries`](@ref)).
"""
function sector(F::Transseries, n::Integer)
    n < 0 && throw(InvalidArgument("sector index must be ≥ 0, got $n"))
    n < length(F.sectors) ? F.sectors[n + 1] : _zero_sector(F)
end

_zero_sector(F::Transseries) =
    FormalSeries([zero(F.sectors[1].coeffs[1])], F.var)

"""
    sectors(F::Transseries) -> Vector{FormalSeries}

The sector list ``[Φ_0, Φ_1, …]`` (a copy; transseries are immutable).
"""
sectors(F::Transseries) = copy(F.sectors)

variable(F::Transseries) = F.var
is_exact(F::Transseries) = all(is_exact, F.sectors)
Base.eltype(::Type{Transseries{T,A}}) where {T,A} = T

Base.:(==)(F::Transseries, G::Transseries) =
    F.action == G.action && F.var == G.var && F.sectors == G.sectors

# -- arithmetic --------------------------------------------------------------------

function _check_compatible(F::Transseries, G::Transseries)
    F.var == G.var || throw(IncompatibleSeries(:var, F.var, G.var))
    F.action == G.action || throw(IncompatibleSeries(:action, F.action, G.action))
    nothing
end

function Base.:+(F::Transseries{T}, G::Transseries{T}) where {T}
    _check_compatible(F, G)
    nF, nG = length(F.sectors), length(G.sectors)
    secs = [i ≤ nF && i ≤ nG ? F.sectors[i] + G.sectors[i] :
            (i ≤ nF ? F.sectors[i] : G.sectors[i]) for i in 1:max(nF, nG)]
    Transseries(F.action, secs, F.var)
end

Base.:-(F::Transseries) = Transseries(F.action, [-Φ for Φ in F.sectors], F.var)
Base.:-(F::Transseries{T}, G::Transseries{T}) where {T} = F + (-G)

"""
    F * G

Weight-graded Cauchy product ``(F·G)_n = \\sum_{i=0}^{n} Φ_i Ψ_{n-i}`` (the
exponential weights multiply as ``e^{-iA/ħ} e^{-(n-i)A/ħ} = e^{-nA/ħ}``), truncated
to `min(n_sectors(F), n_sectors(G))` sectors - the highest weight at which every
contribution is available.
"""
function Base.:*(F::Transseries{T}, G::Transseries{T}) where {T}
    _check_compatible(F, G)
    N = min(length(F.sectors), length(G.sectors))
    secs = [sum(F.sectors[i + 1] * G.sectors[n - i + 1] for i in 0:n) for n in 0:(N - 1)]
    Transseries(F.action, secs, F.var)
end

Base.:*(λ::Number, F::Transseries) =
    Transseries(F.action, [λ * Φ for Φ in F.sectors], F.var)
Base.:*(F::Transseries, λ::Number) = λ * F

"""
    truncate(F::Transseries; sectors = n_sectors(F), terms = nothing) -> Transseries

The transseries truncated to its first `sectors` sectors, each truncated to `terms`
coefficients (`nothing` keeps every sector's own length).
"""
function Base.truncate(F::Transseries; sectors::Integer = n_sectors(F),
                       terms::Union{Nothing,Integer} = nothing)
    1 ≤ sectors ≤ length(F.sectors) ||
        throw(InvalidArgument("sector count must be in 1:$(length(F.sectors)), got $sectors"))
    secs = F.sectors[1:sectors]
    terms === nothing || (secs = [truncate(Φ, terms) for Φ in secs])
    Transseries(F.action, secs, F.var)
end

# -- summation ---------------------------------------------------------------------

"""
    transseries_sum(F::Transseries, σ, ħ; θ = 0, side = :plus, tilt = 1//100,
                    order = nothing, rtol = eps^(3/4)) -> Complex

Lateral Borel sum of the transseries at parameter `σ`:
``\\sum_n σ^n e^{-nA/ħ}\\, s_θ^{±}(Φ_n)(ħ)``, where each sector is summed by
[`lateral_sum`](@ref) of its Borel transform along the ray `θ` (constant sectors,
e.g. a pure one-instanton term ``Φ_1 = 1``, are evaluated exactly instead).

The two lateral sums are related by the Stokes automorphism acting on the
parameter: `transseries_sum(F, σ, ħ; side = :minus)` equals
`transseries_sum(F, σ + S₁, ħ; side = :plus)` - see [`stokes_automorphism`](@ref).
"""
function transseries_sum(F::Transseries, σ::Number, ħ::Number; θ::Real = 0,
                         side::Symbol = :plus, kwargs...)
    T = _real_float_type(ħ)
    total = zero(Complex{T}) * (one(σ) * one(F.action))
    for n in 0:(length(F.sectors) - 1)
        Φ = F.sectors[n + 1]
        all(iszero, Φ.coeffs) && continue
        weight = σ^n * exp(-(n * F.action) / ħ)
        iszero(weight) && continue
        val = n_terms(Φ) == 1 ? evaluate(Φ, ħ) :
              lateral_sum(borel(Φ), ħ; θ, side, kwargs...)
        total += weight * val
    end
    total
end
