# Alien calculus on one-parameter transseries: Stokes-constant extraction from the
# Borel plane (pole residues) and from coefficient asymptotics (large order), the
# bridge equation as a sector map, the Stokes automorphism as the σ-shift, and
# median summation.
#
# Sign conventions (all pinned by the Euler oracle tests, coherent with the M1
# convention disc = s₊ − s₋ = +2πi e^{-1/ħ} for Σ n! ħ^{n+1}):
#   S₁ = 2πi · Res_{ζ=A} B̂(ζ)          (full Borel transform, incl. ζ^{β-1}/Γ(β))
#   s₋(F)(σ) = s₊(F)(σ + S₁)            (Stokes automorphism = parameter shift)
#   disc(Φ₀) = s₊ − s₋ = −S₁ e^{-A/ħ} s(Φ₁)
#   s_med(F)(σ) = s₊(F)(σ + S₁/2)       (real on the ray for real resurgent data)

_first_stokes(stokes::Number) = stokes
function _first_stokes(stokes)
    isempty(stokes) && throw(InvalidArgument("need ≥ 1 Stokes constant, got none"))
    stokes[1]
end

function _nth_stokes(stokes::Number, l::Integer)
    l == 1 || throw(InvalidArgument("a single Stokes constant only defines Δ_{1A}; " *
                                    "pass a vector (S₁, S₂, …) for l = $l"))
    stokes
end
function _nth_stokes(stokes, l::Integer)
    length(stokes) ≥ l ||
        throw(InvalidArgument("need ≥ $l Stokes constants, got $(length(stokes))"))
    stokes[l]
end

"""
    stokes_constant(B::BorelSeries; location = nothing, order = nothing) -> Complex

Stokes constant from the Borel plane: ``S = 2πi · \\mathrm{Res}_{ζ=ζ_0} \\hat{B}(ζ)``
at the Padé pole `ζ₀` nearest `location` (default: the pole closest to the origin,
the leading singularity). The residue of the *full* Borel transform includes the
``ζ^{β-1}/Γ(β)`` prefactor split off by the rising-factorial normalization.

Sign convention (pinned by the Euler oracle): for ``Φ_0 = Σ_k k!\\,ħ^{k+1}`` the
Borel transform is ``1/(1-ζ)`` with residue ``-1`` at ``ζ = A = 1``, so
``S_1 = -2πi``, coherent with [`stokes_discontinuity`](@ref) via
``\\mathrm{disc} = s_+ - s_- = -S_1 e^{-A/ħ} s(Φ_1)``.

Exact (to working precision) for **simple poles**, i.e. when the one-instanton
sector is a constant. At a *branch point* (a nontrivial ``Φ_1``, e.g. Airy) the
Padé pole string emulates the cut and the nearest residue does **not** converge to
``S_1`` - use the large-order method `stokes_constant(Φ::FormalSeries)` there.
"""
function stokes_constant(B::BorelSeries; location::Union{Nothing,Number} = nothing,
                         order::Union{Nothing,Integer} = nothing)
    rs = residues(pade(B; order))
    isempty(rs) && throw(InvalidArgument("the Padé approximant has no poles - " *
                                         "no Borel singularity to read a Stokes constant from"))
    ζ0, ρ = location === nothing ? rs[1] :
            rs[argmin([abs(ζ - location) for (ζ, _) in rs])]
    β = BigFloat(B.beta)
    prefactor = isone(B.beta) ? one(ζ0) : ζ0^(β - 1) / _gamma(β)
    2 * BigFloat(π) * im * prefactor * ρ
end

"""
    stokes_constant(Φ::FormalSeries; order = 4) -> Real

Modulus of the Stokes constant from coefficient asymptotics: the resurgence
large-order relation ``a_n \\sim \\frac{S_1}{2πi}\\, \\frac{Γ(p_n+b)}{A^{p_n+b}}
(c_0^{(1)} + O(1/n))`` gives, for a one-instanton sector normalized to
``c_0^{(1)} = 1``,

```math
|S_1| = 2π\\, S_{\\mathrm{fit}}\\, A^{b} ,
```

with ``(A, b, S_{\\mathrm{fit}})`` from [`large_order_fit`](@ref) (`order` sweeps of
Richardson extrapolation). Valid at branch points where the pole-residue method is
not (the phase of ``S_1`` is not visible in coefficient moduli; it follows from the
direction of the singularity - e.g. purely imaginary for real coefficients).
"""
function stokes_constant(Φ::FormalSeries; order::Integer = 4)
    fit = large_order_fit(Φ; order)
    2 * oftype(fit.S, π) * fit.S * fit.A^fit.b
end

"""
    alien_derivative(F::Transseries, l::Integer = 1; stokes) -> Transseries

The alien derivative ``Δ_{lA}`` of a one-parameter transseries via the **bridge
equation** ``\\dot{Δ}_{lA} F = S_l\\, σ^{1-l}\\, ∂_σ F``, which sector-wise reads

```math
Δ_{lA}\\, Φ_m = S_l\\, (m + l)\\, Φ_{m+l} .
```

`stokes` is the Stokes constant ``S_1`` (a number, for `l = 1`) or a vector
``(S_1, S_2, …)``. The result has `n_sectors(F) - l` sectors; when
`l ≥ n_sectors(F)` every visible sector is annihilated and the zero transseries
(one zero sector) is returned.

# Example
```julia
F = Transseries(:euler, 20)
alien_derivative(F; stokes = -2π*im)     # Δ_A F: single sector S₁·Φ₁ = -2πi
```
"""
function alien_derivative(F::Transseries, l::Integer = 1; stokes)
    l ≥ 1 || throw(InvalidArgument("alien derivative index l must be ≥ 1, got $l"))
    S = _nth_stokes(stokes, l)
    N = n_sectors(F)
    if l ≥ N
        z = FormalSeries([zero(S * F.sectors[1].coeffs[1])], F.var)
        return Transseries(F.action, [z], F.var)
    end
    Transseries(F.action, [S * (m + l) * F.sectors[m + l + 1] for m in 0:(N - 1 - l)],
                F.var)
end

"""
    stokes_automorphism(F::Transseries, σ::Number; stokes) -> Number

The Stokes automorphism ``\\mathfrak{S} = \\exp(\\sum_{l≥1} \\dot{Δ}_{lA})`` acting
on a one-parameter transseries. By the bridge equation it exponentiates to the
**parameter shift** ``σ ↦ σ + S_1``: the two lateral sums of the same transseries
represent the same function at shifted parameters,

```math
s_-(F)(σ) = s_+(F)(σ + S_1) ,
```

which is the connection formula across the Stokes ray (pinned by the Euler oracle
test). This function returns the shifted parameter `σ + S₁`; apply it between
`side = :minus` and `side = :plus` calls of [`transseries_sum`](@ref).
"""
stokes_automorphism(::Transseries, σ::Number; stokes) = σ + _first_stokes(stokes)

"""
    median_sum(F::Transseries, σ, ħ; stokes, θ = 0, kwargs...) -> Complex

Median summation of a transseries: the lateral sum at the half-shifted parameter,
``s_{\\mathrm{med}}(F)(σ) = s_+(F)(σ + S_1/2)`` - halfway (in the sense of
``\\mathfrak{S}^{1/2}``) between the two lateral representatives, so
``s_-(F)(σ - S_1/2)`` gives the same value. For real resurgent data (real
coefficients, real ``σ``, real ``ħ`` on the Stokes ray) the median sum is **real**.
"""
function median_sum(F::Transseries, σ::Number, ħ::Number; stokes, θ::Real = 0,
                    kwargs...)
    transseries_sum(F, σ + _first_stokes(stokes) / 2, ħ; θ, side = :plus, kwargs...)
end

"""
    median_sum(Φ::FormalSeries, ħ; beta = power_offset(Φ), θ = 0, kwargs...) -> Complex

Median summation of a single series: the average of the two lateral sums,
``(s_+(Φ) + s_-(Φ))/2``. This is the leading-sector case of the transseries method
(they agree exactly when the one-instanton sector is constant and higher sectors
vanish, e.g. Euler); real on the ray for real coefficients and real ``ħ``.
"""
function median_sum(Φ::FormalSeries, ħ::Number;
                    beta::Union{Rational{Int},Integer} = power_offset(Φ),
                    θ::Real = 0, kwargs...)
    B = borel(Φ; beta)
    (lateral_sum(B, ħ; θ, side = :plus, kwargs...) +
     lateral_sum(B, ħ; θ, side = :minus, kwargs...)) / 2
end
