# Alien calculus on general K-parameter transseries: alien operators Δ_ω at any lattice
# charge, the pointed (weight-preserving) derivatives Δ̇_ω = e^{-ω/ħ}Δ_ω that generate the
# Stokes automorphism, and the automorphism in its two guises — the summation-time
# per-direction σ-shift (user-facing, reproducing the M2 connection formula) and the
# operator exp(Σ_ω Δ̇_ω) acting on a transseries (bookkeeping). Sign conventions inherit
# from alien.jl (Euler-pinned S₁ = −2πi, s₋(F)(σ) = s₊(F)(σ + S₁)).
#
# The general bridge equation for a single-direction charge ℓ = c·e_j reads
#   Δ_{ℓ·A} Φ_n = S_ℓ · μ_ℓ(n) · Φ_{n+ℓ},   μ_ℓ(n) = n_j + c   (source component j),
# reducing to Δ_{lA}Φ_m = S_l (m+l) Φ_{m+l} at K = 1. Mixed charges (≥ 2 nonzero
# components) have a model-dependent multiplicity and require an explicit `multiplicity`.

_unit(K::Int, j::Int) = ntuple(i -> i == j ? 1 : 0, K)

# default bridge multiplicity for a single-direction charge (result index n)
function _bridge_multiplicity(n::NTuple{K,Int}, ℓ::NTuple{K,Int}) where {K}
    j = findfirst(!iszero, ℓ)
    n[j] + ℓ[j]
end

# resolve the Stokes constant S_ℓ for a charge: scalar (any charge), dict, or callable
_resolve_stokes(s::Number, ℓ) = s
function _resolve_stokes(s::AbstractDict, ℓ)
    haskey(s, ℓ) ||
        throw(InvalidArgument("no Stokes constant provided for alien charge ℓ = $ℓ"))
    s[ℓ]
end
_resolve_stokes(s, ℓ) = s(ℓ)

"""
    alien_derivative(F::MultiTransseries, ℓ = e_1; stokes, multiplicity = nothing)
        -> MultiTransseries

The alien derivative ``Δ_{ℓ·A}`` at integer lattice charge `ℓ` (an `NTuple{K}`), via the
general bridge equation ``Δ_{ℓ·A} Φ_n = S_ℓ\\, μ_ℓ(n)\\, Φ_{n+ℓ}``. Weight-lowering: the
result sector at `n` reads the input sector at `n+ℓ` (dropped when `n ⋡ 0`), so the
result carries the exponential weights shifted down by `ℓ`.

`stokes` supplies ``S_ℓ`` — a scalar (used for every charge), a `Dict` `ℓ ↦ S_ℓ`, or a
callable `ℓ -> S_ℓ`. For a **single-direction** charge `ℓ = c·e_j` the multiplicity
``μ_ℓ(n) = n_j + c`` is used automatically; a **mixed** charge (≥ 2 nonzero components)
is model-dependent and requires an explicit `multiplicity = (n, ℓ) -> Number`.

Passing a `Vector` of charges composes the derivatives ``Δ_{ℓ_r} ∘ ⋯ ∘ Δ_{ℓ_1}`` (first
element applied innermost). Reduces to the one-parameter [`alien_derivative`](@ref) under
the rank-1 embedding.
"""
function alien_derivative(F::MultiTransseries{T,A,K}, ℓ::NTuple{K,Integer};
                          stokes, multiplicity = nothing) where {T,A,K}
    all(iszero, ℓ) && throw(InvalidArgument("alien charge ℓ must be nonzero"))
    ℓi = NTuple{K,Int}(ℓ)
    nz = count(!iszero, ℓi)
    μ = if multiplicity !== nothing
        multiplicity
    elseif nz == 1
        _bridge_multiplicity
    else
        throw(InvalidArgument("mixed-charge alien derivative (ℓ = $ℓi) needs an explicit " *
                              "`multiplicity = (n, ℓ) -> Number`; the bridge multiplicity " *
                              "is model-dependent for charges with ≥ 2 nonzero components"))
    end
    S = _resolve_stokes(stokes, ℓi)
    # src = n + ℓ is the sector read from; keep results with n = src − ℓ ≥ 0. The value
    # type is inferred (S may promote rational sectors to a wider type).
    kept = [(src .- ℓi) => (S * μ(src .- ℓi, ℓi)) * Φ
            for (src, Φ) in F.sectors if all(≥(0), src .- ℓi)]
    d = if isempty(kept)
        z = zero(S * first(values(F.sectors)).coeffs[1])
        Dict(ntuple(_ -> 0, K) => FormalSeries([z], F.var))
    else
        Dict(kept)
    end
    MultiTransseries(F.actions, d; var = F.var)
end

alien_derivative(F::MultiTransseries{T,A,K}; stokes, multiplicity = nothing) where {T,A,K} =
    alien_derivative(F, _unit(K, 1); stokes, multiplicity)

function alien_derivative(F::MultiTransseries, ℓs::AbstractVector; stokes,
                          multiplicity = nothing)
    G = F
    for ℓ in ℓs
        G = alien_derivative(G, ℓ; stokes, multiplicity)
    end
    G
end

"""
    pointed_alien_derivative(F::MultiTransseries, ℓ = e_1; stokes, multiplicity = nothing)
        -> MultiTransseries

The **pointed** alien derivative ``\\dot{Δ}_{ℓ·A} = e^{-(ℓ·A)/ħ}\\, Δ_{ℓ·A}`` —
weight-preserving (the reintroduced weight `e^{-(ℓ·A)/ħ}` cancels the lowering, so a
sector at `n+ℓ` maps back to `n+ℓ`). These are the building blocks of the Stokes
automorphism `exp(Σ_ω Δ̇_ω)` and commute for distinct fundamental directions.
"""
function pointed_alien_derivative(F::MultiTransseries{T,A,K}, ℓ::NTuple{K,Integer};
                                  stokes, multiplicity = nothing) where {T,A,K}
    D = alien_derivative(F, ℓ; stokes, multiplicity)
    ℓi = NTuple{K,Int}(ℓ)
    d = Dict((n .+ ℓi) => Φ for (n, Φ) in D.sectors)   # value type inherited from D
    MultiTransseries(F.actions, d; var = F.var)
end

pointed_alien_derivative(F::MultiTransseries{T,A,K}; stokes,
                         multiplicity = nothing) where {T,A,K} =
    pointed_alien_derivative(F, _unit(K, 1); stokes, multiplicity)

# -- Stokes automorphism -----------------------------------------------------------

# is action a on the ray at angle θ (mod 2π)?
function _on_ray(a::Number, θ::Real)
    φ = mod(angle(complex(float(a))) - θ, 2π)
    isapprox(φ, 0; atol = 1e-10) || isapprox(φ, 2π; atol = 1e-10)
end

"""
    stokes_automorphism(F::MultiTransseries, σ; stokes, θ = 0) -> NTuple

The Stokes automorphism as the **parameter shift** it induces: for each fundamental
direction `j` whose action `A_j` lies on the ray `θ`, `σ_j ↦ σ_j + S_{e_j}`; directions
off the ray are unchanged. The shift `S_{e_j}` is the leading Stokes constant of the
one-instanton sector in direction `j` (equivalently the leading coefficient of
[`pointed_alien_derivative`](@ref)`(F, e_j)`), so the σ-shift is a *derived* consequence
of the alien machinery rather than a definition. At `K = 1` this returns `σ + S₁`,
reproducing the M2 connection formula `s₋(F)(σ) = s₊(F)(σ + S₁)`.
"""
function stokes_automorphism(F::MultiTransseries{T,A,K}, σ::Tuple; stokes,
                             θ::Real = 0) where {T,A,K}
    length(σ) == K ||
        throw(InvalidArgument("parameter vector σ has length $(length(σ)), expected $K"))
    ntuple(K) do j
        _on_ray(F.actions[j], θ) ? σ[j] + _resolve_stokes(stokes, _unit(K, j)) : σ[j]
    end
end

stokes_automorphism(F::MultiTransseries{T,A,1}, σ::Number; kwargs...) where {T,A} =
    stokes_automorphism(F, (σ,); kwargs...)[1]

# one application of the generator 𝒩 = Σ_{on-ray directions j} S_{e_j} Δ_{e_j}
function _stokes_generator(G::MultiTransseries{T,A,K}, stokes, θ::Real) where {T,A,K}
    acc = nothing
    for j in 1:K
        _on_ray(G.actions[j], θ) || continue
        term = alien_derivative(G, _unit(K, j); stokes)
        acc = acc === nothing ? term : acc + term
    end
    acc
end

"""
    stokes_automorphism(F::MultiTransseries; stokes, θ = 0, order = nothing)
        -> MultiTransseries

The Stokes automorphism as an **operator** ``\\mathfrak{S}_θ = \\exp\\big(\\sum_j S_{e_j}
Δ_{e_j·A}\\big) = \\sum_{m≥0} \\tfrac{1}{m!}\\, \\mathcal N^m F`` with generator
``\\mathcal N = \\sum_j S_{e_j} Δ_{e_j·A}`` over the fundamental directions on the ray
`θ`. The plain alien derivatives lower the total weight, so ``\\mathcal N`` is nilpotent
on the bounded support and the series terminates (`order` caps it if given). This is the
closed derivation-algebra object on `MultiTransseries`; the summed action of its
exponential is the parameter shift of the other method.
"""
function stokes_automorphism(F::MultiTransseries{T,A,K}; stokes, θ::Real = 0,
                             order::Union{Nothing,Integer} = nothing) where {T,A,K}
    acc = F
    term = F
    cap = order === nothing ? sum(_bounding_box(F)) + 1 : order
    for m in 1:cap
        term = _stokes_generator(term, stokes, θ)
        (term === nothing || all(all(iszero, Φ.coeffs) for Φ in values(term.sectors))) &&
            break
        acc = acc + (1 // factorial(big(m))) * term
    end
    acc
end

"""
    median_sum(F::MultiTransseries, σ, ħ; stokes, θ = 0, kwargs...) -> Complex

Median summation of a `K`-parameter transseries: the lateral sum at the per-direction
half-shifted parameter ``σ_j + S_{e_j}/2`` on the ray `θ`. Real for real resurgent data
(real coefficients, real `σ`, real `ħ`), as in the one-parameter [`median_sum`](@ref).
"""
function median_sum(F::MultiTransseries{T,A,K}, σ::Tuple, ħ::Number; stokes,
                    θ::Real = 0, kwargs...) where {T,A,K}
    length(σ) == K ||
        throw(InvalidArgument("parameter vector σ has length $(length(σ)), expected $K"))
    shifted = ntuple(K) do j
        _on_ray(F.actions[j], θ) ? σ[j] + _resolve_stokes(stokes, _unit(K, j)) / 2 : σ[j]
    end
    transseries_sum(F, shifted, ħ; θ, side = :plus, kwargs...)
end

median_sum(F::MultiTransseries{T,A,1}, σ::Number, ħ::Number; kwargs...) where {T,A} =
    median_sum(F, (σ,), ħ; kwargs...)
