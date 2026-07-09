# Superasymptotics and level-1 hyperasymptotics: optimal truncation at the smallest
# term (error ~ e^{-|A|/|ħ|}), Dingle terminant functions, and the Berry–Howls
# level-1 re-expansion of the remainder through the adjacent sector.

"""
    dingle_terminant(ν, x; side = :plus, rtol = eps^(3/4)) -> Complex

The (basic, Dingle-normalized) terminant — the Borel sum of ``Σ_{k≥0} Γ(ν+k)x^k``:

```math
Λ_ν(x) = \\int_0^∞ t^{ν-1} e^{-t}\\, \\frac{dt}{1 - x t} , \\qquad ν > 0 ,
```

evaluated by quadrature in the precision of `x`. ``Λ_ν(0) = Γ(ν)``. For `x` on the
positive real axis the pole ``t = 1/x`` sits on the contour (the Stokes case): the
ray is tilted by `side` (`:plus` above, `:minus` below — lateral terminants), and
the two sides differ by the residue jump
``Λ_ν^+(x) - Λ_ν^-(x) = 2πi\\, x^{-ν} e^{-1/x}``.
"""
function dingle_terminant(ν::Real, x::Number; side::Symbol = :plus,
                          rtol::Union{Nothing,Real} = nothing)
    ν > 0 || throw(InvalidArgument("dingle_terminant requires ν > 0, got $ν"))
    side in (:plus, :minus) ||
        throw(InvalidArgument("side must be :plus or :minus, got :$side"))
    T = _real_float_type(x)
    rt = rtol === nothing ? eps(T)^(3//4) : rtol
    νT = T(ν)
    δ = zero(T)
    if !iszero(x) && abs(angle(x)) < sqrt(eps(T))
        δ = side === :plus ? T(1//100) : -T(1//100)
    end
    eiδ = cis(δ)
    integrand = u -> (u * eiδ)^(νT - 1) * exp(-u * eiδ) / (1 - x * u * eiδ) * eiδ
    val, _ = quadgk(integrand, zero(T), T(Inf); rtol = rt)
    val
end

"""
    optimal_truncation(Φ::FormalSeries, ħ; rule = :smallest_term, A = nothing)
        -> NamedTuple{(:N, :value, :error)}

Superasymptotics: truncate the divergent series at its smallest term. Returns the
number of kept terms ``N^*`` (so `value` ``= Σ_{n<N^*} a_n ħ^{p_n}``) and the
magnitude of the first omitted term as `error` — the classic
``\\sim e^{-|A|/|ħ|}`` superasymptotic error estimate. Rules:

- `:smallest_term` (default) — ``N^* = \\arg\\min_n |a_n ħ^{p_n}|`` over the stored
  range (if the terms are still shrinking at the truncation order, the stored
  range is simply exhausted);
- `:action` — ``N^* = \\mathrm{round}(|A/ħ|) - \\mathrm{power\\_offset}`` clamped
  to the stored range, with the action modulus `A` taken from the keyword or from
  [`large_order_fit`](@ref).

Quadrature-free.
"""
function optimal_truncation(Φ::FormalSeries, ħ::Number; rule::Symbol = :smallest_term,
                            A::Union{Nothing,Number} = nothing)
    N = n_terms(Φ)
    terms = abs.(diff(vcat(zero(first(partial_sums(Φ, ħ))), partial_sums(Φ, ħ))))
    if rule === :smallest_term
        Nstar = argmin(terms) - 1                      # 0-based index of omitted term
    elseif rule === :action
        Amod = abs(A === nothing ? large_order_fit(Φ).A : A)
        Nstar = clamp(round(Int, Amod / abs(ħ) - float(power_offset(Φ))), 0, N - 1)
    else
        throw(InvalidArgument("rule must be :smallest_term or :action, got :$rule"))
    end
    s = partial_sums(Φ, ħ)
    value = Nstar == 0 ? zero(s[1]) : s[Nstar]
    (N = Nstar, value = value, error = terms[Nstar + 1])
end

"""
    hyper_sum(Φ::FormalSeries, ħ; adjacent, action, stokes, N = nothing,
              M = nothing, side = :plus, rtol = eps^(3/4)) -> Complex

Level-1 hyperasymptotics (Berry–Howls): truncate `Φ` at ``N_0`` terms and re-expand
the remainder through the **adjacent sector** ``Φ_1 = Σ_m a^{(1)}_m ħ^{m+β_1}`` at
the dominant `action` ``A`` with Stokes constant `stokes` ``= S_1`` (this package's
convention: ``S_1 = 2πi\\,\\mathrm{Res}`` for a simple Borel pole, see
`stokes_constant`), using [`dingle_terminant`](@ref)s:

```math
φ(ħ) ≈ \\sum_{n<N_0} a_n ħ^{p_n} \\;-\\; \\frac{S_1}{2πi}\\,(ħ/A)^{p_{N_0}}
       \\sum_{m<M} a^{(1)}_m\\, A^{m+β_1}\\, Λ_{p_{N_0}-m-β_1}(ħ/A) .
```

Defaults: ``N_0`` from `optimal_truncation(Φ, ħ; rule = :action, A = action)`; the
terminant series is truncated at its own smallest term (or at `M`, or when the
terminant order ``ν`` would leave the ``ν > 0`` domain). For Euler-type simple
poles the adjacent sector is a constant and the re-expansion is *exact*; the error
of the level-1 result improves on the superasymptotic ``e^{-|A|/|ħ|}`` by the
documented Berry–Howls factor. `side` selects the lateral terminant on Stokes rays
(``ħ/A`` positive real).

Pass plain series data — extract sectors from a `Transseries` yourself
(`sector(F, 1)`, `action(F)`).
"""
function hyper_sum(Φ::FormalSeries, ħ::Number; adjacent::FormalSeries,
                   action::Number, stokes::Number,
                   N::Union{Nothing,Integer} = nothing,
                   M::Union{Nothing,Integer} = nothing,
                   side::Symbol = :plus,
                   rtol::Union{Nothing,Real} = nothing)
    T = _real_float_type(ħ)
    N0 = N === nothing ? optimal_truncation(Φ, ħ; rule = :action, A = action).N : Int(N)
    0 ≤ N0 ≤ n_terms(Φ) ||
        throw(InvalidArgument("truncation N = $N0 outside the stored range " *
                              "0:$(n_terms(Φ))"))
    ħc = Complex{T}(ħ)
    base = N0 == 0 ? zero(ħc) : partial_sums(Φ, ħc)[N0]
    β1 = power_offset(adjacent)
    pN0 = N0 + power_offset(Φ)
    x = ħc / Complex{T}(action)
    pref = -Complex{T}(stokes) / (2 * T(π) * im)
    Apow = Complex{T}(action)^T(β1)                       # A^{m+β₁}
    Mmax = M === nothing ? n_terms(adjacent) : min(Int(M), n_terms(adjacent))
    tail = zero(ħc)
    prevmag = T(Inf)
    for m in 0:(Mmax - 1)
        ν = T(pN0 - m - β1)
        ν > 0 || break
        term = Complex{T}(adjacent[m]) * Apow * dingle_terminant(ν, x; side, rtol)
        mag = abs(term)
        # default truncation of the level-1 series: its own smallest term
        M === nothing && m ≥ 1 && mag > prevmag && break
        tail += term
        prevmag = mag
        Apow *= action
    end
    base + pref * x^T(pN0) * tail
end
