# Topological-string resurgence: the resolved conifold and the universal conifold gap.
#
# This file is a *consumer*. Every free energy here is a quoted closed form - Bernoulli
# numbers and polylogarithms - not the output of topological recursion, holomorphic-anomaly
# recursion, or a Gromov-Witten / Gopakumar-Vafa computation. The only enumerative input is
# the single resolved-conifold Gopakumar-Vafa invariant n⁰₁ = 1, which enters through the
# quoted multicover resummation [GV98]. Nothing here generates invariants.
#
# Convention ledger
#
#   λ-grid       Free energies are series in the string coupling λ on the *unit* λ grid:
#                power offset 2, odd coefficients exactly zero, coefficient of λ^{2g-2}
#                equal to F_g. This is the normalization in which `borel` divides genus g
#                by Γ(2g-2) = (2g-3)!, which is exactly what makes the gap series' Borel
#                transform elementary. Pinned by `conifold_gap_borel` against
#                `borel(conifold_gap_series(...))`, and by the closed form's independent
#                reconstruction from the sinh series in `test_topological_string.jl`.
#
#   gap sign     `conifold_gap_series` carries NO (-1)^{g-1}:
#                F_g = B_{2g} / (2g(2g-2)) · t^{2-2g}, Borel singularities at ζ = 2πi m t
#                (imaginary for real t). The resolved conifold's own conifold limit is this
#                series at t ↦ i t, which rotates the singularities onto the real axis.
#                Both are first-class; the rotation is pinned by the Q → 1 polylog residue
#                oracle `lim (1-Q)^{n+1} Li_{-n}(Q) = n!`.
#
#   Q vs t       Q = e^{-t}. `resolved_conifold_series` takes the Kähler parameter as `Q`
#                (exact rational coefficients for rational Q); the action lattices take `t`.
#
#   Laurent data `conifold_gap_stokes_data` reports the Borel transform's Laurent
#                coefficients (c₂, c₁) at ζ_m, in the ζ-plane (not the rescaled z-plane).
#                Combined with the package's disc = s₊ - s₋ convention (`laplace.jl`), a
#                pole contributes  disc_m = [2πi c₂/λ - 2πi c₁] e^{-ζ_m/λ}; that sign is
#                pinned numerically against `stokes_discontinuity` in the test file.
#
#   RC lattice   The resolved conifold's Borel singularities are placed at
#                ζ_{m,n} = 2π m (t + 2πi n), m ≠ 0, n ∈ ℤ. This is *derived here*, by
#                pinching the poles of 1/(2 sin(x/2))² at x = 2πm against the poles of the
#                d-sum's geometric factor at x = ζ/(t + 2πi n), and it is verified
#                numerically against `borel_pade_poles`; it is not copied from [ASTT22].

"""
    conifold_gap_series(n::Integer; t = 1//1, var = :λ) -> FormalSeries

The conifold-gap free energy as a formal series in the string coupling `λ`, with `n`
genus terms (genus `2` through `n+1`):

```math
F(λ, t) = \\sum_{g≥2} \\frac{B_{2g}}{2g(2g-2)}\\, t^{2-2g}\\, λ^{2g-2} .
```

This is the universal behaviour of the topological-string free energies of *any*
Calabi-Yau near a conifold point, with `t` the flat coordinate vanishing there
(Huang-Klemm-Quackenbush 2007, "the gap condition"); it is equally the c = 1 string at
self-dual radius, and the leading `t → 0` term of the resolved conifold
([`resolved_conifold_series`](@ref)).

The result lives on the unit `λ` grid with power offset 2 and exactly-zero odd
coefficients, so that [`borel`](@ref) divides genus `g` by `Γ(2g-2) = (2g-3)!`. In that
normalization the Borel transform is *elementary*,

```math
\\mathcal{B}[F](ζ) = \\frac1t\\left[\\frac{-1}{4z\\sinh^2(z/2)} + \\frac1{z^3}
                    - \\frac1{12 z}\\right], \\qquad z = ζ/t ,
```

see [`conifold_gap_borel`](@ref): a tower of *double* poles at ``ζ = 2πi m t``,
``m ≠ 0``. Coefficients are exact `Rational{BigInt}` for rational `t`.

```jldoctest
julia> F = conifold_gap_series(3);

julia> coefficients(F)
5-element Vector{Rational{BigInt}}:
 -1//240
   0
  1//1008
   0
 -1//1440
```
"""
function conifold_gap_series(n::Integer; t = 1//1, var::Symbol = :λ)
    n ≥ 1 || throw(InvalidArgument("need n ≥ 1 genus terms, got $n"))
    tt = t isa Integer ? Rational{BigInt}(t) : t
    iszero(tt) && throw(InvalidArgument("conifold_gap_series needs t ≠ 0"))
    B = _bernoulli(2n + 2)
    first = _gap_coefficient(B, 2, tt)
    c = zeros(typeof(first), 2n - 1)
    for g in 2:(n + 1)
        c[2g - 3] = _gap_coefficient(B, g, tt)
    end
    FormalSeries(c, var; power_offset = 2//1)
end

# F_g = B_{2g} / (2g(2g-2)) · t^{2-2g}, with B[k+1] = B_k.
_gap_coefficient(B, g::Integer, t) = (B[2g + 1] / (2g * (2g - 2))) * t^(2 - 2g)

"""
    resolved_conifold_series(n::Integer; Q, var = :λ) -> FormalSeries

The all-genus Gromov-Witten free energy of the resolved conifold as a formal series in
the string coupling `λ`, with `n` genus terms (genus `2` through `n+1`) and Kähler
parameter `Q = e^{-t}`:

```math
F(λ, t) = \\sum_{g≥2} (-1)^{g-1}\\,\\frac{B_{2g}}{2g\\,(2g-2)!}\\,
          \\mathrm{Li}_{3-2g}(Q)\\; λ^{2g-2} .
```

These `F_g` are *quoted*, not computed: they are the genus expansion of the
Gopakumar-Vafa resummation ``\\sum_{d≥1} Q^d / (d\\,(2\\sin(dλ/2))^2)`` (Gopakumar-Vafa
1998) of the single invariant ``n^0_1 = 1``. For `3-2g ≤ -1` the polylogarithm is a
rational function of `Q`, so the coefficients are exact `Rational{BigInt}` whenever `Q`
is rational.

Layout matches [`conifold_gap_series`](@ref) (unit `λ` grid, power offset 2, zero odd
coefficients); as `t → 0` the coefficients reduce to that gap series at `t ↦ i t`, since
``\\lim_{Q→1}(1-Q)^{n+1}\\mathrm{Li}_{-n}(Q) = n!``. The Borel singularities sit on the
lattice returned by [`resolved_conifold_actions`](@ref), and the Stokes jumps across the
corresponding rays are where Donaldson-Thomas invariants enter (Alim-Saha-Teschner-Tulli
2023, Alim 2024).
"""
function resolved_conifold_series(n::Integer; Q, var::Symbol = :λ)
    n ≥ 1 || throw(InvalidArgument("need n ≥ 1 genus terms, got $n"))
    QQ = Q isa Integer ? Rational{BigInt}(Q) : Q
    isone(QQ) && throw(InvalidArgument("resolved_conifold_series needs Q ≠ 1"))
    B = _bernoulli(2n + 2)
    first = _resolved_coefficient(B, 2, QQ)
    c = zeros(typeof(first), 2n - 1)
    for g in 2:(n + 1)
        c[2g - 3] = _resolved_coefficient(B, g, QQ)
    end
    FormalSeries(c, var; power_offset = 2//1)
end

# F_g = (-1)^{g-1} B_{2g} / (2g(2g-2)!) · Li_{3-2g}(Q)
function _resolved_coefficient(B, g::Integer, Q)
    κ = (isodd(g) ? 1 : -1) * B[2g + 1] / (2g * factorial(big(2g - 2)))
    κ * _polylog_negative(2g - 3, Q)
end

"""
    conifold_gap_borel(ζ; t = 1) -> Number

The exact Borel transform of [`conifold_gap_series`](@ref), in closed form:

```math
\\mathcal{B}[F](ζ) = \\frac1t\\left[\\frac{-1}{4z\\sinh^2(z/2)} + \\frac1{z^3}
                    - \\frac1{12 z}\\right], \\qquad z = ζ/t .
```

The two subtracted terms remove the genus-0 and genus-1 poles at the origin, leaving a
function analytic at `ζ = 0` whose only singularities are *double* poles at
``ζ = 2πi m t``, ``m ≠ 0`` - the exact statement behind the gap condition's large-order
behaviour. For `|ζ/t| < 1` the Bernoulli series is used instead of the closed form, which
there suffers catastrophic cancellation.

This is the oracle for the whole file: `borel(conifold_gap_series(n; t))` must reproduce
its Taylor coefficients exactly, and `laplace_sum` through a Padé approximant must
reproduce its Laplace integral numerically.
"""
function conifold_gap_borel(ζ::Number; t::Number = 1)
    z = ζ / t
    if abs(z) < 1
        T = typeof(float(z))
        z2 = z * z
        acc = zero(T)
        for k in length(_GAP_BOREL_COEFFS):-1:1                # Horner in z²
            acc = acc * z2 + T(_GAP_BOREL_COEFFS[k])
        end
        return acc * z / t                                     # the odd power z^{2g-3}
    end
    (-1 / (4z * sinh(z / 2)^2) + 1 / z^3 - 1 / (12z)) / t
end

"""
    conifold_gap_actions(; t = 1, mmax = 3) -> Vector{Complex}

The Borel singularities ``ζ_m = 2πi m t`` of [`conifold_gap_series`](@ref) for
``0 < |m| ≤ mmax``, sorted by modulus. For real `t` these sit on the imaginary axis, so
the Stokes rays are ``θ = ±π/2``.
"""
function conifold_gap_actions(; t = 1, mmax::Integer = 3)
    mmax ≥ 1 || throw(InvalidArgument("mmax must be ≥ 1, got $mmax"))
    T = typeof(float(real(t)))
    ζ = [2 * T(π) * im * m * t for m in Iterators.flatten((1:mmax, -1:-1:-mmax))]
    sort!(ζ; by = abs)
end

"""
    conifold_gap_stokes_data(m::Integer; t = 1) -> NamedTuple{(:action, :double, :residue)}

The Laurent data of [`conifold_gap_borel`](@ref) at its `m`-th singularity
``ζ_m = 2πi m t``:

```math
\\mathcal{B}[F](ζ) = \\frac{c_2}{(ζ-ζ_m)^2} + \\frac{c_1}{ζ-ζ_m} + O(1),
\\qquad c_2 = \\frac{-t}{2πi m},\\quad c_1 = \\frac{-1}{4π^2m^2} .
```

The double pole is what distinguishes a topological-string free energy from the Euler /
Airy oracles, and the ``1/m^2`` in `c_1` is the multicover weight of the single
underlying state. With the package's `disc = s₊ - s₋` convention ([`stokes_discontinuity`](@ref))
these give the exact discontinuity across the ray through ``ζ_m``,

```math
\\mathrm{disc}_m = \\left[\\frac{2πi\\,c_2}{λ} - 2πi\\,c_1\\right] e^{-ζ_m/λ}
                 = \\left[\\frac{-t}{mλ} + \\frac{i}{2πm^2}\\right] e^{-2πi m t/λ} ,
```

summed over all `m` of the same sign on that ray.
"""
function conifold_gap_stokes_data(m::Integer; t = 1)
    iszero(m) && throw(InvalidArgument("conifold_gap_stokes_data needs m ≠ 0"))
    T = typeof(float(real(t)))
    twoπi = 2 * T(π) * im
    ζm = twoπi * m * t
    (action = ζm, double = -t / (twoπi * m), residue = -1 / (4 * T(π)^2 * m^2))
end

"""
    resolved_conifold_actions(; t = 1, mmax = 2, nmax = 1) -> Vector{Complex}

The Borel singularities of [`resolved_conifold_series`](@ref) at Kähler parameter
`t = -log Q`:

```math
ζ_{m,n} = 2π m\\,(t + 2πi n), \\qquad m ≠ 0,\\ n ∈ ℤ ,
```

for ``0 < |m| ≤ mmax`` and ``|n| ≤ nmax``, sorted by modulus. The `n = 0` tower
``ζ = 2π m t`` is the D2-brane tower whose ``m``-fold covers carry the ``1/m^2``
multicover weight; the ``n ≠ 0`` shifts are the D0-brane charge. For real `t` the
``n = 0`` singularities lie on the positive and negative real axes, so ``θ = 0`` is a
Stokes ray - this is the geometry in which the Stokes constants are read as
Donaldson-Thomas invariants (Alim-Saha-Teschner-Tulli 2023).

The lattice is derived in this package (see the ledger at the top of
`src/topological_string.jl`) and verified numerically against
[`borel_pade_poles`](@ref); it is not read off from the literature.
"""
function resolved_conifold_actions(; t = 1, mmax::Integer = 2, nmax::Integer = 1)
    mmax ≥ 1 || throw(InvalidArgument("mmax must be ≥ 1, got $mmax"))
    nmax ≥ 0 || throw(InvalidArgument("nmax must be ≥ 0, got $nmax"))
    T = typeof(float(real(t)))
    twoπ = 2 * T(π)
    ζ = [twoπ * m * (t + twoπ * im * n)
         for m in Iterators.flatten((1:mmax, -1:-1:-mmax)), n in (-nmax):nmax]
    sort!(vec(ζ); by = abs)
end

# ── Exact special-number helpers ────────────────────────────────────────────────────
# Bernoulli numbers and negative-index polylogarithms, both as exact rationals, so the
# free energies above stay in Rational{BigInt}. Neither is worth a dependency: the
# package's [deps] stay at QuadGK + PrecompileTools.

# B_0 … B_N as Rational{BigInt} with B[k+1] = B_k, from Σ_{k=0}^{m} C(m+1,k) B_k = 0
# (B_1 = -1/2 branch; only the even ones are ever used here).
function _bernoulli(N::Integer)
    N ≥ 0 || throw(InvalidArgument("need N ≥ 0, got $N"))
    B = Vector{Rational{BigInt}}(undef, N + 1)
    B[1] = 1
    for m in 1:N
        acc = zero(Rational{BigInt})
        for k in 0:(m - 1)
            acc += binomial(big(m + 1), big(k)) * B[k + 1]
        end
        B[m + 1] = -acc // (m + 1)
    end
    B
end

# Li_{-n}(Q) = P_n(Q) / (1-Q)^{n+1} for n ≥ 0, where P_n is the Eulerian polynomial
# generated by P_0 = Q and P_n = Q[(1-Q) P_{n-1}' + n P_{n-1}] (the recursion is
# Li_{-n} = (Q d/dQ) Li_{-(n-1)} applied to the rational form). Exact for exact Q.
function _polylog_negative(n::Integer, Q)
    n ≥ 0 || throw(InvalidArgument("_polylog_negative needs n ≥ 0, got $n"))
    P = _eulerian_polynomial(n)
    acc = zero(Q * one(eltype(P)))
    for k in length(P):-1:1
        acc = acc * Q + P[k]
    end
    acc / (1 - Q)^(n + 1)
end

# coefficients of P_n, index i holding the coefficient of Q^{i-1}
function _eulerian_polynomial(n::Integer)
    P = BigInt[0, 1]                                     # P_0 = Q
    for k in 1:n
        d = length(P) - 1
        dP = BigInt[i * P[i + 1] for i in 1:d]           # P'
        # S = (1-Q) P' + k P, then P_k = Q·S
        S = zeros(BigInt, max(length(dP) + 1, length(P)))
        for i in eachindex(dP)
            S[i] += dP[i]
            S[i + 1] -= dP[i]
        end
        for i in eachindex(P)
            S[i] += k * P[i]
        end
        P = pushfirst!(S, zero(BigInt))
        while length(P) > 1 && iszero(P[end])
            pop!(P)
        end
    end
    P
end

# B_{2g}/(2g(2g-2)!) for g = 2 … 60: the Taylor coefficients of the closed-form Borel
# transform at z^{2g-3}. Built once at load - the terms fall off like (2π)^{-2g}, so 59 of
# them carry the |z| < 1 branch of `conifold_gap_borel` past 256-bit working precision.
const _GAP_BOREL_COEFFS = let B = _bernoulli(122)
    [B[2g + 1] / (2g * factorial(big(2g - 2))) for g in 2:60]
end
