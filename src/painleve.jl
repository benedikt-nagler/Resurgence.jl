# Painlevé I as the canonical nonlinear-ODE oracle for the sector-recursion solver
# (`ode.jl`).
#
# Convention ledger - two normalizations of the same equation, both first-class (the
# actions are pinned by `test_painleve.jl`, NOT copied from the roadmap, which had
# A = 8√6/5 against the :gikm equation):
#
#   :gikm   −(1/6) u″ + u² = z         action A = 8√3/5   (A² = 192/25)
#   :string −(1/12) u″ + u² = z        action A = 8√6/5   (A² = 384/25)
#
# [GIKM] (Garoufalidis–Its–Kapaev–Mariño, arXiv:1002.3634) is the oracle spine: it has the
# exactly-rational perturbative recursion (`FormalSeries(:painleve1, n)`) and a documented
# exact Stokes constant S₁ = −i·3^{1/4}/(2√π). :string is the (2,3) minimal-string /
# matrix-model normalization; for −c·u″ + u² = z the action is A = (4/5)√(2/c).
#
# The solver works in the variable s = z^{-1/4}, which turns PI into an ODE with an
# irregular singularity at s = 0 whose transmonomials are e^{-kA s^{-5}} = e^{-kA/ħ}
# (ħ = s⁵ = z^{-5/4}) and whose coefficients are all *integer* powers of s. In s,
#   ∂_z² = (1/16)(s^{10}∂_s² + 5 s⁹ ∂_s),   z = s^{-4},
# and the perturbative solution u₀ = √z Σ aₙ z^{-5n/2} = Σ aₙ s^{-2+10n} lives on the
# unit-s grid (step 10). Sector k's fluctuation ψ_k sits at leading power −2 + 5k/2,
# step 5.

# -- the action --------------------------------------------------------------------

_pi_convention_c(convention::Symbol) =
    convention === :gikm ? 1 // 6 :
    convention === :string ? 1 // 12 :
    throw(InvalidArgument("unknown Painlevé I convention :$convention " *
                          "(expected :gikm or :string)"))

"""
    painleve1_action_squared(; convention = :gikm) -> Rational{Int}

The exact squared instanton action `A²` of Painlevé I: `192//25` in the GIKM
normalization `−(1/6)u″ + u² = z`, `384//25` in the `:string` normalization
`−(1/12)u″ + u² = z`. Rational even though `A` itself is irrational - this is the exact
oracle for [`painleve1_action`](@ref). Derived from the leading balance
`A = (4/5)√(2/c)` for `−c·u″ + u² = z`.
"""
function painleve1_action_squared(; convention::Symbol = :gikm)
    c = _pi_convention_c(convention)
    # A = (4/5)√(2/c) ⇒ A² = (16/25)(2/c) = 32/(25 c)
    (32 // 25) / c
end

"""
    painleve1_action(; convention = :gikm, T = BigFloat) -> T

The instanton action `A = 8√3/5 ≈ 2.7713` (`:gikm`) or `A = 8√6/5 ≈ 3.9192`
(`:string`), as a `T`-valued float. `A` is inexact, so a resonant Painlevé I transseries
must be handed the resonance kernel explicitly - `lattice = [(1, 1)]` for the `(A, −A)`
lattice (see [`resonance_lattice`](@ref)). The exact square is
[`painleve1_action_squared`](@ref).
"""
function painleve1_action(; convention::Symbol = :gikm, T::Type = BigFloat)
    A2 = painleve1_action_squared(; convention)
    sqrt(T(numerator(A2)) / T(denominator(A2)))
end

# -- the perturbative sector on the s-grid, and the ODE residual oracle -------------

# u₀ = Σ aₙ s^{-2+10n} on the unit-s grid, coefficient type T (length ≥ 10(nterms-1)+1).
function _painleve1_u0(nterms::Integer; T::Type = Rational{BigInt}, pad::Integer = 0,
                       convention::Symbol = :gikm)
    a = _painleve1(nterms, _pi_convention_c(convention))
    N = 10 * (nterms - 1) + 1 + pad
    c = zeros(T, N)
    for i in 0:(nterms - 1)
        c[10i + 1] = T(a[i + 1])
    end
    FormalSeries(c, :s; power_offset = -2 // 1)
end

"""
    _painleve1_residual(u0; convention = :gikm) -> FormalSeries

The Painlevé I residual `−c(s^{10}u₀″ + 5 s⁹ u₀′)/16 + u₀² − s^{-4}` on the s-grid, which
vanishes to truncation order for the true perturbative `u₀` (the "every theorem is a
test" oracle for [`FormalSeries(:painleve1, n)`](@ref)). Exact in the coefficient type of
`u0`. Here `c = 1/6` (`:gikm`) or `1/12` (`:string`) and `∂_z² = (s^{10}∂_s² + 5s⁹∂_s)/16`.
"""
function _painleve1_residual(u0::FormalSeries{T}; convention::Symbol = :gikm) where {T}
    c = _pi_convention_c(convention)
    u1 = derivative(u0)
    u2 = derivative(u1)
    d2z = (1 // 16) * (_shift_power(u2, 10) + 5 * _shift_power(u1, 9))   # ∂_z² u₀
    # z = s^{-4}, padded to full length so the subtraction below keeps every order
    zc = zeros(T, n_terms(u0)); zc[1] = one(T)
    one_over_z = FormalSeries(zc, :s; power_offset = -4 // 1)
    (-c) * d2z + u0 * u0 - one_over_z
end

# -- the linearized (Jacobian) operator per instanton sector ------------------------

# L_k ψ = −c·(∂_z²) conjugated by e^{-kA s^{-5}}, plus the 2u₀ from linearizing u²:
#   e^{kA s^{-5}} ∂_z² (e^{-kA s^{-5}} ψ)
#     = (1/16)[ s^{10}ψ″ + (10kA s⁴ + 5 s⁹)ψ′ + (−5kA s³ + 25k²A² s^{-2}) ψ ],
# so L_k ψ = −(c/16)[ … ] + 2 u₀ ψ. The C^k e^{-kA/ħ} grading of −c u″ + u² = z gives
# L_k ψ_k = −Σ_{i=1}^{k-1} ψ_i ψ_{k-i}. The s^{-2} diagonal is 2(1−k²) (nonzero for k≥2);
# for k=1 it vanishes - the transmonomial is exactly balanced - and the free constant is
# the transseries parameter (u_{0,1}=1).
function _painleve1_apply_L(ψ::FormalSeries{T}, k::Integer, A, u0::FormalSeries{T},
                            c::Rational{Int}) where {T}
    ψ1 = derivative(ψ)
    ψ2 = derivative(ψ1)
    t = _shift_power(ψ2, 10)
    t = t + (10 * k * A) * _shift_power(ψ1, 4)
    t = t + 5 * _shift_power(ψ1, 9)
    t = t + (-5 * k * A) * _shift_power(ψ, 3)
    t = t + (25 * k^2 * A^2) * _shift_power(ψ, -2)
    (-c / 16) * t + 2 * (u0 * ψ)
end

# g = −Σ_{i=1}^{k-1} ψ_i ψ_{k-i}, on the s-grid at leading power β_k − 2, length N.
function _painleve1_source(sofar::Vector{FormalSeries{T}}, k::Integer,
                           βk::Rational{Int}, N::Integer) where {T}
    acc = FormalSeries(zeros(T, N), :s; power_offset = βk - 2)
    for i in 1:(k - 1)
        acc = acc - sofar[i + 1] * sofar[k - i + 1]
    end
    acc
end

_painleve1_beta(k::Integer) = -2 + (5 * k) // 2

# one Newton step: solve L_k ψ_k = source_k for the k-th fluctuation
function _painleve1_solve_sector(sofar::Vector{FormalSeries{T}}, k::Integer,
                                 A, u0::FormalSeries{T}, c::Rational{Int},
                                 N::Integer) where {T}
    βk = _painleve1_beta(k)
    g = _painleve1_source(sofar, k, βk, N)
    δ = k == 1 ? 3 : -2                 # diagonal shift: differential (k=1) vs algebraic
    free = k == 1 ? [0] : Int[]         # k=1 carries the free transseries parameter
    apply_L = ψ -> _painleve1_apply_L(ψ, k, A, u0, c)
    _solve_linear_sector(apply_L, g, βk, δ, N; free, freeval = one(T))
end

# -- the assembled transseries ------------------------------------------------------

"""
    painleve1(n; sectors = 2, convention = :gikm, T = BigFloat) -> Transseries

The Painlevé I transseries solution, sectors solved by the [`transseries_solve`](@ref)
recursion in the variable `s = z^{-1/4}`: sector 0 is the perturbative `u₀` and sector
`k ≥ 1` is the one-Newton-step fluctuation of the transmonomial `e^{-kA s^{-5}}`
(`= e^{-kA z^{5/4}}`). Each sector carries `n` fluctuation orders. The action is
`A = ` [`painleve1_action`](@ref)`(; convention)`; the coefficient type `T` defaults to
`BigFloat`.

The returned [`Transseries`](@ref) is in the `s` variable (so its `action` multiplies
`s^{-5} = z^{5/4}`), with sector 0 the exact perturbative series cast to `T`. The
resonant two-parameter `(A, −A)` structure - where `s`-independent logs appear on the
diagonal - is built with `lattice = [(1, 1)]`; see the tests and
[`resonance_lattice`](@ref).
"""
function painleve1(n::Integer; sectors::Integer = 2, convention::Symbol = :gikm,
                   T::Type = BigFloat)
    n ≥ 1 || throw(InvalidArgument("need n ≥ 1 fluctuation orders, got $n"))
    sectors ≥ 1 || throw(InvalidArgument("need sectors ≥ 1, got $sectors"))
    c = _pi_convention_c(convention)
    A = painleve1_action(; convention, T = BigFloat)
    N = 10 * n + 15                     # unit-s grid: u₀ has step 10, sectors step 5
    u0 = _painleve1_u0(n; T = T, pad = N - (10 * (n - 1) + 1), convention)
    solve = (sofar, k) -> _painleve1_solve_sector(sofar, k, A, u0, c, N)
    transseries_solve(u0, A, solve; sectors, var = :s)
end
