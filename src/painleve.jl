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
# Two further choices are ours, both in the two-parameter resonant layer below:
#
#   log variable    the tower is in log s, the *series* variable, not log ħ or log z
#                   (log ħ = 5 log s, log z = −4 log s). Pinned by the exact top-log
#                   coefficient of sector (2,1), κ = −10/√3, in `test_painleve.jl`;
#                   a different log rescales κ by the same factor.
#   large-order 2   the factor 2 in `painleve1_stokes_constant`'s fit: Φ₀ has only even
#                   powers of ħ, so both Borel singularities ±A contribute equally at
#                   even order. Pinned by the GIKM value S₁ = −i·3^{1/4}/(2√π); without
#                   it the fit returns 2S₁.
#   Stokes sign     `painleve1_stokes_constant` returns the *GIKM* constant, the one
#                   quoted in the literature. alien.jl's own convention is
#                   disc(Φ₀) = −S₁ e^{-A/ħ} s(Φ₁), so the two differ by a sign:
#                   S₁_alien = −S₁_GIKM = +i·3^{1/4}/(2√π). Feeding this function's
#                   output straight into `alien_derivative`/`stokes_automorphism` needs
#                   that flip. Pinned by the `:lateral` route, which measures
#                   disc/(e^{-A/ħ}Φ₁) directly and lands on +S₁_GIKM.
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

# -- the two-parameter (resonant) transseries ---------------------------------------

# The action lattice (A, −A) is resonant: (n, m) and (n+1, m+1) share the weight (n−m)A,
# and the sector operator 𝓛_k depends on the multi-index only through k = n − m. Writing
# L = log s, the log commutators of 𝓛_k (see `_solve_log_sector`) are, from
# ∂_s(w L^p) = w′L^p + p w L^{p-1}/s,
#   C(w) = −(c/16)[2 s⁹ w′ + 4 s⁸ w + 10kA s³ w],     D(w) = −(c/16) s⁸ w.
# Only k = ±1 has a kernel (the diagonal 2(1−k²) vanishes and the balance turns
# differential), at the power s^{1/2}, and the homogeneous solution there is *the*
# one-instanton fluctuation Φ_{(1,0)} - so the log tower is generated by the sector that
# defines the transseries parameter.

const _PI_KERNEL_POWER = 1 // 2

_painleve1_apply_C(ψ::FormalSeries, k::Integer, A, c::Rational{Int}) =
    (-c / 16) * (2 * _shift_power(derivative(ψ), 9) + 4 * _shift_power(ψ, 8) +
                 (10 * k * A) * _shift_power(ψ, 3))

_painleve1_apply_D(ψ::FormalSeries, c::Rational{Int}) = (-c / 16) * _shift_power(ψ, 8)

# −Σ' Φ_j Φ_{idx−j} over 0 < j < idx (the j = 0 and j = idx terms are the 2u₀ψ in 𝓛),
# padded onto the length-N unit grid. Returns (source, its leading power).
function _painleve1_two_source(Φ_of, idx::NTuple{2,Int}, N::Integer, ::Type{T}) where {T}
    prods = LogSeries{T}[]
    for n1 in 0:idx[1], n2 in 0:idx[2]
        j = (n1, n2)
        jc = (idx[1] - n1, idx[2] - n2)
        (j == (0, 0) || jc == (0, 0)) && continue
        p = Φ_of(j) * Φ_of(jc)
        _is_zero_series(p) || push!(prods, p)
    end
    isempty(prods) && return (nothing, nothing)
    og = minimum(Φ -> minimum(b -> b.power_offset, log_blocks(Φ)), prods)
    pad = LogSeries(FormalSeries(zeros(T, N), :s; power_offset = og))
    (-reduce(+, prods; init = pad), og)
end

function _painleve1_two_sector(sofar::Dict{NTuple{2,Int},LogSeries{T}},
                               idx::NTuple{2,Int}, A, u0::FormalSeries{T},
                               c::Rational{Int}, N::Integer) where {T}
    k = idx[1] - idx[2]
    δ = abs(k) == 1 ? 3 : -2                # differential (|k|=1) vs algebraic diagonal
    zero_ls = LogSeries(FormalSeries(zeros(T, 1), :s))
    src, og = _painleve1_two_source(j -> get(sofar, j, zero_ls), idx, N, T)
    # a sourceless sector is pure kernel, hence starts at the kernel power
    β = src === nothing ? _PI_KERNEL_POWER : og - δ
    g = src === nothing ?
        LogSeries(FormalSeries(zeros(T, N), :s; power_offset = β + δ)) : src
    kernel, homog = nothing, nothing
    if abs(k) == 1
        m = _PI_KERNEL_POWER - β
        if denominator(m) == 1 && 0 ≤ numerator(m) < N
            kernel = Int(numerator(m))
            h = get(sofar, k == 1 ? (1, 0) : (0, 1), nothing)
            h === nothing || (homog = FormalSeries(zeros(T, N), :s; power_offset = β) +
                                      log_block(h, 0))
        end
    end
    freeval = (idx == (1, 0) || idx == (0, 1)) ? one(T) : zero(T)
    _solve_log_sector(ψ -> _painleve1_apply_L(ψ, k, A, u0, c),
                      ψ -> _painleve1_apply_C(ψ, k, A, c),
                      ψ -> _painleve1_apply_D(ψ, c),
                      g, β, δ, N; kernel, homog, freeval)
end

"""
    painleve1_two_parameter(n; orders = 3, convention = :gikm, T = BigFloat)
        -> MultiTransseries

The **two-parameter** Painlevé I transseries
``u = \\sum_{p,q ≥ 0} σ_1^p σ_2^q e^{-(p-q)A/ħ}\\, Φ_{(p,q)}(s)``, solved by
[`multi_transseries_solve`](@ref) in the variable `s = z^{-1/4}` (`ħ = s^5 = z^{-5/4}`),
for every multi-index of total degree up to `orders`. The argument `n` is the number of
fluctuation orders each sector carries.

The action lattice is `(A, −A)` with `A = ` [`painleve1_action`](@ref), which is
**resonant**: the kernel `ℤ(1,1)` (see [`resonance_lattice`](@ref)) makes distinct sectors
share an exponential weight, and sectors acquire powers of `log s`. Sectors are therefore
[`LogSeries`](@ref); the `(0,0)` and `(1,0)`/`(0,1)` sectors reproduce
[`painleve1`](@ref)'s.

The tower starts at `(2,1)` and `(1,2)`, the first multi-indices whose source lands on the
kernel of the one-instanton operator, and its degree is bounded by
[`resonance_depth`](@ref)` = min(p,q)`. **The resonant diagonal `p = q` is log-free at
every order**: the growth comes from the kernel of the `|p−q| = 1` operator, and the
zero-weight operator has none - so the bound is saturated off the diagonal and strict on
it. Degree 2 first appears at `(3,2)`.

Logs are in the **series variable** `s`; convert with `log ħ = 5 log s`,
`log z = −4 log s`.
"""
function painleve1_two_parameter(n::Integer; orders::Integer = 3,
                                 convention::Symbol = :gikm, T::Type = BigFloat)
    n ≥ 1 || throw(InvalidArgument("need n ≥ 1 fluctuation orders, got $n"))
    orders ≥ 0 || throw(InvalidArgument("need orders ≥ 0, got $orders"))
    c = _pi_convention_c(convention)
    A = painleve1_action(; convention, T = BigFloat)
    N = 10 * n + 20 + 10 * orders
    u0 = _painleve1_u0(n; T = T, pad = N - (10 * (n - 1) + 1), convention)
    solve = (sofar, idx) -> _painleve1_two_sector(sofar, idx, A, u0, c, N)
    multi_transseries_solve(LogSeries(u0), (A, -A), solve; orders, var = :s)
end

"""
    _painleve1_two_residual(mt, idx; convention = :gikm) -> LogSeries

The Painlevé I residual of [`painleve1_two_parameter`](@ref) graded at multi-index `idx`:
``\\hat{𝓛}_k Φ_{idx} + \\sum' Φ_j Φ_{idx-j}`` with `k = idx[1] - idx[2]`, where
``\\hat{𝓛}_k`` is the sector operator including its log commutators. It vanishes to
truncation order for every sector - the "every theorem is a test" oracle for the
two-parameter solve.
"""
function _painleve1_two_residual(mt::MultiTransseries{T,A,2}, idx::NTuple{2,Int};
                                 convention::Symbol = :gikm) where {T,A}
    c = _pi_convention_c(convention)
    act = actions(mt)[1]
    u0 = log_block(sector(mt, (0, 0)), 0)
    # the perturbative sector is graded against z, not against a source: use the
    # one-parameter residual, which is the same equation at k = 0 minus u₀²  + s^{-4}
    idx == (0, 0) && return LogSeries(_painleve1_residual(u0; convention))
    k = idx[1] - idx[2]
    ψ = sector(mt, idx)
    P = log_degree(ψ)
    blocks = map(0:P) do p
        r = _painleve1_apply_L(log_block(ψ, p), k, act, u0, c)
        p + 1 ≤ P && (r = r + (p + 1) * _painleve1_apply_C(log_block(ψ, p + 1), k, act, c))
        p + 2 ≤ P &&
            (r = r + ((p + 2) * (p + 1)) * _painleve1_apply_D(log_block(ψ, p + 2), c))
        r
    end
    res = LogSeries(blocks, :s)
    src, _ = _painleve1_two_source(j -> sector(mt, j), idx, n_terms(u0), T)
    src === nothing ? res : res - src
end

# -- the Kapaev Stokes constant -----------------------------------------------------

# Γ(m + 1/2) = (2m)!/(4^m m!) √π - exact rational times √π, so the large-order fit needs
# no SpecialFunctions dependency.
_gamma_half(m::Integer, ::Type{T}) where {T} =
    (factorial(big(2m)) // (big(4)^m * factorial(big(m)))) * sqrt(T(π))

# the perturbative series and the one-instanton fluctuation on a common grid, both as
# ħ-series with the shared s^{-2} prefactor stripped: Φ₀ = Σ aₙ ħ^{2n}, Φ₁ = ħ^{1/2} Σ b_h ħ^h
function _painleve1_stokes_data(nterms::Integer, fluctuations::Integer,
                                convention::Symbol, ::Type{T}) where {T}
    c = _pi_convention_c(convention)
    A = painleve1_action(; convention, T)
    N = 10 * nterms + 5 * fluctuations + 25
    u0 = _painleve1_u0(nterms; T = T, pad = N - (10 * (nterms - 1) + 1), convention)
    zero_src = FormalSeries(zeros(T, N), :s; power_offset = _PI_KERNEL_POWER)
    Φ1 = _solve_linear_sector(ψ -> _painleve1_apply_L(ψ, 1, A, u0, c), zero_src,
                              _PI_KERNEL_POWER, 3, N; free = [0], freeval = one(T))
    a = T[u0[10n] for n in 0:(nterms - 1)]
    b = T[Φ1[5h] for h in 0:fluctuations]
    (A, a, b)
end

"""
    painleve1_stokes_constant(; method = :large_order, nterms = 40, fluctuations = 8,
                              order = 6, ħ = 3//10, convention = :gikm, T = BigFloat)
        -> Complex{T}

The Painlevé I Stokes constant `S₁`, computed numerically by either of two independent
routes. The GIKM value is `S₁ = −i·3^{1/4}/(2√π) ≈ −0.371258 i`
(Garoufalidis–Its–Kapaev–Mariño 2010).

`method = :large_order` (the accurate one, 12+ digits) reads `S₁` off the growth of the
perturbative coefficients, using the exact one-instanton fluctuation as the subleading
tail:

``a_n ≃ 2\\,\\frac{S_1}{2πi} \\sum_{h=0}^{H} b_h\\,
\\frac{Γ(2n - 1/2 - h)}{A^{2n - 1/2 - h}}``

where `a_n` are the coefficients of `Φ₀ = Σ a_n ħ^{2n}` and `b_h` those of
`Φ₁ = ħ^{1/2} Σ b_h ħ^h`. `nterms` sets the perturbative order, `fluctuations` the tail
depth `H`, and `order` the Richardson order accelerating the resulting sequence in `n`.

`method = :lateral` instead *measures* the Stokes jump: it Borel–Padé sums `Φ₀` on both
sides of the positive real axis and divides the [`stokes_discontinuity`](@ref) by the
optimally-truncated one-instanton contribution `e^{-A/ħ} Φ₁(ħ)` at the given `ħ`. This is
the two-lateral route and is genuinely independent of the large-order fit - it reaches
about `1e-6`, limited by the Padé approximant rather than by precision, so `T = Float64`
is the sensible choice (extended-precision Padé is far slower and no more accurate here).

The two routes differ by exactly the factor `2` in the formula above, and that is the
check on it: `Φ₀` has only even powers of `ħ`, so **both** Borel singularities `±A`
feed its large-order growth, while the discontinuity across `arg ζ = 0` sees only `+A`.
"""
function painleve1_stokes_constant(; method::Symbol = :large_order,
                                   nterms::Integer = 40, fluctuations::Integer = 8,
                                   order::Integer = 6, ħ = 3 // 10,
                                   convention::Symbol = :gikm, T::Type = BigFloat)
    fluctuations ≥ 0 || throw(InvalidArgument("need fluctuations ≥ 0, got $fluctuations"))
    if method === :large_order
        window = order + 8
        nterms ≥ window + 2 ||
            throw(InvalidArgument("need nterms ≥ $(window + 2) for Richardson order " *
                                  "$order, got $nterms"))
        A, a, b = _painleve1_stokes_data(nterms, fluctuations, convention, T)
        seq = map((nterms - window):(nterms - 1)) do n
            # Γ(2n − 1/2 − h) = Γ((2n − 1 − h) + 1/2); A^{2n−1/2−h} = A^{2n−h}/√A
            tail = sum(b[h + 1] * _gamma_half(2n - 1 - h, T) * sqrt(A) / A^(2n - h)
                       for h in 0:fluctuations)
            2 * T(π) * im * a[n + 1] / (2 * tail)
        end
        return richardson(seq, order)
    elseif method === :lateral
        nterms ≥ 2 || throw(InvalidArgument("need nterms ≥ 2, got $nterms"))
        h = T(ħ)
        h > 0 || throw(InvalidArgument("need ħ > 0 for the lateral route, got $ħ"))
        A, a, b = _painleve1_stokes_data(nterms, fluctuations, convention, T)
        # Φ₀ on the unit ħ grid: only even powers are populated
        coeffs = zeros(T, 2 * nterms)
        for n in 0:(nterms - 1)
            coeffs[2n + 1] = a[n + 1]
        end
        disc = stokes_discontinuity(borel(FormalSeries(coeffs, :ħ)), h)
        # e^{-A/ħ} Φ₁(ħ), truncated at its smallest term
        terms = [b[j + 1] * h^j for j in 0:fluctuations]
        kopt = argmin(abs.(terms))
        return disc / (sqrt(h) * sum(terms[1:kopt]) * exp(-A / h))
    end
    throw(InvalidArgument("unknown method :$method (expected :large_order or :lateral)"))
end
