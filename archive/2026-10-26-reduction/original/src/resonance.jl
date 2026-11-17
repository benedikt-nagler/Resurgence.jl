# Resonance: when the weight map n ↦ n·A stops being injective, and the log sectors it
# forces. The mechanism, which drives every design choice here:
#
# Conjugating a sector's own transmonomial away, every sector equation is linear in the
# channel operator (ħ²∂_ħ + c) with detuning c = n·A − A_j, since e^{w/ħ}(ħ²∂_ħ)e^{-w/ħ}
# = ħ²∂_ħ + w. Then:
#
#   c ≠ 0 - the recursion c·u_m + (β+m−1)·u_{m−1} = g_m is triangular with nonzero
#           diagonal, so there is a unique power-series solution and NO logs. This is the
#           theorem that lets the non-resonant M2/M6a layers ignore logs entirely.
#   c = 0 - ħ²∂_ħ maps ħ^k ↦ k·ħ^{k+1}, so its image misses ħ¹ exactly. The cokernel is
#           ℂ·ħ and its preimage is log ħ (ħ²∂_ħ log ħ = ħ). A source with a nonzero ħ¹
#           coefficient is the secular term, and it forces a log.
#
# Products and alien derivatives never generate logs; only inversion does.
#
# Two consequences worth stating because they become tests:
#  * zero detuning is necessary but NOT sufficient - if power_offset(g) ∉ ℤ the power grid
#    never hits ħ¹, so ħ²∂_ħ stays invertible at c = 0 and a fractional-power sector is
#    log-free.
#  * the detuning is constant along the kernel Λ = {k ∈ ℤ^K : k·A = 0}, so log degree grows
#    by one per resonant inversion - hence log_degree(Φ_n) ≤ resonance_depth(n), which for
#    A = (A,−A) with k = (1,1) is the familiar min(n₁,n₂) bound.

# -- the resonance lattice ---------------------------------------------------------

# Rational rows of the ℚ-linear map k ↦ k·A: one row for real actions, two (re, im) for
# complex ones. `nothing` when an action is not exactly rational.
function _action_rows(actions::NTuple{K,Any}) where {K}
    rows = Vector{Vector{Rational{BigInt}}}()
    re = Vector{Rational{BigInt}}(undef, K)
    im = Vector{Rational{BigInt}}(undef, K)
    anyim = false
    for (i, a) in pairs(actions)
        r, c = real(a), imag(a)
        (_is_rational(r) && _is_rational(c)) || return nothing
        re[i] = Rational{BigInt}(r)
        im[i] = Rational{BigInt}(c)
        iszero(c) || (anyim = true)
    end
    push!(rows, re)
    anyim && push!(rows, im)
    rows
end

_is_rational(x::Union{Integer,Rational}) = true
_is_rational(::Any) = false

# Reduced row echelon form over ℚ; returns (rref, pivot columns).
function _rref!(M::Matrix{Rational{BigInt}})
    nr, nc = size(M)
    pivots = Int[]
    r = 1
    for c in 1:nc
        r > nr && break
        p = findfirst(i -> !iszero(M[i, c]), r:nr)
        p === nothing && continue
        p += r - 1
        M[r, :], M[p, :] = M[p, :], M[r, :]
        M[r, :] ./= M[r, c]
        for i in 1:nr
            i == r && continue
            iszero(M[i, c]) || (M[i, :] .-= M[i, c] .* M[r, :])
        end
        push!(pivots, c)
        r += 1
    end
    M, pivots
end

"""
    resonance_lattice(actions; lattice = nothing) -> Vector{NTuple{K,Int}}
    resonance_lattice(mt::MultiTransseries) -> Vector{NTuple{K,Int}}

A ℤ-basis of the **resonance lattice** ``Λ = \\{k ∈ ℤ^K : k·A = 0\\}`` - the kernel of the
weight map. `Λ ≠ 0` is exactly the condition for distinct sectors to share a transmonomial,
i.e. for [`is_resonant`](@ref) to hold and for logs to appear; `Λ = 0` (always the case at
`K = 1`, since `A_1 ≠ 0`) means the transseries is log-free.

Each basis vector is sign-canonicalized (first nonzero entry positive). Note that a kernel
vector need not be `≽ 0`: `A = (1, 2)` has `Λ = ℤ·(2,-1)`, so sectors `(2,0)` and `(0,1)`
collide even though no *positive* combination vanishes. [`resonance_depth`](@ref) needs a
`≽ 0` generator and says so when there is none.

**Exact actions only.** Rational (or complex-rational) actions give the kernel exactly;
inexact ones - Painlevé I's `A = 8√6/5` among them - throw, since deciding a ℚ-linear
relation between floats is not a well-posed question. Pass the kernel yourself via
`lattice` in that case (every consumer here takes the same keyword).
"""
function resonance_lattice(actions::Tuple; lattice = nothing)
    lattice === nothing || return _check_lattice(lattice, actions)
    K = length(actions)
    rows = _action_rows(actions)
    rows === nothing &&
        throw(InvalidArgument("resonance_lattice needs exactly-rational actions (got " *
                              "$(actions)); a ℚ-linear relation between inexact actions " *
                              "is not decidable - pass the kernel explicitly, e.g. " *
                              "`lattice = [(1, 1)]`"))
    M = Matrix{Rational{BigInt}}(undef, length(rows), K)
    for (i, row) in pairs(rows); M[i, :] = row; end
    R, pivots = _rref!(M)
    free = setdiff(1:K, pivots)
    isempty(free) && return NTuple{K,Int}[]
    basis = NTuple{K,Int}[]
    for f in free
        # the nullspace vector with 1 in free coordinate f, 0 in the other free ones
        v = zeros(Rational{BigInt}, K)
        v[f] = 1
        for (r, p) in pairs(pivots)
            v[p] = -R[r, f]
        end
        push!(basis, _integralize(v, length(free), actions))
    end
    basis
end

resonance_lattice(mt::MultiTransseries; lattice = nothing) =
    resonance_lattice(actions(mt); lattice)

# Clear denominators and reduce. With 1 in every free coordinate, an *integral* v-basis is
# provably the saturated kernel: any w ∈ ℤ-kernel is w = Σ_f w[f]·v_f with integer
# coefficients. A 1-dimensional kernel is also safe (clearing denominators then dividing
# by the gcd lands on the primitive generator). Higher-dimensional kernels with fractional
# pivots could yield a finite-index sublattice, so we refuse rather than lie.
function _integralize(v::Vector{Rational{BigInt}}, dim::Int, actions)
    if !all(x -> denominator(x) == 1, v) && dim > 1
        throw(InvalidArgument("the resonance lattice of actions $(actions) has dimension " *
                              "$dim with non-integral echelon pivots; saturating it is not " *
                              "implemented - pass `lattice = [...]` explicitly"))
    end
    d = lcm(denominator.(v))
    w = numerator.(v .* d)
    g = foldl(gcd, w; init = big(0))
    iszero(g) || (w = w .÷ g)
    i = findfirst(!iszero, w)
    i === nothing || w[i] > 0 || (w = -w)               # canonical sign
    NTuple{length(v),Int}(Int.(w))
end

function _check_lattice(lattice, actions)
    K = length(actions)
    basis = NTuple{K,Int}[]
    for k in lattice
        length(k) == K ||
            throw(InvalidArgument("lattice vector $k has rank $(length(k)), expected $K"))
        ki = NTuple{K,Int}(k)
        all(iszero, ki) && throw(InvalidArgument("lattice vectors must be nonzero"))
        push!(basis, ki)
    end
    basis
end

"""
    is_resonant(actions::Tuple; lattice = nothing) -> Bool

`true` when the action lattice itself is resonant - the weight map `k ↦ k·A` has a
nonzero integer kernel (see [`resonance_lattice`](@ref)), so *some* pair of distinct
sectors must share a transmonomial and logs are possible.

This is the **lattice-level** predicate. The [`is_resonant(::MultiTransseries)`](@ref)
method is the weaker support-level one: it asks whether two sectors that are actually
*stored* collide.
"""
is_resonant(actions::Tuple; lattice = nothing) =
    !isempty(resonance_lattice(actions; lattice))

"""
    resonance_depth(mt::MultiTransseries, n; lattice = nothing) -> Int

How many times the sector `n` can be walked back along the resonance lattice and stay in
`ℤ_{≥0}^K`: ``\\max\\{m : n - m·k ≽ 0\\}`` for the `≽ 0` generator `k` of
[`resonance_lattice`](@ref). Since each resonant inversion adds one power of `log`, this
bounds the log degree of `Φ_n`, and the bound is attained (see [`resonant_solve`](@ref)).

For `A = (A, -A)` with `k = (1,1)` this is `min(n₁, n₂)`. Returns `0` for a non-resonant
lattice. Throws when the lattice has no `≽ 0` generator - e.g. `A = (1, 2)` (kernel
`(2,-1)`) or duplicated actions `A = (a, a)` (kernel `(1,-1)`): sectors collide there, but
no positive walk generates them, and the correct depth is model-dependent.
"""
function resonance_depth(mt::MultiTransseries{T,A,K}, n::NTuple{K,Integer};
                         lattice = nothing) where {T,A,K}
    basis = resonance_lattice(mt; lattice)
    isempty(basis) && return 0
    pos = filter(k -> all(≥(0), k), basis)
    isempty(pos) &&
        throw(InvalidArgument("the resonance lattice $(basis) of actions $(actions(mt)) " *
                              "has no ≽ 0 generator, so no positive walk reaches a " *
                              "resonant sector; the log-degree bound is model-dependent " *
                              "here - supply it yourself"))
    maximum(pos) do k
        idx = findall(!iszero, k)
        minimum(i -> n[i] ÷ k[i], idx)
    end
end

resonance_depth(mt::MultiTransseries, n::Integer...; lattice = nothing) =
    resonance_depth(mt, n; lattice)

"""
    charges(mt::MultiTransseries, ω) -> Vector{NTuple{K,Int}}

The **fiber** over a singularity location `ω`: every nonzero charge `ℓ ≽ 0` inside the
bounding box of `mt`'s support with `ℓ·A = ω`, sorted lexicographically. At a non-resonant
lattice this is a single charge (or empty); at a resonant one it is the whole tower - for
`A = (A, -A)` the fiber over `ω = A` is `(1,0), (2,1), (3,2), …`, which is precisely why a
resonant `Δ_ω` carries independent Stokes data per charge (see
[`alien_derivative`](@ref)).
"""
function charges(mt::MultiTransseries{T,A,K}, ω::Number) where {T,A,K}
    box = _bounding_box(mt)
    out = NTuple{K,Int}[]
    for ℓ in _box_indices(box)
        all(iszero, ℓ) && continue
        weight(mt, ℓ) == ω && push!(out, ℓ)
    end
    sort!(out)
end

# -- the channel solver ------------------------------------------------------------

# Solve (ħ²∂ + c)u = g for a plain power series.
#   c ≠ 0: u has g's power offset β; c·u_m + (β+m−1)·u_{m−1} = g_m is triangular.
#   c = 0: ħ²∂(ħ^{β-1+m}) = (β−1+m)·ħ^{β+m}, so u has offset β−1 and u_m = g_m/(β−1+m).
#          The one m with β+m = 1 is the obstruction: the caller must have cleared g's ħ¹
#          coefficient (that is the log's job), and u_m there is the free constant, set 0.
function _solve_channel(g::FormalSeries{T}, c::Number) where {T}
    β = g.power_offset
    if !iszero(c)
        # note: `similar` would leave undef refs for non-isbits coefficients
        # (Rational{BigInt} among them), so seed `prev` from the type, not from `u`
        Tu = typeof(zero(T) / c)
        u = Vector{Tu}(undef, length(g.coeffs))
        prev = zero(Tu)
        for m in eachindex(g.coeffs)
            u[m] = (g.coeffs[m] - _power_factor(β + (m - 2)) * prev) / c
            prev = u[m]
        end
        return FormalSeries(u, g.var; power_offset = β)
    end
    coeffs = map(eachindex(g.coeffs)) do m
        f = β - 1 + (m - 1)
        iszero(f) ? zero(g.coeffs[m]) / 1 : g.coeffs[m] / _power_factor(f)
    end
    FormalSeries(collect(coeffs), g.var; power_offset = β - 1)
end

"""
    resonant_solve(g, detuning) -> LogSeries

Solve the transseries channel equation ``(ħ^2 ∂_ħ + c)\\,u = g`` for `u`, where `g` is a
[`FormalSeries`](@ref) or [`LogSeries`](@ref) and `c = detuning`. This is *the* primitive
that creates log sectors - see the discussion at the top of `resonance.jl`.

- `detuning ≠ 0`: the solution is unique and the log degree is **preserved**
  (`log_degree(u) == log_degree(g)`); with `g` log-free, so is `u`.
- `detuning == 0` (resonance): the log degree **grows by one**, with the top block known in
  closed form as ``u^{[P+1]} = ([ħ^1] g^{[P]})/(P+1)``. The remaining integration constant
  (the `ħ^0` coefficient of `u^{[0]}`, which `ħ^2 ∂_ħ` annihilates) is fixed to zero.

At `detuning == 0` the growth is conditional on the source actually hitting the cokernel:
if `power_offset(g) ∉ ℤ` the power grid misses `ħ^1` entirely and `u` stays log-free.

# Examples
```julia
resonant_solve(FormalSeries([0//1, 1//1]), 0//1)   # (ħ²∂)u = ħ  ⇒  u = log ħ
resonant_solve(FormalSeries([1//1]), 2//1)         # detuned  ⇒  log-free
```
"""
resonant_solve(g::FormalSeries, detuning::Number) = resonant_solve(LogSeries(g), detuning)

function resonant_solve(g::LogSeries{T}, detuning::Number) where {T}
    P = log_degree(g)
    if !iszero(detuning)
        # Top-down: block p of (ħ²∂+c)u = g reads (ħ²∂+c)u^{[p]} = g^{[p]} − (p+1)·ħ·u^{[p+1]},
        # and the top block u^{[P+1]} solves (ħ²∂+c)u = 0, hence vanishes. No growth.
        blocks = Vector{Any}(undef, P + 1)
        for p in P:-1:0
            rhs = log_block(g, p)
            p < P && (rhs = rhs - (p + 1) * _shift_power(blocks[p + 2], 1))
            blocks[p + 1] = _solve_channel(rhs, detuning)
        end
        return LogSeries(_unify(blocks), g.var)
    end
    # Resonant. Descending, block p's obstruction ([ħ¹] of its RHS must vanish, since
    # ħ²∂ misses ħ¹) fixes the free constant κ of the block *above* it:
    #   [ħ¹](g^{[p]} − (p+1)·ħ·u^{[p+1]}) = [ħ¹]g^{[p]} − (p+1)·κ = 0.
    # u^{[0]}'s own constant is never fixed by anything below - it is the genuine
    # integration constant, and _solve_channel leaves it zero.
    blocks = Vector{Any}(undef, P + 2)
    blocks[P + 2] = _zero_like(log_block(g, P))
    for p in P:-1:0
        gp = log_block(g, p)
        κ = _coeff_at(gp, 1 // 1) / (p + 1)
        blocks[p + 2] = blocks[p + 2] + FormalSeries([κ], g.var)
        # A zero block above contributes nothing; skip it rather than let the ħ· shift
        # collide with a fractional power_offset (the log-free fractional case, where κ
        # is zero precisely because the power grid misses ħ¹).
        rhs = _is_zero_series(blocks[p + 2]) ? gp :
              gp - (p + 1) * _shift_power(blocks[p + 2], 1)
        blocks[p + 1] = _solve_channel(rhs, detuning)
    end
    LogSeries(_unify(blocks), g.var)
end

# a common coefficient type across the solved blocks (κ and the divisions can widen it)
function _unify(blocks::Vector{Any})
    Tp = mapreduce(Φ -> eltype(typeof(Φ)), promote_type, blocks)
    FormalSeries{Tp}[_promote_series(Φ, Tp) for Φ in blocks]
end
