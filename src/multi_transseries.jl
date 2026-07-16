# General k-parameter transseries F(σ, ħ) = Σ_{n ∈ ℤ_{≥0}^K} (∏_i σ_i^{n_i}) e^{-(n·A)/ħ}
# Φ_n(ħ) over an action lattice A = (A_1,…,A_K). Sectors live in a Dict keyed by the
# multi-index n (an opaque, structured key - no contiguous-integer assumption), so the
# support is sparse and absent-inside-the-bounding-box keys are exact zeros. This is the
# additive M6 generalization of the locked one-parameter `Transseries`: the K = 1 case
# embeds a `Transseries` and round-trips back to it. The parameter vector σ is not
# stored; it enters at summation time, where the Stokes automorphism acts (multi_alien.jl).

"""
    MultiTransseries{T,A,K,S} <: AbstractSeries

A `K`-parameter transseries ``F(σ, ħ) = \\sum_{n ∈ ℤ_{≥0}^K} \\Big(\\prod_i σ_i^{n_i}\\Big)
e^{-(n·A)/ħ}\\, Φ_n(ħ)`` with action lattice `actions::NTuple{K,A}` and sectors
`sectors::Dict{NTuple{K,Int},S}` (multi-index `n ↦ Φ_n`). Absent keys inside the
componentwise bounding box of the stored support are **exact zeros**; keys outside it are
unknown (truncated). The parameter vector `σ` is not stored - it is supplied to
[`transseries_sum`](@ref), and the Stokes automorphism acts on it (see
[`stokes_automorphism`](@ref)).

The sector type `S` is `FormalSeries{T}` for the ordinary non-resonant case and
[`LogSeries`](@ref)`{T}` when any sector carries `log ħ` blocks - which resonance forces,
see [`resonant_solve`](@ref). It is homogeneous (a log-free sector inside a resonant
transseries is a one-block `LogSeries`), so the `Dict` stays type-stable, and it is
inferred from the sectors you pass; `MultiTransseries{T,A,K}` remains a valid type to
dispatch on.

`K = 1` embeds the one-parameter [`Transseries`](@ref): `MultiTransseries(F)` and
`Transseries(mt)` are inverse (a mandatory round-trip).

# Examples
```julia
mt = MultiTransseries((1//1, 2//1), Dict((0,0) => Φ00, (1,0) => Φ10))
sector(mt, (1, 0))                                # Φ_{(1,0)}
weight(mt, (1, 1))                                # n·A = 1·1 + 1·2 = 3
MultiTransseries(Transseries(:euler, 12))         # rank-1 embedding
```
"""
struct MultiTransseries{T,A<:Number,K,S<:AbstractSeries} <: AbstractSeries
    actions::NTuple{K,A}
    sectors::Dict{NTuple{K,Int},S}
    var::Symbol

    function MultiTransseries{T,A,K,S}(actions::NTuple{K,A},
                                       sectors::Dict{NTuple{K,Int},S},
                                       var::Symbol) where {T,A<:Number,K,S<:AbstractSeries}
        K ≥ 1 || throw(InvalidArgument("a MultiTransseries needs rank K ≥ 1, got $K"))
        any(iszero, actions) &&
            throw(InvalidArgument("a MultiTransseries needs nonzero actions, got A = $actions"))
        isempty(sectors) &&
            throw(InvalidArgument("a MultiTransseries needs ≥ 1 sector"))
        canon = Dict{NTuple{K,Int},S}()
        for (idx, Φ) in sectors
            any(<(0), idx) &&
                throw(InvalidArgument("sector index components must be ≥ 0, got $idx"))
            variable(Φ) == var || throw(IncompatibleSeries(:var, var, variable(Φ)))
            _is_zero_series(Φ) && continue              # drop structural zeros
            canon[idx] = Φ
        end
        if isempty(canon)                               # keep a representable origin zero
            canon[ntuple(_ -> 0, K)] = _zero_sector_like(first(values(sectors)), var)
        end
        new{T,A,K,S}(actions, canon, var)
    end
end

# a zero sector of the same type as an existing one (S is homogeneous across the Dict)
_zero_sector_like(Φ::FormalSeries, var::Symbol) = FormalSeries([zero(Φ.coeffs[1])], var)
_zero_sector_like(L::LogSeries, var::Symbol) =
    LogSeries([_zero_sector_like(L.blocks[1], var)], var)

# Type parameters K, A, T, S are inferred from the values (an `NTuple{K,A}` signature would
# leave A unbound for the empty-tuple K = 0 case; here rank/action/coefficient/sector types
# are all recovered at runtime and passed explicitly to the inner constructor).
function MultiTransseries(actions::Tuple, sectors::AbstractDict;
                          var::Union{Nothing,Symbol} = nothing)
    isempty(actions) &&
        throw(InvalidArgument("a MultiTransseries needs rank K ≥ 1, got 0 actions"))
    isempty(sectors) && throw(InvalidArgument("a MultiTransseries needs ≥ 1 sector"))
    K = length(actions)
    A = mapreduce(typeof, promote_type, actions)
    A <: Number || throw(InvalidArgument("actions must be numbers, got element type $A"))
    acts = map(a -> convert(A, a), actions)
    T = mapreduce(Φ -> eltype(typeof(Φ)), promote_type, values(sectors))
    v = var === nothing ? variable(first(values(sectors))) : var
    # one log block anywhere makes every sector a LogSeries - S is homogeneous
    S = any(Φ -> Φ isa LogSeries, values(sectors)) ? LogSeries{T} : FormalSeries{T}
    typed = Dict{NTuple{K,Int},S}(NTuple{K,Int}(k) => _promote_sector(Φ, S)
                                  for (k, Φ) in sectors)
    MultiTransseries{T,A,K,S}(acts, typed, v)
end

# widen a sector to the transseries' homogeneous sector type
_promote_sector(Φ::FormalSeries, ::Type{FormalSeries{Tp}}) where {Tp} =
    _promote_series(Φ, Tp)
_promote_sector(Φ::FormalSeries, ::Type{LogSeries{Tp}}) where {Tp} =
    LogSeries(_promote_series(Φ, Tp))
_promote_sector(L::LogSeries, ::Type{LogSeries{Tp}}) where {Tp} = _promote_series(L, Tp)

# -- embedding of the one-parameter Transseries ------------------------------------

"""
    MultiTransseries(F::Transseries) -> MultiTransseries{T,A,1}

Embed a one-parameter [`Transseries`](@ref) as the rank-1 case (keys `(n,)`). Inverse
of [`Transseries(::MultiTransseries)`](@ref); the round-trip is exact.
"""
function MultiTransseries(F::Transseries{T,A}) where {T,A}
    d = Dict{NTuple{1,Int},FormalSeries{T}}()
    for n in 0:(n_sectors(F) - 1)
        d[(n,)] = sector(F, n)
    end
    MultiTransseries((action(F),), d; var = variable(F))
end

"""
    Transseries(mt::MultiTransseries{T,A,1}) -> Transseries

Project a rank-1 `MultiTransseries` back to a one-parameter [`Transseries`](@ref)
(inverse of the embedding). Throws `InvalidArgument` for `K ≠ 1`, and
`IncompatibleSeries` when a sector carries genuine `log ħ` blocks (the one-parameter
`Transseries` has no representation for those).
"""
function Transseries(mt::MultiTransseries{T,A,K}) where {T,A,K}
    K == 1 || throw(InvalidArgument("only a rank-1 MultiTransseries projects to a " *
                                    "one-parameter Transseries, got rank $K"))
    idxs = sector_indices(mt)
    N = idxs[end][1] + 1
    secs = [_as_formal(get(mt.sectors, (n,), _zero_sector(mt))) for n in 0:(N - 1)]
    Transseries(mt.actions[1], secs, mt.var)
end

_as_formal(Φ::FormalSeries) = Φ
_as_formal(L::LogSeries) = FormalSeries(L)

# -- accessors ---------------------------------------------------------------------

"""
    actions(mt::MultiTransseries) -> NTuple{K}

The action lattice ``A = (A_1,…,A_K)`` of the exponential weights ``e^{-(n·A)/ħ}``.
"""
actions(mt::MultiTransseries) = mt.actions

"""
    n_actions(mt::MultiTransseries) -> Int

The lattice rank ``K`` (number of independent instanton actions / parameters).
"""
n_actions(::MultiTransseries{T,A,K}) where {T,A,K} = K

"""
    n_sectors(mt::MultiTransseries) -> Int

Number of stored (nonzero) sectors.
"""
n_sectors(mt::MultiTransseries) = length(mt.sectors)

"""
    sector(mt::MultiTransseries, n) -> FormalSeries or LogSeries

The sector ``Φ_n`` at multi-index `n` (an `NTuple` or varargs of nonnegative integers);
an absent key returns the zero series, mirroring `Φ[n]` on [`FormalSeries`](@ref). The
returned type is the transseries' sector type `S` (see [`MultiTransseries`](@ref)).
"""
function sector(mt::MultiTransseries{T,A,K}, n::NTuple{N,Integer}) where {T,A,K,N}
    N == K || throw(InvalidArgument("sector index has rank $N, expected $K"))
    any(<(0), n) &&
        throw(InvalidArgument("sector index components must be ≥ 0, got $n"))
    get(mt.sectors, NTuple{K,Int}(n), _zero_sector(mt))
end
sector(mt::MultiTransseries, n::Integer...) = sector(mt, n)

_zero_sector(mt::MultiTransseries) =
    _zero_sector_like(first(values(mt.sectors)), mt.var)

"""
    sectors(mt::MultiTransseries) -> Dict

The sector map ``n ↦ Φ_n`` (a copy; transseries are immutable).
"""
sectors(mt::MultiTransseries) = copy(mt.sectors)

"""
    sector_indices(mt::MultiTransseries) -> Vector{NTuple{K,Int}}

The stored multi-indices, sorted lexicographically (deterministic iteration order).
"""
sector_indices(mt::MultiTransseries) = sort!(collect(keys(mt.sectors)))

"""
    weight(mt::MultiTransseries, n) -> Number

The transmonomial exponent ``n·A = \\sum_i n_i A_i`` of the sector at multi-index `n`.
"""
weight(mt::MultiTransseries{T,A,K}, n::NTuple{K,Integer}) where {T,A,K} =
    sum(n .* mt.actions)

variable(mt::MultiTransseries) = mt.var
is_exact(mt::MultiTransseries) = all(is_exact, values(mt.sectors))
Base.eltype(::Type{MultiTransseries{T,A,K,S}}) where {T,A,K,S} = T

"""
    log_degree(mt::MultiTransseries[, n]) -> Int

The `log ħ` degree of the sector at multi-index `n`, or the maximum over all stored
sectors when `n` is omitted. Always `0` for a non-resonant transseries (sector type
`FormalSeries`). Bounded by [`resonance_depth`](@ref).
"""
log_degree(mt::MultiTransseries) = maximum(_log_degree, values(mt.sectors))
log_degree(mt::MultiTransseries{T,A,K}, n::NTuple{K,Integer}) where {T,A,K} =
    _log_degree(sector(mt, n))
log_degree(mt::MultiTransseries, n::Integer...) = log_degree(mt, n)

_log_degree(::FormalSeries) = 0
_log_degree(L::LogSeries) = log_degree(L)

Base.:(==)(F::MultiTransseries, G::MultiTransseries) =
    F.var == G.var && F.actions == G.actions && F.sectors == G.sectors

"""
    is_resonant(mt::MultiTransseries) -> Bool

`true` when two distinct stored sectors share the same transmonomial weight ``n·A`` -
the **support-level** resonance condition (distinct instanton numbers meeting on one ray).

This asks only about the sectors actually stored. The stronger, structural question - can
*any* pair of sectors collide, i.e. does the weight map have a nonzero kernel - is
[`is_resonant(::Tuple)`](@ref) on the actions, whose kernel is
[`resonance_lattice`](@ref).
"""
function is_resonant(mt::MultiTransseries)
    seen = Set{Any}()
    for idx in keys(mt.sectors)
        w = weight(mt, idx)
        w in seen && return true
        push!(seen, w)
    end
    false
end

# -- arithmetic --------------------------------------------------------------------

function _check_compatible(F::MultiTransseries{T,A,K},
                           G::MultiTransseries{S,B,L}) where {T,A,K,S,B,L}
    K == L || throw(IncompatibleSeries(:rank, K, L))
    F.var == G.var || throw(IncompatibleSeries(:var, F.var, G.var))
    F.actions == G.actions || throw(IncompatibleSeries(:actions, F.actions, G.actions))
    nothing
end

_promote_series(Φ::FormalSeries{T}, ::Type{Tp}) where {T,Tp} =
    T === Tp ? Φ :
    FormalSeries(convert(Vector{Tp}, Φ.coeffs), Φ.var; power_offset = Φ.power_offset)

# the common sector type of two operands: LogSeries wins (a log-free sector embeds in it)
_common_sector_type(::Type{<:FormalSeries}, ::Type{<:FormalSeries}, ::Type{Tp}) where {Tp} =
    FormalSeries{Tp}
_common_sector_type(::Type{<:AbstractSeries}, ::Type{<:AbstractSeries}, ::Type{Tp}) where {Tp} =
    LogSeries{Tp}

# addition promotes the coefficient type (e.g. a rational F plus a complex alien
# derivative in the operator-form Stokes automorphism) and the sector type. Per-arg ranks
# K1, K2 so a mismatch reaches `_check_compatible` (a typed `:rank` error) rather than a
# MethodError; past that check K1 == K2.
function Base.:+(F::MultiTransseries{T1,A1,K1,S1},
                 G::MultiTransseries{T2,A2,K2,S2}) where {T1,A1,K1,S1,T2,A2,K2,S2}
    _check_compatible(F, G)
    Tp = promote_type(T1, T2)
    S = _common_sector_type(S1, S2, Tp)
    d = Dict{NTuple{K1,Int},S}()
    for idx in union(keys(F.sectors), keys(G.sectors))
        inF, inG = haskey(F.sectors, idx), haskey(G.sectors, idx)
        d[idx] = inF && inG ?
                 _promote_sector(F.sectors[idx], S) + _promote_sector(G.sectors[idx], S) :
                 inF ? _promote_sector(F.sectors[idx], S) : _promote_sector(G.sectors[idx], S)
    end
    MultiTransseries(F.actions, d; var = F.var)
end

Base.:-(F::MultiTransseries{T,A,K,S}) where {T,A,K,S} =
    MultiTransseries(F.actions, Dict{NTuple{K,Int},S}(idx => -Φ
                     for (idx, Φ) in F.sectors); var = F.var)
Base.:-(F::MultiTransseries, G::MultiTransseries) = F + (-G)

# value type inferred: λ may promote the coefficient type (e.g. rational → complex)
Base.:*(λ::Number, F::MultiTransseries) =
    MultiTransseries(F.actions, Dict(idx => λ * Φ for (idx, Φ) in F.sectors); var = F.var)
Base.:*(F::MultiTransseries, λ::Number) = λ * F

# all multi-indices 0 ≼ i ≼ cap (componentwise), as a flat vector
_box_indices(cap::NTuple{K,Int}) where {K} =
    vec(collect(Iterators.product((0:c for c in cap)...)))::Vector{NTuple{K,Int}}

_bounding_box(mt::MultiTransseries{T,A,K}) where {T,A,K} =
    ntuple(i -> maximum(idx[i] for idx in keys(mt.sectors)), K)

"""
    F * G

Weight-graded Cauchy product ``(F·G)_n = \\sum_{i+j=n} Φ_i Ψ_j`` (the exponential
weights add: ``e^{-(i·A)/ħ} e^{-(j·A)/ħ} = e^{-(n·A)/ħ}``), truncated to the
componentwise `min` of the two operands' bounding boxes - the frontier where every
contribution is available. Reduces to the one-parameter `min(n_sectors,…)` at `K = 1`.
"""
function Base.:*(F::MultiTransseries{T1,A1,K1,S1},
                 G::MultiTransseries{T2,A2,K2,S2}) where {T1,A1,K1,S1,T2,A2,K2,S2}
    _check_compatible(F, G)
    Tp = promote_type(T1, T2)
    S = _common_sector_type(S1, S2, Tp)
    cap = min.(_bounding_box(F), _bounding_box(G))
    d = Dict{NTuple{K1,Int},S}()
    for n in _box_indices(cap)
        terms = S[]
        for i in _box_indices(n)
            j = n .- i
            fi = get(F.sectors, i, nothing)
            gj = get(G.sectors, j, nothing)
            (fi === nothing || gj === nothing) && continue
            push!(terms, _promote_sector(fi, S) * _promote_sector(gj, S))
        end
        isempty(terms) || (d[n] = reduce(+, terms))
    end
    if isempty(d)
        d[ntuple(_ -> 0, K1)] = _promote_sector(FormalSeries([zero(Tp)], F.var), S)
    end
    MultiTransseries(F.actions, d; var = F.var)
end

"""
    truncate(mt::MultiTransseries; degree = nothing, terms = nothing) -> MultiTransseries

The transseries truncated to the sectors with multi-index `≼ degree` componentwise
(`degree` an `Integer`, applied to every direction, or an `NTuple`; `nothing` keeps all),
each sector truncated to `terms` coefficients (`nothing` keeps every sector's length).
"""
function Base.truncate(mt::MultiTransseries{T,A,K,S};
                       degree::Union{Nothing,Integer,NTuple{K,Integer}} = nothing,
                       terms::Union{Nothing,Integer} = nothing) where {T,A,K,S}
    cap = degree === nothing ? nothing :
          degree isa Integer ? ntuple(_ -> Int(degree), K) : NTuple{K,Int}(degree)
    d = Dict{NTuple{K,Int},S}()
    for (idx, Φ) in mt.sectors
        cap === nothing || all(idx .≤ cap) || continue
        d[idx] = terms === nothing ? Φ : _truncate_sector(Φ, terms)
    end
    if isempty(d)
        d[ntuple(_ -> 0, K)] = _zero_sector(mt)
    end
    MultiTransseries(mt.actions, d; var = mt.var)
end

_truncate_sector(Φ::FormalSeries, terms::Integer) = truncate(Φ, min(terms, n_terms(Φ)))
_truncate_sector(L::LogSeries, terms::Integer) = truncate(L; terms)

# -- summation ---------------------------------------------------------------------

"""
    transseries_sum(F::MultiTransseries, σ, ħ; θ = 0, side = :plus, log_ħ = log(complex(ħ)),
                    kwargs...) -> Complex

Lateral Borel sum of the `K`-parameter transseries at parameter vector `σ` (an `NTuple`
of length `K`; a scalar is accepted for `K = 1`):
``\\sum_n \\big(\\prod_i σ_i^{n_i}\\big) e^{-(n·A)/ħ}\\, s_θ^{±}(Φ_n)(ħ)``. Each
non-constant sector is summed by [`lateral_sum`](@ref) of its Borel transform along the
ray `θ`; constant sectors are evaluated exactly. Several Borel singularities per ray are
handled sector by sector. See [`stokes_automorphism`](@ref) for the side↔parameter-shift
relation.

Resonant sectors carry `log ħ` blocks ([`LogSeries`](@ref)). `log ħ` is inert under
Borel–Laplace - it does not depend on the Borel variable ζ - so each block is summed
independently and weighted by `log_ħ^p`. `log_ħ` defaults to the **principal** branch,
which is correct for `θ ∈ (-π, π]`; a transseries continued past that (across several
Stokes rays) lives on the universal cover of `ħ = 0`, and you pass the continued value
explicitly. The monodromy `log_ħ ↦ log_ħ + 2πi` genuinely changes the sum - that *is*
resonance.
"""
function transseries_sum(F::MultiTransseries{T,A,K}, σ::Tuple, ħ::Number;
                         θ::Real = 0, side::Symbol = :plus,
                         log_ħ = log(complex(ħ)), kwargs...) where {T,A,K}
    length(σ) == K ||
        throw(InvalidArgument("parameter vector σ has length $(length(σ)), expected $K"))
    RT = _real_float_type(ħ)
    scale = one(Complex{RT})
    for s in σ; scale *= one(s); end
    for a in F.actions; scale *= one(a); end
    total = zero(scale)
    for (idx, Φ) in F.sectors
        _is_zero_series(Φ) && continue
        w = prod(σ[i]^idx[i] for i in 1:K) * exp(-weight(F, idx) / ħ)
        iszero(w) && continue
        total += w * _sum_sector(Φ, ħ, log_ħ; θ, side, kwargs...)
    end
    total
end

_sum_sector(Φ::FormalSeries, ħ, log_ħ; θ, side, kwargs...) =
    n_terms(Φ) == 1 ? evaluate(Φ, ħ) : lateral_sum(borel(Φ), ħ; θ, side, kwargs...)

function _sum_sector(L::LogSeries, ħ, log_ħ; θ, side, kwargs...)
    total = _sum_sector(L.blocks[1], ħ, log_ħ; θ, side, kwargs...) * one(log_ħ)
    for p in 1:log_degree(L)
        Φ = L.blocks[p + 1]
        _is_zero_series(Φ) && continue
        total += log_ħ^p * _sum_sector(Φ, ħ, log_ħ; θ, side, kwargs...)
    end
    total
end

transseries_sum(F::MultiTransseries{T,A,1}, σ::Number, ħ::Number; kwargs...) where {T,A} =
    transseries_sum(F, (σ,), ħ; kwargs...)
