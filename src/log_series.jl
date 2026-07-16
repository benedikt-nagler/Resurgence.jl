# A sector block polynomial in log(x): Φ(x) = Σ_p log(x)^p Φ^{[p]}(x). This is the
# transmonomial algebra's nilpotent generator - needed exactly when the transseries is
# resonant (n ↦ n·A not injective), see resonance.jl for why. log(x) is inert under
# Borel–Laplace (it does not depend on the Borel variable ζ), so every operation here is
# blockwise in p and the summation seam sums each block separately.
#
# The two indices of a resonant transseries have *different* truncation semantics, which
# is why they live in different containers: the sector multi-index n is truncated (higher
# sectors are unknown), while the log degree p is complete and exact - a block polynomial
# is fully known once known, and log degrees *add* under products.

"""
    LogSeries{T} <: AbstractSeries

A series polynomial in `log(x)` with [`FormalSeries`](@ref) coefficients:
``Φ(x) = \\sum_{p=0}^{P} \\log(x)^p\\, Φ^{[p]}(x)``, stored as `blocks::Vector{FormalSeries{T}}`
with `blocks[p+1]` ``= Φ^{[p]}``. Trailing zero blocks are dropped, so
[`log_degree`](@ref) ``= P`` is the true degree; a log-free `LogSeries` has one block and
[`FormalSeries`](@ref)`(L)` projects it back (the round-trip with `LogSeries(::FormalSeries)`
is exact).

Blocks each carry their own `power_offset` - [`resonant_solve`](@ref) generally produces
blocks on different power grids.

**Convention:** the basis is `log(x)` itself, not `log(1/x)` or `log(x/A)`. A rescaling
`log(c·x^r)` is a unipotent triangular change of basis in the log degree, hence pure
convention; this one is fixed so that ``x^2 ∂_x \\log x = x``.

# Examples
```julia
L = LogSeries([FormalSeries([0//1]), FormalSeries([1//1])])   # log(ħ)
log_degree(L)                                                  # 1
log_block(L, 1)                                                # the coefficient of log(ħ)
evaluate(L, 0.5)                                               # log(0.5) + 0.0
```
"""
struct LogSeries{T} <: AbstractSeries
    blocks::Vector{FormalSeries{T}}
    var::Symbol

    function LogSeries{T}(blocks::Vector{FormalSeries{T}}, var::Symbol) where {T}
        isempty(blocks) && throw(InvalidArgument("a LogSeries needs ≥ 1 block"))
        for Φ in blocks
            variable(Φ) == var || throw(IncompatibleSeries(:var, var, variable(Φ)))
        end
        last = findlast(!_is_zero_series, blocks)     # drop trailing zero blocks
        keep = last === nothing ? blocks[1:1] : blocks[1:last]
        new{T}(copy(keep), var)
    end
end

function LogSeries(blocks::Vector{FormalSeries{T}}, var::Symbol) where {T}
    LogSeries{T}(blocks, var)
end

# the var defaults to the first block's - hence the emptiness check comes first
function LogSeries(blocks::Vector{FormalSeries{T}}) where {T}
    isempty(blocks) && throw(InvalidArgument("a LogSeries needs ≥ 1 block"))
    LogSeries{T}(blocks, variable(blocks[1]))
end

"""
    LogSeries(Φ::FormalSeries) -> LogSeries

Embed a `FormalSeries` as the log-degree-0 case. Inverse of
[`FormalSeries(::LogSeries)`](@ref).
"""
LogSeries(Φ::FormalSeries{T}) where {T} = LogSeries{T}([Φ], variable(Φ))

"""
    FormalSeries(L::LogSeries) -> FormalSeries

Project a log-free `LogSeries` back to a `FormalSeries` (inverse of the embedding).
Throws `IncompatibleSeries` when `log_degree(L) > 0` - a genuine log block has no
`FormalSeries` representation.
"""
function FormalSeries(L::LogSeries)
    log_degree(L) == 0 || throw(IncompatibleSeries(:log_degree, log_degree(L), 0))
    L.blocks[1]
end

# -- accessors ---------------------------------------------------------------------

"""
    log_degree(L::LogSeries) -> Int

The degree `P` in `log(x)` (trailing zero blocks are dropped, so this is exact).
"""
log_degree(L::LogSeries) = length(L.blocks) - 1

"""
    log_blocks(L::LogSeries) -> Vector{FormalSeries}

The blocks ``[Φ^{[0]}, Φ^{[1]}, …]`` (a copy; series are immutable).
"""
log_blocks(L::LogSeries) = copy(L.blocks)

"""
    log_block(L::LogSeries, p) -> FormalSeries

The block ``Φ^{[p]}`` multiplying ``\\log(x)^p`` (0-based); `p > log_degree(L)` returns
the zero series, mirroring `Φ[n]` on [`FormalSeries`](@ref).
"""
function log_block(L::LogSeries, p::Integer)
    p < 0 && throw(InvalidArgument("log-block index must be ≥ 0, got $p"))
    p ≤ log_degree(L) ? L.blocks[p + 1] : _zero_like(L.blocks[1])
end

_zero_like(Φ::FormalSeries) = FormalSeries([zero(Φ.coeffs[1])], Φ.var)

variable(L::LogSeries) = L.var
n_terms(L::LogSeries) = maximum(n_terms, L.blocks)
is_exact(L::LogSeries) = all(is_exact, L.blocks)
Base.eltype(::Type{LogSeries{T}}) where {T} = T
_is_zero_series(L::LogSeries) = all(_is_zero_series, L.blocks)

Base.:(==)(L::LogSeries, M::LogSeries) = L.var == M.var && L.blocks == M.blocks

# -- arithmetic --------------------------------------------------------------------

_promote_series(L::LogSeries{T}, ::Type{Tp}) where {T,Tp} =
    T === Tp ? L : LogSeries([_promote_series(Φ, Tp) for Φ in L.blocks], L.var)

function Base.:+(L::LogSeries{T}, M::LogSeries{T}) where {T}
    L.var == M.var || throw(IncompatibleSeries(:var, L.var, M.var))
    P = max(log_degree(L), log_degree(M))
    LogSeries([log_block(L, p) + log_block(M, p) for p in 0:P], L.var)
end

Base.:-(L::LogSeries) = LogSeries([-Φ for Φ in L.blocks], L.var)
Base.:-(L::LogSeries, M::LogSeries) = L + (-M)

"""
    L * M

Product of two log-block polynomials: Cauchy in the log degree (degrees **add** - a log
block is exact, unlike the truncated `x` direction) with each block product the truncated
Cauchy product of [`FormalSeries`](@ref).
"""
function Base.:*(L::LogSeries{T}, M::LogSeries{T}) where {T}
    L.var == M.var || throw(IncompatibleSeries(:var, L.var, M.var))
    P = log_degree(L) + log_degree(M)
    blocks = map(0:P) do p
        terms = [log_block(L, i) * log_block(M, p - i)
                 for i in max(0, p - log_degree(M)):min(p, log_degree(L))]
        reduce(+, terms)
    end
    LogSeries(blocks, L.var)
end

Base.:*(λ::Number, L::LogSeries) = LogSeries([λ * Φ for Φ in L.blocks], L.var)
Base.:*(L::LogSeries, λ::Number) = λ * L

"""
    truncate(L::LogSeries; degree = nothing, terms = nothing) -> LogSeries

Truncate to log degree `≤ degree` (`nothing` keeps every block) with each block cut to
`terms` coefficients (`nothing` keeps each block's length).
"""
function Base.truncate(L::LogSeries; degree::Union{Nothing,Integer} = nothing,
                       terms::Union{Nothing,Integer} = nothing)
    P = degree === nothing ? log_degree(L) : min(Int(degree), log_degree(L))
    P ≥ 0 || throw(InvalidArgument("log-degree truncation must be ≥ 0, got $degree"))
    blocks = [terms === nothing ? L.blocks[p + 1] :
              truncate(L.blocks[p + 1], min(terms, n_terms(L.blocks[p + 1]))) for p in 0:P]
    LogSeries(blocks, L.var)
end

# -- calculus ----------------------------------------------------------------------

_shift_power(L::LogSeries, s::Integer) = LogSeries([_shift_power(Φ, s) for Φ in L.blocks],
                                                   L.var)

"""
    derivative(L::LogSeries) -> LogSeries

The derivative ``∂_x``, by Leibniz on each block:
``∂_x(\\log^p x · Φ) = \\log^p x · Φ' + p \\log^{p-1} x · x^{-1}Φ``.
"""
function derivative(L::LogSeries)
    P = log_degree(L)
    blocks = map(0:P) do p
        d = derivative(L.blocks[p + 1])
        p == P ? d : d + (p + 1) * _shift_power(L.blocks[p + 2], -1)
    end
    LogSeries(blocks, L.var)
end

"""
    evaluate(L::LogSeries, x; log_x = log(complex(x))) -> Number

Evaluate the truncated log-block polynomial at `x`, formally in each block (see
[`evaluate(::FormalSeries, x)`](@ref)).

`log_x` defaults to the **principal** branch. A transseries lives on a sector of the
universal cover of `x = 0`, so a continued `arg x` may leave `(-π, π]`; pass the
continued value explicitly in that case. The monodromy `log_x ↦ log_x + 2πi` mixes the
blocks - that is physical, not an artifact.
"""
function evaluate(L::LogSeries, x; log_x = log(complex(x)))
    total = evaluate(L.blocks[1], x) * one(log_x)
    for p in 1:log_degree(L)
        total += log_x^p * evaluate(L.blocks[p + 1], x)
    end
    total
end
