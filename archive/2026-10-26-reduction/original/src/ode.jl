# Sector-recursion Newton solver for transseries solutions of nonlinear ODEs. The driver
# here is model-agnostic: it walks the instanton sectors in order, and at each sector calls
# a `solve_sector(sofar, k)` closure that returns Φ_k - the fluctuation of the k-th
# transmonomial e^{-kA/ħ}. The model (its perturbative seed, its linearized operator, and
# the quadratic source assembled from lower sectors) lives in that closure; the canonical
# instance is Painlevé I, see `painleve.jl`.
#
# "Newton" is literal: substituting u = Σ_k C^k e^{-kA…} Φ_k into a nonlinear ODE and
# grading by C^k gives, for each k ≥ 1, one *linear* solve of the ODE's Jacobian
# (the linearization about the perturbative solution) against a source built from the
# already-known lower sectors - exactly a Newton step per sector.

"""
    transseries_solve(seed, action, solve_sector; sectors = 2, var = variable(seed))
        -> Transseries

Assemble a one-parameter transseries ``u = \\sum_{k≥0} C^k e^{-kA/ħ}\\, Φ_k`` from its
perturbative seed ``Φ_0 =`` `seed` (a [`FormalSeries`](@ref)) by the sector recursion: for
each `k = 1, …, sectors-1`, `Φ_k = solve_sector(Φ₀…Φ_{k-1}, k)`, where `solve_sector` is a
closure `(sofar::Vector{FormalSeries}, k::Int) -> FormalSeries` that performs the model's
per-sector linear (Newton) solve. The result is a [`Transseries`](@ref) with the given
`action`.

This is the model-agnostic driver; the per-sector linear solve - the ODE's linearized
operator and the quadratic source assembled from `sofar` - is the caller's closure. See
[`painleve1`](@ref) for the Painlevé I instance and [`_solve_linear_sector`](@ref) for the
order-by-order primitive most closures build on.
"""
function transseries_solve(seed::FormalSeries{T}, action::Number, solve_sector;
                           sectors::Integer = 2,
                           var::Symbol = variable(seed)) where {T}
    sectors ≥ 1 || throw(InvalidArgument("need sectors ≥ 1, got $sectors"))
    secs = FormalSeries{T}[seed]
    for k in 1:(sectors - 1)
        Φ = solve_sector(secs, k)
        variable(Φ) == var ||
            throw(TransseriesSolveError("sector $k came back in variable $(variable(Φ)), " *
                                        "expected $var"))
        push!(secs, Φ)
    end
    Transseries(action, secs, var)
end

"""
    multi_transseries_solve(seed, actions, solve_sector; orders = 2, var = variable(seed))
        -> MultiTransseries

The `K`-parameter generalization of [`transseries_solve`](@ref): assemble
``u = \\sum_{n ∈ ℤ_{≥0}^K} (\\prod_i σ_i^{n_i}) e^{-(n·A)/ħ}\\, Φ_n`` from its perturbative
seed `Φ_0 =` `seed` over the action lattice `actions`, for every multi-index of total
degree `1 … orders`. `solve_sector` is a closure
`(sofar::Dict{NTuple{K,Int}}, n::NTuple{K,Int}) -> AbstractSeries` returning `Φ_n`; it is
called in order of increasing `sum(n)`, so every strictly smaller multi-index is already
in `sofar` (an absent key is an exact zero).

Unlike the rank-1 driver this admits **resonance**: when the actions are ℤ-dependent,
distinct multi-indices share an exponential weight and sectors acquire `log` blocks (see
[`resonance_lattice`](@ref)). Return a [`LogSeries`](@ref) from `solve_sector` in that
case - [`_solve_log_sector`](@ref) is the order-by-order primitive that produces one - and
the assembled [`MultiTransseries`](@ref) carries the tower. The canonical instance is
[`painleve1_two_parameter`](@ref).
"""
function multi_transseries_solve(seed::AbstractSeries, actions::Tuple, solve_sector;
                                 orders::Integer = 2,
                                 var::Symbol = variable(seed))
    orders ≥ 0 || throw(InvalidArgument("need orders ≥ 0, got $orders"))
    K = length(actions)
    K ≥ 1 || throw(InvalidArgument("need rank K ≥ 1, got 0 actions"))
    sectors = Dict{NTuple{K,Int},typeof(seed)}(ntuple(_ -> 0, K) => seed)
    for total in 1:orders, idx in _multi_indices(K, total)
        Φ = solve_sector(sectors, idx)
        variable(Φ) == var ||
            throw(TransseriesSolveError("sector $idx came back in variable " *
                                        "$(variable(Φ)), expected $var"))
        sectors[idx] = Φ
    end
    MultiTransseries(actions, sectors; var)
end

# every NTuple{K,Int} with nonnegative entries summing to `total`, lexicographic
function _multi_indices(K::Int, total::Int)
    K == 1 && return [(total,)]
    out = NTuple{K,Int}[]
    for first in total:-1:0, rest in _multi_indices(K - 1, total - first)
        push!(out, (first, rest...))
    end
    out
end

"""
    _solve_linear_sector(apply_L, g, β, δ, N; free = Int[], freeval = 1) -> FormalSeries

Order-by-order solution of a linear sector equation `apply_L(Φ) = g` for the fluctuation
`Φ = Σ_m b_m s^{β+m}` on the unit grid (length `N`, leading power `β`). `apply_L` is the
sector's linearized operator (a `FormalSeries -> FormalSeries` closure). `δ` is the
operator's diagonal shift - `apply_L` sends a coefficient at power `p` to a leading output
at `p + δ` - so the coefficient at grid index `m` is fixed by the residual at power
`β + m + δ`. Grid positions whose probed diagonal vanishes are skipped (this absorbs any
step > 1 sparsity automatically), and indices in `free` are integration constants: set to
`freeval` and never solved (the homogeneous mode carrying the transseries parameter).
"""
function _solve_linear_sector(apply_L, g::FormalSeries{T}, β::Rational{Int},
                              δ::Integer, N::Integer;
                              free = Int[], freeval = one(T)) where {T}
    c = zeros(T, N)
    for i in free
        0 ≤ i < N && (c[i + 1] = T(freeval))
    end
    Φ = FormalSeries(c, variable(g); power_offset = β)
    unit(m) = (e = zeros(T, N); e[m + 1] = one(T);
               FormalSeries(e, variable(g); power_offset = β))
    for m in 0:(N - 1)
        m in free && continue
        pout = β + m + δ
        d = _coeff_at(apply_L(unit(m)), pout)
        iszero(d) && continue
        r = _coeff_at(g - apply_L(Φ), pout)
        c[m + 1] += r / d
        Φ = FormalSeries(c, variable(g); power_offset = β)
    end
    Φ
end

# -- the resonant (log-tower) sector solve ------------------------------------------

# A sector operator 𝓛 acting on a log block picks up commutators, because ∂ hits log x:
# with L = log x and ψ = w L^p,
#   𝓛(w L^p) = L^p 𝓛(w) + p L^{p-1} C(w) + p(p-1) L^{p-2} D(w),
# where C and D are the model's first- and second-order log commutators. Collecting L^p,
# the block equations of 𝓛ψ = g are, descending in p,
#   𝓛ψ^{[p]} = g^{[p]} - (p+1) C(ψ^{[p+1]}) - (p+2)(p+1) D(ψ^{[p+2]}),
# which is why the tower is solved from the top down. Growth happens at the *kernel*: when
# 𝓛 has a homogeneous solution on the grid, its position carries a vanishing diagonal, and
# a nonzero residual there cannot be solved - it is instead cancelled by the C-term, which
# fixes the free multiple of the homogeneous solution in the block *above*. That is the
# same move `resonant_solve` makes for the channel equation (ħ²∂ + c)u = g.

"""
    _solve_log_sector(apply_L, apply_C, apply_D, g, β, δ, N; kernel = nothing,
                      homog = nothing, freeval = 0) -> LogSeries

Order-by-order solution of a resonant sector equation `apply_L(ψ) = g` for a `ψ` carrying
a `log` tower, on the unit grid of length `N` with leading power `β` and diagonal shift
`δ` (as in [`_solve_linear_sector`](@ref), which solves each individual log block).
`apply_C` and `apply_D` are the operator's first- and second-order log commutators; pass
`apply_D = nothing` when the operator is first order.

`kernel` is the grid index at which `apply_L`'s diagonal vanishes and `homog` the
corresponding homogeneous solution normalized to `1` there; with `kernel = nothing` the
solve cannot grow the tower and only propagates the logs already present in `g`. Otherwise
an unsolvable residual at `kernel` fixes the free multiple of `homog` one log block up -
so the returned degree is `log_degree(g)` or `log_degree(g) + 1`. `freeval` is the value
left at `kernel` in block `0`, the sector's own integration constant (the transseries
parameter, hence `0` in every sector but the one that defines it).
"""
function _solve_log_sector(apply_L, apply_C, apply_D, g::LogSeries{T},
                           β::Rational{Int}, δ::Integer, N::Integer;
                           kernel::Union{Nothing,Integer} = nothing,
                           homog::Union{Nothing,FormalSeries{T}} = nothing,
                           freeval = zero(T)) where {T}
    kernel === nothing && homog !== nothing &&
        throw(InvalidArgument("_solve_log_sector needs a `kernel` index for `homog`"))
    var = variable(g)
    P = log_degree(g)
    free = kernel === nothing ? Int[] : [Int(kernel)]
    # the source and every block equation live on the output grid, at power β + δ. An
    # absent log block comes back as a bare zero whose offset need not sit on that grid.
    grid = FormalSeries(zeros(T, N), var; power_offset = β + δ)
    pad = Φ -> _is_zero_series(Φ) ? grid : grid + Φ
    # blocks[p+1] = ψ^{[p]}; two guard blocks above the top keep the recursion total
    zeroblock = FormalSeries(zeros(T, N), var; power_offset = β)
    blocks = FormalSeries{T}[zeroblock for _ in 1:(P + 4)]      # immutable, so shared
    rhs_of = p -> begin
        r = pad(log_block(g, p))
        r = r - (p + 1) * pad(apply_C(blocks[p + 2]))
        apply_D === nothing ? r :
        r - ((p + 2) * (p + 1)) * pad(apply_D(blocks[p + 3]))
    end
    for p in (P + 1):-1:0
        if kernel !== nothing && p ≤ P
            # the residual at the vanishing diagonal fixes block p+1's free constant
            pout = β + kernel + δ
            r = _coeff_at(rhs_of(p), pout)
            if !iszero(r)
                homog === nothing &&
                    throw(TransseriesSolveError("log block $p is obstructed at power " *
                                                "$pout with residual $r, but no " *
                                                "homogeneous solution was supplied"))
                ch = _coeff_at(apply_C(homog), pout)
                iszero(ch) &&
                    throw(TransseriesSolveError("log block $p is obstructed at power " *
                                                "$pout with residual $r, but the " *
                                                "commutator cannot absorb it"))
                blocks[p + 2] = blocks[p + 2] + (r / ((p + 1) * ch)) * homog
            end
        end
        blocks[p + 1] = _solve_linear_sector(apply_L, rhs_of(p), β, δ, N;
                                             free, freeval = p == 0 ? freeval : zero(T))
    end
    LogSeries(blocks[1:(P + 2)], var)
end
