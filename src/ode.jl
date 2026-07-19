# Sector-recursion Newton solver for transseries solutions of nonlinear ODEs. The driver
# here is model-agnostic: it walks the instanton sectors in order, and at each sector calls
# a `solve_sector(sofar, k)` closure that returns Φ_k — the fluctuation of the k-th
# transmonomial e^{-kA/ħ}. The model (its perturbative seed, its linearized operator, and
# the quadratic source assembled from lower sectors) lives in that closure; the canonical
# instance is Painlevé I, see `painleve.jl`.
#
# "Newton" is literal: substituting u = Σ_k C^k e^{-kA…} Φ_k into a nonlinear ODE and
# grading by C^k gives, for each k ≥ 1, one *linear* solve of the ODE's Jacobian
# (the linearization about the perturbative solution) against a source built from the
# already-known lower sectors — exactly a Newton step per sector.

"""
    transseries_solve(seed, action, solve_sector; sectors = 2, var = variable(seed))
        -> Transseries

Assemble a one-parameter transseries ``u = \\sum_{k≥0} C^k e^{-kA/ħ}\\, Φ_k`` from its
perturbative seed ``Φ_0 =`` `seed` (a [`FormalSeries`](@ref)) by the sector recursion: for
each `k = 1, …, sectors-1`, `Φ_k = solve_sector(Φ₀…Φ_{k-1}, k)`, where `solve_sector` is a
closure `(sofar::Vector{FormalSeries}, k::Int) -> FormalSeries` that performs the model's
per-sector linear (Newton) solve. The result is a [`Transseries`](@ref) with the given
`action`.

This is the model-agnostic driver; the per-sector linear solve — the ODE's linearized
operator and the quadratic source assembled from `sofar` — is the caller's closure. See
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
    _solve_linear_sector(apply_L, g, β, δ, N; free = Int[], freeval = 1) -> FormalSeries

Order-by-order solution of a linear sector equation `apply_L(Φ) = g` for the fluctuation
`Φ = Σ_m b_m s^{β+m}` on the unit grid (length `N`, leading power `β`). `apply_L` is the
sector's linearized operator (a `FormalSeries -> FormalSeries` closure). `δ` is the
operator's diagonal shift — `apply_L` sends a coefficient at power `p` to a leading output
at `p + δ` — so the coefficient at grid index `m` is fixed by the residual at power
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
