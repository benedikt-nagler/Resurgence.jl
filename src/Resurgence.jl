"""
    Resurgence

Foundation package for resurgence theory: generic formal power series, exact Borel
transforms, Padé approximants over any field, Borel–Padé–Laplace summation with
lateral sums and Stokes discontinuities, and large-order analysis. Deliberately free
of any cluster-algebra or WKB dependency — that bridge lives in `ExactWKB.jl`.
"""
module Resurgence

using LinearAlgebra: LinearAlgebra
using QuadGK: quadgk
using PrecompileTools: @setup_workload, @compile_workload

export ResurgenceError, InvalidArgument, IncompatibleSeries, DegeneratePade, PoleOnRay
export NoSingularityFound
export AbstractSeries, FormalSeries, BorelSeries, PadeApproximant
export AbstractBorelApproximant
export n_terms, coefficients, variable, power_offset, is_exact, evaluate
export borel, inverse_borel
export pade, poles, residues, borel_pade_poles, numerator_degree, denominator_degree
export laplace_sum, lateral_sum, stokes_discontinuity, borel_sum
export richardson, large_order_fit, accelerate, partial_sums
export ConformalMap, ConformalPade, conformal_map, conformal_borel, inverse
export subtract_singularity
export dingle_terminant, optimal_truncation, hyper_sum
export Transseries, action, n_sectors, sector, sectors, transseries_sum
export stokes_constant, alien_derivative, stokes_automorphism, median_sum
export plot_borel_plane, aaa_approximant, aaa_borel

include("errors.jl")
include("formal_series.jl")
include("transseries.jl")
include("named_series.jl")
include("borel.jl")
include("pade.jl")
include("laplace.jl")
include("large_order.jl")
include("acceleration.jl")
include("conformal.jl")
include("singularities.jl")
include("hyperasymptotics.jl")
include("alien.jl")
include("show.jl")

"""
    plot_borel_plane(B::BorelSeries; order = nothing, rays = [0.0], kwargs...)

Scatter the Borel–Padé poles of `B` in the complex ζ-plane with the Laplace rays
overlaid. Implemented in the Makie package extension — load a Makie backend
(e.g. `using CairoMakie`) first.
"""
function plot_borel_plane(args...; kwargs...)
    if isempty(args) || args[1] isa Union{BorelSeries,PadeApproximant}
        error("plot_borel_plane requires a Makie backend: run `using CairoMakie` " *
              "(or GLMakie) and retry")
    else
        throw(MethodError(plot_borel_plane, args))
    end
end

"""
    aaa_approximant(zs, fs; tol = 1e-13, mmax = 100, var = :ζ) -> AAAApproximant
    aaa_borel(B::BorelSeries; samples = 128, radius = nothing, tol = 1e-13,
              mmax = 100) -> AAAApproximant

AAA rational approximation of the reduced Borel function — an alternative
pole/branch-cut hunter, robust where equispaced-coefficient Padé struggles.
Implemented in the BaryRational package extension: run `using BaryRational` first.
`aaa_approximant` fits samples `(zs, fs)` directly; `aaa_borel` samples the
truncated reduced Borel sum of `B` on disks inside its radius of convergence. The
result is an [`AbstractBorelApproximant`](@ref): it plugs into
`laplace_sum(B, aaa_borel(B), ħ)` and supports [`poles`](@ref)/[`residues`](@ref).
"""
function aaa_approximant(args...; kwargs...)
    error("aaa_approximant requires BaryRational: run `using BaryRational` and retry")
end

@doc (@doc aaa_approximant)
function aaa_borel(args...; kwargs...)
    if isempty(args) || args[1] isa BorelSeries
        error("aaa_borel requires BaryRational: run `using BaryRational` and retry")
    else
        throw(MethodError(aaa_borel, args))
    end
end

@setup_workload begin
    @compile_workload begin
        Φ = FormalSeries(:euler, 12)
        B = borel(Φ)
        r = pade(B; order = 1)   # the [n/n] Padé of the rational Euler Borel
        poles(r)                 # function is degenerate for n ≥ 2 …
        pade(B; order = 3, reduce = true)   # … unless reduction is requested
        c = conformal_borel(B; zeta0 = -1 // 1)
        poles(c)
        accelerate([1.0, 0.5, 0.75, 0.625, 0.6875]; method = :wynn)
        Φ + Φ
        Φ * Φ
        sprint(show, Φ)
        F = Transseries(:euler, 12)
        alien_derivative(F; stokes = -2 * BigFloat(π) * im)
        F + F
        sprint(show, F)
    end
end

end # module Resurgence
