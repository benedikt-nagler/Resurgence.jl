# Makie weakdep extension: Borel-plane, large-order and optimal-truncation plots.
#
# Structure note (2026-08-10): every plot here is a Makie `@recipe`, not a
# Figure-returning function.  That is what gives each one a `!` form, a
# `fig[i, j]` target and Observable-reactive inputs, so the plots compose into a
# caller's layout and can be animated.  The exported `Resurgence.plot_*` names
# are thin forwards onto the generated recipe functions, which keeps the public
# API owned by `Resurgence` rather than by this extension module.
#
# Consequence of the recipe form: a recipe draws into a parent it does not own,
# so it cannot call `axislegend`.  Series are still labelled, so a caller gets a
# legend with one `axislegend(ax)`.  Axis labels/title/scale come back through
# `Makie.preferred_axis_attributes`, which the caller can override.
module ResurgenceMakieExt

using Makie: Makie, Figure, Axis, scatter!, lines!, hlines!, axislegend, @recipe
using Resurgence: Resurgence, BorelSeries, PadeApproximant, AbstractBorelApproximant,
                  FormalSeries, pade, poles, residues,
                  large_order_fit, optimal_truncation, partial_sums

const _SERIES_COLORS  = [:crimson, :dodgerblue, :seagreen, :darkorange, :purple]
const _SERIES_MARKERS = [:circle, :rect, :diamond, :utriangle, :cross]

_series_color(k)  = _SERIES_COLORS[mod1(k, length(_SERIES_COLORS))]
_series_marker(k) = _SERIES_MARKERS[mod1(k, length(_SERIES_MARKERS))]

# Anything a plot target may be: `fig[1, 2]`, an `Axis`, a `Scene`.
const _Target = Union{Makie.GridPosition, Makie.GridSubposition, Makie.Block, Makie.Scene}

# A short human label for an approximant, used in the legend / title.
_approx_label(r::PadeApproximant) =
    "[$(Resurgence.numerator_degree(r))/$(Resurgence.denominator_degree(r))] Padé"
_approx_label(r) = string(nameof(typeof(r)))

# Marker sizes proportional to residue magnitude (log-spaced, robust to a
# single pole or to all-equal magnitudes).  Returns a constant size when no
# residue data is usable.
function _marker_sizes(r; refine, smin = 8.0, smax = 26.0)
    mags = Float64[]
    try
        mags = [abs(ComplexF64(res)) for (_, res) in residues(r; refine)]
    catch
        return smin
    end
    (isempty(mags) || !all(isfinite, mags)) && return smin
    logs = log10.(mags .+ eps())
    lo, hi = extrema(logs)
    hi - lo < eps() && return (smin + smax) / 2
    smin .+ (smax - smin) .* (logs .- lo) ./ (hi - lo)
end

# ─── Borel plane ──────────────────────────────────────────────────────────────

"""
    plot_borel_plane(B::BorelSeries; order = nothing, rays = [0.0], …)
    plot_borel_plane(r::AbstractBorelApproximant; …)
    plot_borel_plane(rs::AbstractVector{<:AbstractBorelApproximant}; labels = nothing, …)

Scatter the Borel-plane poles of a Borel-Padé (or AAA) approximant in the
complex ζ-plane, with the Laplace rays at the angles in `rays` overlaid.

When `scale_by_residue` is set, marker area grows with `log₁₀|residue|`, so the
dominant singularities read as the largest dots - the standard "which pole
matters?" diagnostic.

Passing a *vector* of approximants overlays them in one axis (e.g. a Padé and an
AAA fit of the same series) with distinct colours and markers, for a
side-by-side comparison of where each method places the singularities.

This is a Makie recipe: `plot_borel_plane!(ax, r)` draws into an existing axis
and `plot_borel_plane(fig[1, 2], r)` into a layout slot. Call `axislegend(ax)`
for the legend.
"""
@recipe BorelPlane (approximants,) begin
    "Padé order used when the argument is a `BorelSeries` (`nothing` = default)."
    order = nothing
    "Laplace-ray angles (radians) drawn as dashed rays from the origin."
    rays = [0.0]
    "Refine pole locations before plotting."
    refine = true
    "Scale marker area with `log₁₀|residue|`."
    scale_by_residue = true
    "Per-approximant legend labels; `nothing` derives them from the types."
    labels = nothing
end

# Keep our own argument types intact - no coordinate conversion applies to them.
Makie.convert_arguments(::Type{<:BorelPlane}, x) = (x,)

@static if isdefined(Makie, :preferred_axis_attributes)
    Makie.preferred_axis_attributes(::Type{Axis}, ::BorelPlane) =
        (; xlabel = "Re ζ", ylabel = "Im ζ", title = "Borel plane: poles")
end

function _draw_rays!(p, ζs, rays)
    R = isempty(ζs) ? 1.0 : 1.5 * maximum(abs, ζs)
    for θ in rays
        lines!(p, [0.0, R * cos(θ)], [0.0, R * sin(θ)];
               color = :gray, linestyle = :dash)
    end
end

# Core: draw a fixed list of approximants.  The *number* of series is fixed when
# the plot is created; pole positions stay reactive through `refine`.
function _draw_borel!(p, rs, labs)
    refine = p.refine[]
    allζ = ComplexF64[]
    for r in rs
        append!(allζ, ComplexF64.(poles(r; refine)))
    end
    _draw_rays!(p, allζ, p.rays[])
    for (k, r) in enumerate(rs)
        ζs = ComplexF64.(poles(r; refine))
        isempty(ζs) && continue
        ms = p.scale_by_residue[] ? _marker_sizes(r; refine) : 10.0
        scatter!(p, real.(ζs), imag.(ζs);
                 color = _series_color(k), marker = _series_marker(k),
                 markersize = ms, label = string(labs[k]))
    end
    p
end

function Makie.plot!(p::BorelPlane{<:Tuple{<:AbstractVector{<:AbstractBorelApproximant}}})
    rs = p.approximants[]
    isempty(rs) && throw(ArgumentError("plot_borel_plane: no approximants to plot"))
    labs = p.labels[] === nothing ? [_approx_label(r) for r in rs] : p.labels[]
    _draw_borel!(p, rs, labs)
end

function Makie.plot!(p::BorelPlane{<:Tuple{<:AbstractBorelApproximant}})
    r = p.approximants[]
    labs = p.labels[] === nothing ? [_approx_label(r)] : p.labels[]
    _draw_borel!(p, [r], labs)
end

function Makie.plot!(p::BorelPlane{<:Tuple{<:BorelSeries}})
    r = pade(p.approximants[]; order = p.order[])
    labs = p.labels[] === nothing ? [_approx_label(r)] : p.labels[]
    _draw_borel!(p, [r], labs)
end

# ─── Large-order growth (ratio test) ──────────────────────────────────────────

"""
    plot_large_order(Φ::FormalSeries; asymptote = true)

Plot the ratio-test diagnostic `|aₙ/aₙ₋₁| / n` against `n`.  For a series with
factorial growth `aₙ ~ S · Γ(n+b) / Aⁿ` this tends to `1/A`, the reciprocal
distance to the leading Borel singularity.  When `asymptote` is set, the fitted
`1/A` from [`large_order_fit`](@ref) is drawn as a dashed line.

This is a Makie recipe: see [`plot_borel_plane`](@ref) for the `!` and
layout-slot forms.
"""
@recipe LargeOrder (series,) begin
    "Draw the fitted `1/A` asymptote."
    asymptote = true
end

Makie.convert_arguments(::Type{<:LargeOrder}, x) = (x,)

@static if isdefined(Makie, :preferred_axis_attributes)
    Makie.preferred_axis_attributes(::Type{Axis}, ::LargeOrder) =
        (; xlabel = "n", ylabel = "|aₙ/aₙ₋₁| / n", title = "Large-order ratio test")
end

function Makie.plot!(p::LargeOrder{<:Tuple{<:FormalSeries}})
    Φ = p.series[]
    c = Float64.(Φ.coeffs)
    ns = Int[]
    diag = Float64[]
    for n in 2:length(c)
        (iszero(c[n - 1]) || !isfinite(c[n])) && continue
        push!(ns, n - 1)                       # ratio uses index n-1 (1-based)
        push!(diag, abs(c[n] / c[n - 1]) / (n - 1))
    end
    isempty(ns) || scatter!(p, ns, diag; color = :crimson, label = "ratio")
    if p.asymptote[]
        try
            A = Float64(abs(large_order_fit(Φ).A))
            A > 0 && hlines!(p, [1 / A];
                             color = :black, linestyle = :dash,
                             label = "1/A = $(round(1 / A; sigdigits = 4))")
        catch
        end
    end
    p
end

# ─── Optimal truncation (superasymptotic U-curve) ─────────────────────────────

# Term magnitudes |aₙ ħ^{pₙ}| reconstructed from the partial sums, and the
# optimal truncation order N⋆ (index of the smallest term).
function _truncation_terms(Φ::FormalSeries, ħ)
    ps = partial_sums(Φ, ħ)
    terms = abs.(diff(vcat(zero(first(ps)), ps)))
    Nstar = try
        optimal_truncation(Φ, ħ).N
    catch
        argmin(terms) - 1
    end
    Float64.(terms), Int(Nstar)
end

"""
    plot_optimal_truncation(Φ::FormalSeries, ħ::Number)
    plot_optimal_truncation(Φ::FormalSeries, ħs::AbstractVector)

Plot the superasymptotic error curve: term magnitudes `|aₙ ħ^{pₙ}|` against the
truncation order `N` on a logarithmic `y`-axis.  The U-shape - terms shrinking,
then diverging factorially - is the hallmark of an asymptotic series; its
minimum is the optimal truncation order `N⋆` from [`optimal_truncation`](@ref),
marked here.  Passing several couplings overlays one U-curve per `ħ`, showing
`N⋆` grow and the least achievable error shrink as `ħ → 0`.

This is a Makie recipe: see [`plot_borel_plane`](@ref) for the `!` and
layout-slot forms.  The log `y`-scale arrives as a preferred axis attribute, so
a caller drawing into their own `Axis` must set `yscale = log10` themselves.
"""
@recipe OptimalTruncation (series, couplings) begin end

Makie.convert_arguments(::Type{<:OptimalTruncation}, Φ, ħ::Number) = (Φ, [ħ])
Makie.convert_arguments(::Type{<:OptimalTruncation}, Φ, ħs) = (Φ, ħs)

@static if isdefined(Makie, :preferred_axis_attributes)
    Makie.preferred_axis_attributes(::Type{Axis}, ::OptimalTruncation) =
        (; xlabel = "truncation order N", ylabel = "|aₙ ħ^{pₙ}|", yscale = log10,
           title = "Optimal truncation (superasymptotic U-curve)")
end

function Makie.plot!(p::OptimalTruncation{<:Tuple{<:FormalSeries,<:AbstractVector}})
    Φ, ħs = p.series[], p.couplings[]
    isempty(ħs) && throw(ArgumentError("plot_optimal_truncation: no couplings given"))
    for (k, ħ) in enumerate(ħs)
        terms, Nstar = _truncation_terms(Φ, ħ)
        # log scale: keep strictly-positive terms
        idx = [n for n in eachindex(terms) if terms[n] > 0]
        isempty(idx) && continue
        col = _series_color(k)
        lines!(p, idx .- 1, terms[idx]; color = col, label = "ħ = $ħ")
        if 0 <= Nstar <= length(terms) - 1 && terms[Nstar + 1] > 0
            scatter!(p, [Nstar], [terms[Nstar + 1]];
                     color = col, marker = :star5, markersize = 18)
        end
    end
    p
end

# ─── Public names: forward `Resurgence.plot_*` onto the recipes ───────────────

const _BorelArg = Union{BorelSeries, AbstractBorelApproximant,
                        AbstractVector{<:AbstractBorelApproximant}}

Resurgence.plot_borel_plane(a::_BorelArg; kw...) = borelplane(a; kw...)
Resurgence.plot_borel_plane(t::_Target, a::_BorelArg; kw...) = borelplane(t, a; kw...)
Resurgence.plot_borel_plane!(a::_BorelArg; kw...) = borelplane!(a; kw...)
Resurgence.plot_borel_plane!(t::_Target, a::_BorelArg; kw...) = borelplane!(t, a; kw...)

Resurgence.plot_large_order(Φ::FormalSeries; kw...) = largeorder(Φ; kw...)
Resurgence.plot_large_order(t::_Target, Φ::FormalSeries; kw...) = largeorder(t, Φ; kw...)
Resurgence.plot_large_order!(Φ::FormalSeries; kw...) = largeorder!(Φ; kw...)
Resurgence.plot_large_order!(t::_Target, Φ::FormalSeries; kw...) = largeorder!(t, Φ; kw...)

Resurgence.plot_optimal_truncation(Φ::FormalSeries, ħ; kw...) =
    optimaltruncation(Φ, ħ; kw...)
Resurgence.plot_optimal_truncation(t::_Target, Φ::FormalSeries, ħ; kw...) =
    optimaltruncation(t, Φ, ħ; kw...)
Resurgence.plot_optimal_truncation!(Φ::FormalSeries, ħ; kw...) =
    optimaltruncation!(Φ, ħ; kw...)
Resurgence.plot_optimal_truncation!(t::_Target, Φ::FormalSeries, ħ; kw...) =
    optimaltruncation!(t, Φ, ħ; kw...)

end # module ResurgenceMakieExt
