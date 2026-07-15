# Makie weakdep extension: Borel-plane pole scatter with Laplace-ray overlay.
module ResurgenceMakieExt

using Makie: Makie, Figure, Axis, scatter!, lines!, hlines!, axislegend
using Resurgence: Resurgence, BorelSeries, PadeApproximant, AbstractBorelApproximant,
                  FormalSeries, pade, poles, residues,
                  large_order_fit, optimal_truncation, partial_sums

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

function _draw_rays!(ax, ζs, rays)
    R = isempty(ζs) ? 1.0 : 1.5 * maximum(abs, ζs)
    for θ in rays
        lines!(ax, [0.0, R * cos(θ)], [0.0, R * sin(θ)];
               color = :gray, linestyle = :dash)
    end
end

"""
    plot_borel_plane(B::BorelSeries; order = nothing, rays = [0.0],
                     refine = true, scale_by_residue = true)
    plot_borel_plane(r::AbstractBorelApproximant; rays = [0.0], refine = true,
                     scale_by_residue = true)
    plot_borel_plane(rs::AbstractVector{<:AbstractBorelApproximant};
                     rays = [0.0], refine = true, labels = nothing)

Scatter the Borel-plane poles of a Borel–Padé (or AAA) approximant in the
complex ζ-plane, with the Laplace rays at the angles in `rays` overlaid.

When `scale_by_residue` is set, marker area grows with `log₁₀|residue|`, so the
dominant singularities read as the largest dots - the standard "which pole
matters?" diagnostic.

Passing a *vector* of approximants overlays them in one axis (e.g. a Padé and an
AAA fit of the same series) with distinct colours, for a side-by-side comparison
of where each method places the singularities.

Returns a `Makie.Figure`.
"""
function Resurgence.plot_borel_plane(B::BorelSeries;
                                     order::Union{Nothing,Integer} = nothing,
                                     kwargs...)
    Resurgence.plot_borel_plane(pade(B; order); kwargs...)
end

function Resurgence.plot_borel_plane(r::AbstractBorelApproximant;
                                     rays::AbstractVector{<:Real} = [0.0],
                                     refine::Bool = true,
                                     scale_by_residue::Bool = true)
    ζs = ComplexF64.(poles(r; refine))
    fig = Figure()
    ax = Axis(fig[1, 1]; xlabel = "Re ζ", ylabel = "Im ζ",
              title = "Borel plane: $(_approx_label(r)) poles")
    _draw_rays!(ax, ζs, rays)
    if !isempty(ζs)
        ms = scale_by_residue ? _marker_sizes(r; refine) : 10.0
        scatter!(ax, real.(ζs), imag.(ζs);
                 color = :crimson, markersize = ms, label = _approx_label(r))
        axislegend(ax)
    end
    fig
end

function Resurgence.plot_borel_plane(rs::AbstractVector{<:AbstractBorelApproximant};
                                     rays::AbstractVector{<:Real} = [0.0],
                                     refine::Bool = true,
                                     labels::Union{Nothing,AbstractVector} = nothing,
                                     scale_by_residue::Bool = true)
    isempty(rs) && throw(ArgumentError("plot_borel_plane: no approximants to plot"))
    colors  = [:crimson, :dodgerblue, :seagreen, :darkorange, :purple]
    markers = [:circle, :rect, :diamond, :utriangle, :cross]
    labs = labels === nothing ? [_approx_label(r) for r in rs] : labels
    allζ = ComplexF64[]
    for r in rs
        append!(allζ, ComplexF64.(poles(r; refine)))
    end
    fig = Figure()
    ax = Axis(fig[1, 1]; xlabel = "Re ζ", ylabel = "Im ζ",
              title = "Borel plane: approximant comparison")
    _draw_rays!(ax, allζ, rays)
    for (k, r) in enumerate(rs)
        ζs = ComplexF64.(poles(r; refine))
        isempty(ζs) && continue
        ms = scale_by_residue ? _marker_sizes(r; refine) : 10.0
        scatter!(ax, real.(ζs), imag.(ζs);
                 color = colors[mod1(k, length(colors))],
                 marker = markers[mod1(k, length(markers))],
                 markersize = ms, label = string(labs[k]))
    end
    axislegend(ax)
    fig
end

# ─── Large-order growth (ratio test) ──────────────────────────────────────────

"""
    plot_large_order(Φ::FormalSeries; asymptote = true)

Plot the ratio-test diagnostic `|aₙ/aₙ₋₁| / n` against `n`.  For a series with
factorial growth `aₙ ~ S · Γ(n+b) / Aⁿ` this tends to `1/A`, the reciprocal
distance to the leading Borel singularity.  When `asymptote` is set, the fitted
`1/A` from [`large_order_fit`](@ref) is drawn as a dashed line.

Returns a `Makie.Figure`.
"""
function Resurgence.plot_large_order(Φ::FormalSeries; asymptote::Bool = true)
    c = Float64.(Φ.coeffs)
    ns = Int[]
    diag = Float64[]
    for n in 2:length(c)
        (iszero(c[n - 1]) || !isfinite(c[n])) && continue
        push!(ns, n - 1)                       # ratio uses index n-1 (1-based)
        push!(diag, abs(c[n] / c[n - 1]) / (n - 1))
    end
    fig = Figure()
    ax = Axis(fig[1, 1]; xlabel = "n", ylabel = "|aₙ/aₙ₋₁| / n",
              title = "Large-order ratio test")
    isempty(ns) || scatter!(ax, ns, diag; color = :crimson, label = "ratio")
    if asymptote
        try
            A = Float64(abs(large_order_fit(Φ).A))
            A > 0 && hlines!(ax, [1 / A];
                             color = :black, linestyle = :dash,
                             label = "1/A = $(round(1 / A; sigdigits = 4))")
        catch
        end
    end
    axislegend(ax; position = :rt)
    fig
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

Returns a `Makie.Figure`.
"""
Resurgence.plot_optimal_truncation(Φ::FormalSeries, ħ::Number; kwargs...) =
    Resurgence.plot_optimal_truncation(Φ, [ħ]; kwargs...)

function Resurgence.plot_optimal_truncation(Φ::FormalSeries, ħs::AbstractVector)
    isempty(ħs) && throw(ArgumentError("plot_optimal_truncation: no couplings given"))
    colors = [:crimson, :dodgerblue, :seagreen, :darkorange, :purple]
    fig = Figure()
    ax = Axis(fig[1, 1]; xlabel = "truncation order N",
              ylabel = "|aₙ ħ^{pₙ}|", yscale = log10,
              title = "Optimal truncation (superasymptotic U-curve)")
    for (k, ħ) in enumerate(ħs)
        terms, Nstar = _truncation_terms(Φ, ħ)
        # log scale: keep strictly-positive terms
        idx = [n for n in eachindex(terms) if terms[n] > 0]
        isempty(idx) && continue
        col = colors[mod1(k, length(colors))]
        lines!(ax, idx .- 1, terms[idx]; color = col, label = "ħ = $ħ")
        if 0 <= Nstar <= length(terms) - 1 && terms[Nstar + 1] > 0
            scatter!(ax, [Nstar], [terms[Nstar + 1]];
                     color = col, marker = :star5, markersize = 18)
        end
    end
    axislegend(ax; position = :lb)
    fig
end

end # module ResurgenceMakieExt
