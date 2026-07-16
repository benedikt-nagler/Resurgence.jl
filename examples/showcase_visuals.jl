# ═══════════════════════════════════════════════════════════════════════════
#  Resurgence.jl — visual showcase
#
#  Renders every plotting function in the package to PNG using the headless
#  CairoMakie backend.  For live figures swap CairoMakie for GLMakie
#  (`using GLMakie`) and replace the `save(...)` calls with `display(...)`.
#
#  Run (from the package root), using the examples environment:
#
#      julia --project=examples examples/showcase_visuals.jl
#
#  or the test environment, which also has CairoMakie + BaryRational resolved:
#
#      julia --project=test examples/showcase_visuals.jl
#
#  Figures are written to examples/figures/.
# ═══════════════════════════════════════════════════════════════════════════

using Resurgence
using CairoMakie
using BaryRational          # enables the AAA approximant used in the overlay

const OUT = joinpath(@__DIR__, "figures")
mkpath(OUT)
save_fig(name, fig) = (path = joinpath(OUT, name); CairoMakie.save(path, fig);
                       println("  wrote $path"))

println("Resurgence.jl visual showcase → $OUT")

# ── 1. Borel plane: Borel–Padé poles + Laplace ray ─────────────────────────────
# Marker area ∝ log|residue|, so the leading singularity reads as the largest
# dot.  Euler series: single Borel pole at ζ = −1.  (`reduce = true` degrades the
# degenerate Euler [m/m] system gracefully.)
println("\n[1] Borel plane (plot_borel_plane)")
B = borel(FormalSeries(:euler, 16))
r = pade(B; order = 6, reduce = true)
save_fig("01a_borel_euler.png", plot_borel_plane(r; rays = [Float64(π)]))

# Airy: Borel singularity at ζ = −4/3.
B_airy = borel(FormalSeries(:airy, 16))
save_fig("01b_borel_airy.png", plot_borel_plane(pade(B_airy; order = 6, reduce = true)))

# ── 2. Padé-vs-AAA overlay ─────────────────────────────────────────────────────
# The two rational approximants place the Borel singularity in the same spot;
# spurious poles differ.  Distinct colours/markers per method.
println("\n[2] Padé-vs-AAA overlay (plot_borel_plane, vector form)")
aaa = aaa_borel(B)
save_fig("02_borel_pade_vs_aaa.png",
         plot_borel_plane(AbstractBorelApproximant[r, aaa]; labels = ["Padé", "AAA"]))

# ── 3. Large-order ratio test ──────────────────────────────────────────────────
# |aₙ/aₙ₋₁| / n → 1/A.  Dashed line = fitted 1/A from large_order_fit.
println("\n[3] large-order ratio test (plot_large_order)")
save_fig("03a_large_order_euler.png", plot_large_order(FormalSeries(:euler, 24)))
save_fig("03b_large_order_airy.png", plot_large_order(FormalSeries(:airy, 24)))

# ── 4. Optimal truncation (superasymptotic U-curve) ────────────────────────────
# Term magnitudes vs truncation order N; the minimum is N⋆ (starred).  Smaller
# coupling ⇒ later N⋆ and smaller least-error.
println("\n[4] optimal truncation (plot_optimal_truncation)")
Φ = FormalSeries(:euler, 40)
save_fig("04a_truncation_single.png", plot_optimal_truncation(Φ, 0.1))
save_fig("04b_truncation_family.png", plot_optimal_truncation(Φ, [0.05, 0.1, 0.2, 0.4]))

println("\nDone — $(length(readdir(OUT))) figures in $OUT")
