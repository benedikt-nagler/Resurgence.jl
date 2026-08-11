@testset "makie extension" begin
    # the pre-load stubs are asserted in test_extension_stubs.jl, which runs first
    @eval using CairoMakie
    MK = CairoMakie.Makie

    # Every plot_* is a recipe, so the non-mutating call returns FigureAxisPlot.
    B = borel(FormalSeries(:euler, 8))
    fap = plot_borel_plane(B; order = 1, rays = [0.0, Float64(π)])
    @test fap isa MK.FigureAxisPlot
    r = pade(B; order = 1)
    @test plot_borel_plane(r) isa MK.FigureAxisPlot

    # residue-scaled markers (default) and the constant-size fallback both render
    @test plot_borel_plane(r; scale_by_residue = true) isa MK.FigureAxisPlot
    @test plot_borel_plane(r; scale_by_residue = false) isa MK.FigureAxisPlot

    # overlay: a vector of approximants in one axis (Padé vs Padé here)
    r2 = pade(B; order = 2, reduce = true)
    @test plot_borel_plane([r, r2]) isa MK.FigureAxisPlot
    @test plot_borel_plane([r, r2]; labels = ["low", "high"]) isa MK.FigureAxisPlot
    @test_throws ArgumentError plot_borel_plane(PadeApproximant[])

    # Padé-vs-AAA overlay (AAA lives in the BaryRational extension)
    @eval using BaryRational
    a = aaa_borel(B)
    @test plot_borel_plane(a) isa MK.FigureAxisPlot
    @test plot_borel_plane(AbstractBorelApproximant[r, a];
                           labels = ["Padé", "AAA"]) isa MK.FigureAxisPlot

    # large-order ratio test + optimal-truncation U-curve (#9)
    Φ = FormalSeries(:euler, 20)
    @test plot_large_order(Φ) isa MK.FigureAxisPlot
    @test plot_large_order(Φ; asymptote = false) isa MK.FigureAxisPlot
    @test plot_optimal_truncation(Φ, 0.1) isa MK.FigureAxisPlot
    @test plot_optimal_truncation(Φ, [0.05, 0.1, 0.2]) isa MK.FigureAxisPlot
    @test_throws ArgumentError plot_optimal_truncation(Φ, Float64[])

    # ── recipe forms: composition into a caller's own layout ──────────────────
    @testset "recipe composition" begin
        # a layout slot as the target returns AxisPlot, not a new Figure
        fig = MK.Figure()
        @test plot_borel_plane(fig[1, 1], r) isa MK.AxisPlot
        @test plot_large_order(fig[1, 2], Φ) isa MK.AxisPlot
        @test plot_optimal_truncation(fig[2, 1], Φ, 0.1) isa MK.AxisPlot

        # the three panels really landed in one figure
        @test length(fig.content) >= 3

        # mutating forms draw into an existing axis and return the plot object
        ax = MK.Axis(fig[2, 2])
        @test plot_borel_plane!(ax, r) isa MK.AbstractPlot
        @test plot_borel_plane!(ax, [r, r2]) isa MK.AbstractPlot
        @test plot_large_order!(ax, Φ) isa MK.AbstractPlot
        @test plot_optimal_truncation!(ax, Φ, [0.05, 0.1]) isa MK.AbstractPlot

        # the recipe supplies axis decoration, and the caller can override it
        fig2 = MK.Figure()
        ap = plot_borel_plane(fig2[1, 1], r)
        @test ap.axis.xlabel[] == "Re ζ"
        ap2 = plot_borel_plane(fig2[1, 2], r; axis = (; xlabel = "mine"))
        @test ap2.axis.xlabel[] == "mine"
    end
end
