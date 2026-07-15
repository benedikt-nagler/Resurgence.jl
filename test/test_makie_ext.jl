@testset "makie extension" begin
    # before the extension loads, the stubs explain what to do
    @test_throws ErrorException plot_borel_plane(borel(FormalSeries(:euler, 6)))
    @test_throws ErrorException plot_large_order(FormalSeries(:euler, 6))
    @test_throws ErrorException plot_optimal_truncation(FormalSeries(:euler, 6), 0.1)
    @eval using CairoMakie
    B = borel(FormalSeries(:euler, 8))
    fig = plot_borel_plane(B; order = 1, rays = [0.0, Float64(π)])
    @test fig isa CairoMakie.Makie.Figure
    r = pade(B; order = 1)
    @test plot_borel_plane(r) isa CairoMakie.Makie.Figure

    # residue-scaled markers (default) and the constant-size fallback both render
    @test plot_borel_plane(r; scale_by_residue = true) isa CairoMakie.Makie.Figure
    @test plot_borel_plane(r; scale_by_residue = false) isa CairoMakie.Makie.Figure

    # overlay: a vector of approximants in one axis (Padé vs Padé here)
    r2 = pade(B; order = 2, reduce = true)
    @test plot_borel_plane([r, r2]) isa CairoMakie.Makie.Figure
    @test plot_borel_plane([r, r2]; labels = ["low", "high"]) isa
          CairoMakie.Makie.Figure
    @test_throws ArgumentError plot_borel_plane(PadeApproximant[])

    # Padé-vs-AAA overlay (AAA lives in the BaryRational extension)
    @eval using BaryRational
    a = aaa_borel(B)
    @test plot_borel_plane(a) isa CairoMakie.Makie.Figure
    @test plot_borel_plane(AbstractBorelApproximant[r, a];
                           labels = ["Padé", "AAA"]) isa CairoMakie.Makie.Figure

    # large-order ratio test + optimal-truncation U-curve (#9)
    Φ = FormalSeries(:euler, 20)
    @test plot_large_order(Φ) isa CairoMakie.Makie.Figure
    @test plot_large_order(Φ; asymptote = false) isa CairoMakie.Makie.Figure
    @test plot_optimal_truncation(Φ, 0.1) isa CairoMakie.Makie.Figure
    @test plot_optimal_truncation(Φ, [0.05, 0.1, 0.2]) isa CairoMakie.Makie.Figure
    @test_throws ArgumentError plot_optimal_truncation(Φ, Float64[])
end
