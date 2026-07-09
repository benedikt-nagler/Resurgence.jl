@testset "makie extension" begin
    # before the extension loads, the stub explains what to do
    @test_throws ErrorException plot_borel_plane(borel(FormalSeries(:euler, 6)))
    @eval using CairoMakie
    B = borel(FormalSeries(:euler, 8))
    fig = plot_borel_plane(B; order = 1, rays = [0.0, Float64(π)])
    @test fig isa CairoMakie.Makie.Figure
    r = pade(B; order = 1)
    @test plot_borel_plane(r) isa CairoMakie.Makie.Figure
end
