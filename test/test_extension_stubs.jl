# Extension stubs must be checked before anything loads Makie or BaryRational.
# Once an extension is loaded it stays loaded for the rest of the session, so these
# assertions live in their own file, included first, rather than in the per-extension
# test files (where an earlier `@eval using ...` would have already armed the real
# method and the stub would never throw).
@testset "extension stubs" begin
    Φ = FormalSeries(:euler, 9)

    # Makie extension
    @test_throws ErrorException plot_borel_plane(borel(Φ))
    @test_throws ErrorException plot_large_order(Φ)
    @test_throws ErrorException plot_optimal_truncation(Φ, 0.1)

    # BaryRational extension
    @test_throws ErrorException aaa_borel(borel(Φ))
end
