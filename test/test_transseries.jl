@testset "transseries" begin
    Φ0 = FormalSeries([0//1, 1, 2, 6])                 # some 4-term series
    Φ1 = FormalSeries([1//1, -1, 1, -1])
    F = Transseries(1//1, [Φ0, Φ1])

    @testset "construction and accessors" begin
        @test F isa Transseries{Rational{Int},Rational{Int}}
        @test action(F) == 1
        @test n_sectors(F) == 2
        @test variable(F) == :ħ
        @test sector(F, 0) == Φ0 && sector(F, 1) == Φ1
        @test all(iszero, coefficients(sector(F, 5)))   # beyond storage: zero series
        @test_throws InvalidArgument sector(F, -1)
        @test sectors(F) == [Φ0, Φ1]
        @test is_exact(F)
        @test eltype(typeof(F)) == Rational{Int}

        @test_throws InvalidArgument Transseries(0//1, [Φ0])          # zero action
        @test_throws InvalidArgument Transseries(1//1, FormalSeries{Rational{Int}}[])
        Ψ = FormalSeries([1//1], :x)
        @test_throws IncompatibleSeries Transseries(1//1, [Φ0, Ψ])    # var mismatch
    end

    @testset "arithmetic" begin
        G = Transseries(1//1, [Φ1])
        S = F + G
        @test n_sectors(S) == 2
        @test sector(S, 0) == Φ0 + Φ1                  # overlapping sector adds
        @test sector(S, 1) == Φ1                       # missing sector passes through
        @test S - G == F
        @test -(-F) == F
        @test sector(3 * F, 1) == 3 * Φ1 && F * 3 == 3 * F

        # weight-graded Cauchy product: (F·G)_n = Σ_i F_i G_{n-i}
        P = F * F
        @test n_sectors(P) == 2
        @test sector(P, 0) == Φ0 * Φ0
        @test sector(P, 1) == Φ0 * Φ1 + Φ1 * Φ0
        # one-sector G truncates the product to one sector
        @test n_sectors(F * G) == 1

        H = Transseries(2//1, [Φ0])
        @test_throws IncompatibleSeries F + H          # action mismatch
        @test_throws IncompatibleSeries F * H
        K = Transseries(1//1, [FormalSeries([1//1], :x)])
        @test_throws IncompatibleSeries F + K          # var mismatch
    end

    @testset "truncate" begin
        Ft = truncate(F; sectors = 1)
        @test n_sectors(Ft) == 1 && sector(Ft, 0) == Φ0
        Ftt = truncate(F; terms = 2)
        @test n_sectors(Ftt) == 2 && n_terms(sector(Ftt, 0)) == 2
        @test_throws InvalidArgument truncate(F; sectors = 3)
    end

    @testset "transseries_sum: Euler transseries" begin
        E = Transseries(:euler, 12)
        setprecision(BigFloat, 256) do
            ħ = big"0.2"
            B = borel(sector(E, 0))
            s0 = lateral_sum(B, ħ; side = :plus, order = 1)
            # σ = 0 reduces to the lateral sum of Φ₀
            @test transseries_sum(E, big"0.0", ħ; order = 1) ≈ s0 rtol = 1e-40
            # the constant one-instanton sector enters exactly as σ e^{-A/ħ}
            σ = big"0.7"
            s = transseries_sum(E, σ, ħ; order = 1)
            @test s ≈ s0 + σ * exp(-1 / ħ) rtol = 1e-40
            # side is forwarded to every sector
            sm = transseries_sum(E, σ, ħ; side = :minus, order = 1)
            @test sm ≈ lateral_sum(B, ħ; side = :minus, order = 1) + σ * exp(-1 / ħ) rtol = 1e-40
        end
    end
end
