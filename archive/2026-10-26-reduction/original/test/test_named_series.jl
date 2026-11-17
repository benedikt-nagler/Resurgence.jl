@testset "named_series" begin
    @testset "euler" begin
        Φ = FormalSeries(:euler, 5)
        @test Φ isa FormalSeries{Rational{BigInt}}
        @test power_offset(Φ) == 1
        # Σ n!(−ħ)^{n+1}: coefficients (−1)^{n+1} n!
        @test coefficients(Φ) == [-1, 1, -2, 6, -24]
        @test is_exact(Φ)
    end

    @testset "airy" begin
        Φ = FormalSeries(:airy, 12)
        @test power_offset(Φ) == 0
        c = coefficients(Φ)
        @test c[1] == 1 && c[2] == -5//48
        # the underlying u_k = |c_k| (2/3)^k satisfy the classical recursion
        # u_{k+1} = u_k (6k+5)(6k+1)/(72(k+1))
        u = [abs(c[k + 1]) * (2//3)^k for k in 0:11]
        @test all(u[k + 2] == u[k + 1] * (6k + 5) * (6k + 1) // (72 * (k + 1))
                  for k in 0:9)
        @test all(sign(c[k]) == (-1)^(k + 1) for k in 1:12)   # alternating
    end

    @testset "airy_bi" begin
        Φ = FormalSeries(:airy_bi, 12)
        bi = coefficients(Φ)
        ai = coefficients(FormalSeries(:airy, 12))
        # the same u_k (3/2)^k without the alternating sign
        @test all(bi[k] == (-1)^(k + 1) * ai[k] for k in 1:12)
        @test bi[1] == 1 && bi[2] == 5//48
        @test all(c > 0 for c in bi)
    end

    @testset "quartic (Bender–Wu)" begin
        Φ = FormalSeries(:quartic, 6)
        E = coefficients(Φ)
        @test power_offset(Φ) == 0
        # E₀…E₃ of H = p²/2 + x²/2 + g x⁴ (classic Bender–Wu values)
        @test E[1] == 1//2 && E[2] == 3//4 && E[3] == -21//8 && E[4] == 333//16
        # alternating from E₁ on
        @test all(sign(E[k]) == (-1)^k for k in 2:6)
    end

    @testset "phi4 (0d partition function)" begin
        Z = FormalSeries(:phi4, 8; var = :g)
        a = coefficients(Z)
        @test power_offset(Z) == 0
        @test a[1] == 1 && a[2] == -3 && a[3] == 105//2 && a[4] == -3465//2
        # aₙ = (−1)ⁿ (4n−1)!!/n! - the Wick-contraction count - exactly
        dfact(m) = prod(big(m):-2:1)
        @test all(a[n + 1] == (-1)^n * dfact(4n - 1) // factorial(big(n)) for n in 1:7)
        # ratio −(4n+1)(4n+3)/(n+1) → −16n: Borel singularity at ζ = −1/16
        @test a[8] // a[7] == -(4 * 6 + 1) * (4 * 6 + 3) // 7
        # Borel sum recovers the integral it came from (normalization 1/√(2π))
        Z40 = FormalSeries(:phi4, 40; var = :g)
        exact, _ = quadgk(x -> exp(-x^2 / 2 - 0.1x^4) / sqrt(2π), -Inf, Inf)
        @test isapprox(real(borel_sum(Z40, 0.1; order = 15)), exact; rtol = 1e-10)
    end

    @testset "exp" begin
        Φ = FormalSeries(:exp, 5; var = :x)
        @test variable(Φ) == :x
        @test coefficients(Φ) == [1, 1, 1//2, 1//6, 1//24]
    end

    @testset "named transseries" begin
        E = Transseries(:euler, 5)
        @test action(E) == 1 && n_sectors(E) == 2
        # Φ₀ is the ħ → −ħ mirror of FormalSeries(:euler): coefficients +k!
        @test coefficients(sector(E, 0)) == [1, 1, 2, 6, 24]
        @test power_offset(sector(E, 0)) == 1
        @test sector(E, 1) == FormalSeries([1//1])
        A = Transseries(:airy, 8)
        @test action(A) == -4//3
        @test sector(A, 0) == FormalSeries(:airy, 8)
        @test sector(A, 1) == FormalSeries(:airy_bi, 8)
        @test variable(Transseries(:euler, 3; var = :x)) == :x
        @test_throws InvalidArgument Transseries(:nope, 3)
        @test_throws InvalidArgument Transseries(:euler, 0)
    end

    @testset "invalid" begin
        @test_throws InvalidArgument FormalSeries(:nope, 3)
        @test_throws InvalidArgument FormalSeries(:euler, 0)
    end
end
