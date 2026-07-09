@testset "hyperasymptotics" begin
    E₁(x) = first(quadgk(t -> exp(-t) / t, x, typeof(x)(Inf);
                         rtol = eps(typeof(x))^(4 // 5)))

    @testset "terminant basics" begin
        @test_throws InvalidArgument dingle_terminant(-1.0, 0.5)
        @test_throws InvalidArgument dingle_terminant(2.0, 0.5; side = :sideways)
        # Λ_ν(0) = Γ(ν)
        @test dingle_terminant(big"2.5", big"0.0") ≈
              Resurgence._gamma(big"2.5") rtol = 1e-30
        @test dingle_terminant(3.0, 0.0) ≈ 2.0 rtol = 1e-10
        @test dingle_terminant(1.5, 0.0) ≈ SpecialFunctions.gamma(1.5) rtol = 1e-10
    end

    @testset "terminant Stokes jump: Λ⁺ − Λ⁻ = 2πi x^{-ν} e^{-1/x}" begin
        setprecision(BigFloat, 256) do
            for (ν, x) in ((3, big"0.4"), (big"2.5", big"0.25"))
                jump = dingle_terminant(ν, x; side = :plus) -
                       dingle_terminant(ν, x; side = :minus)
                @test jump ≈ 2im * big(π) * x^(-big(ν)) * exp(-1 / x) rtol = 1e-25
            end
        end
    end

    @testset "exact Euler remainder pins the normalization" begin
        # φ_exact − Σ_{n<N} aₙħ^{n+1} == (ħ/A)^{N+1} Λ_{N+1}(ħ/A), A = −1: for the
        # Euler series the large-order relation is exact, so this is an identity
        setprecision(BigFloat, 256) do
            Φ = FormalSeries(:euler, 40)
            for ħ in (big"0.1", big"0.5"), N in (5, 10, 20)
                exact = -exp(1 / ħ) * E₁(1 / ħ)
                remainder = exact - partial_sums(Φ, ħ)[N]
                @test remainder ≈ (-ħ)^(N + 1) *
                                  dingle_terminant(N + 1, -ħ) rtol = 1e-30
            end
        end
    end

    @testset "optimal_truncation" begin
        setprecision(BigFloat, 256) do
            Φ = FormalSeries(:euler, 40)
            ħ = big"0.1"
            exact = -exp(1 / ħ) * E₁(1 / ħ)
            for rule in (:smallest_term, :action)
                t = optimal_truncation(Φ, ħ; rule)
                @test 8 ≤ t.N ≤ 11                       # N* ≈ |A/ħ| = 10
                @test abs(t.value - exact) ≤ 2 * t.error  # error estimate honest
                @test t.error < big"1e-4"                 # ~ e^{-1/ħ} scale
            end
            # explicit action bypasses the fit
            t2 = optimal_truncation(Φ, ħ; rule = :action, A = -1)
            @test t2.N == 9                              # round(|A/ħ|) − offset
            @test_throws InvalidArgument optimal_truncation(Φ, ħ; rule = :nope)
        end
    end

    @testset "hyper_sum: exact on Euler, N-independent" begin
        # adjacent sector is the constant 1, so level 1 is exact: hyper_sum
        # reproduces φ_exact at machine-quadrature precision for any truncation
        setprecision(BigFloat, 256) do
            Φ = FormalSeries(:euler, 40)
            adj = FormalSeries([1 // 1])
            for ħ in (big"0.2", big"0.5")
                exact = -exp(1 / ħ) * E₁(1 / ħ)
                for N in (5, 10, 15)
                    h = hyper_sum(Φ, ħ; adjacent = adj, action = -1,
                                  stokes = -2 * big(π) * im, N)
                    @test h ≈ exact rtol = 1e-25
                end
                # default N (action rule) too
                h = hyper_sum(Φ, ħ; adjacent = adj, action = -1,
                              stokes = -2 * big(π) * im)
                @test h ≈ exact rtol = 1e-25
            end
            # binding inequality: level 1 beats superasymptotics by ≥ 5 orders
            ħ = big"0.2"
            exact = -exp(1 / ħ) * E₁(1 / ħ)
            t = optimal_truncation(Φ, ħ; rule = :action, A = -1)
            h = hyper_sum(Φ, ħ; adjacent = adj, action = -1,
                          stokes = -2 * big(π) * im)
            @test abs(h - exact) < 1e-5 * t.error
        end
    end

    @testset "hyper_sum on a Stokes ray: lateral terminants" begin
        # mirror Euler Σ n! ħ^{n+1}: Borel pole at +1 on the summation ray; the
        # side kwarg must reproduce the corresponding lateral Borel sum
        setprecision(BigFloat, 256) do
            Φ = FormalSeries([Rational{BigInt}(factorial(big(n))) for n in 0:39];
                             power_offset = 1)
            B = borel(Φ)
            adj = FormalSeries([1 // 1])
            ħ = big"0.15"
            for side in (:plus, :minus)
                lat = lateral_sum(B, ħ; side, order = 1)
                h = hyper_sum(Φ, ħ; adjacent = adj, action = 1,
                              stokes = -2 * big(π) * im, side)
                @test h ≈ lat rtol = 1e-25
            end
        end
    end
end
