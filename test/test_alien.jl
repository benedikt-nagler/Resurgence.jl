# Alien calculus oracles. All sign conventions cohere with the M1 pin
# disc = s₊ − s₋ = +2πi e^{-1/ħ} for Σ n! ħ^{n+1}:  S₁ = 2πi·Res = −2πi,
# s₋(F)(σ) = s₊(F)(σ + S₁),  disc(Φ₀) = −S₁ e^{-A/ħ} s(Φ₁),
# s_med(F)(σ) = s₊(F)(σ + S₁/2).
@testset "alien" begin
    @testset "Euler Stokes constant: S₁ = −2πi from the exact residue" begin
        E = Transseries(:euler, 12)
        setprecision(BigFloat, 256) do
            S1 = stokes_constant(borel(sector(E, 0)); order = 1)
            @test S1 ≈ -2im * big(π) rtol = 1e-40
            # the alternating orientation has its pole at ζ = −1 — same constant
            S1m = stokes_constant(borel(FormalSeries(:euler, 12));
                                  order = 1, location = -1)
            @test S1m ≈ -2im * big(π) rtol = 1e-40
            # coherence with the M1 discontinuity pin: disc = −S₁ e^{-A/ħ} s(Φ₁)
            ħ = big"0.2"
            disc = stokes_discontinuity(borel(sector(E, 0)), ħ; order = 1)
            @test disc ≈ -S1 * exp(-1 / ħ) rtol = 1e-30
        end
        # large-order route on the same series: |S₁| = 2π S_fit A^b = 2π
        @test stokes_constant(sector(Transseries(:euler, 30), 0)) ≈ 2 * big(π) rtol = 1e-25
        # no Borel poles → no Stokes constant to read
        @test_throws InvalidArgument stokes_constant(borel(FormalSeries([1//1];
                                                                        power_offset = 1//1)))
    end

    @testset "bridge equation: Δ_{lA} Φ_m = S_l (m+l) Φ_{m+l}" begin
        E = Transseries(:euler, 12)
        S1 = -2im * big(π)
        D = alien_derivative(E; stokes = S1)
        @test action(D) == action(E) && variable(D) == variable(E)
        @test n_sectors(D) == 1
        @test sector(D, 0) == S1 * sector(E, 1)            # Δ_A Φ₀ = S₁ Φ₁
        # Δ_A annihilates the one-instanton tail: Δ_A² F = 0
        D2 = alien_derivative(D; stokes = S1)
        @test n_sectors(D2) == 1 && all(iszero, coefficients(sector(D2, 0)))
        # l = 2 needs S₂ and sees no sector here
        D_2A = alien_derivative(E, 2; stokes = [S1, 0im])
        @test all(iszero, coefficients(sector(D_2A, 0)))
        # general sector map on a 3-sector transseries
        Φ = [FormalSeries([1//1, k]) for k in 1:3]
        G = alien_derivative(Transseries(1//1, Φ), 1; stokes = 2)
        @test n_sectors(G) == 2
        @test sector(G, 0) == 2 * 1 * Φ[2] && sector(G, 1) == 2 * 2 * Φ[3]
        # validation
        @test_throws InvalidArgument alien_derivative(E, 0; stokes = S1)
        @test_throws InvalidArgument alien_derivative(E, 2; stokes = S1)   # scalar, l = 2
        @test_throws InvalidArgument alien_derivative(E, 2; stokes = [S1]) # too short
    end

    @testset "connection formula: s₋(F)(σ) = s₊(F)(σ + S₁)  [Euler]" begin
        E = Transseries(:euler, 12)
        setprecision(BigFloat, 256) do
            S1 = stokes_constant(borel(sector(E, 0)); order = 1)
            for (σ, ħ) in ((big"0.0", big"0.2"), (big"0.7", big"0.1"),
                           (big"-1.3", big"0.3"))
                lhs = transseries_sum(E, σ, ħ; side = :minus, order = 1)
                σ′ = stokes_automorphism(E, σ; stokes = S1)
                @test σ′ == σ + S1
                rhs = transseries_sum(E, σ′, ħ; side = :plus, order = 1)
                @test lhs ≈ rhs rtol = 1e-30
            end
        end
    end

    @testset "median summation is real and matches the lateral average [Euler]" begin
        E = Transseries(:euler, 12)
        setprecision(BigFloat, 256) do
            S1 = -2im * big(π)
            for ħ in (big"0.1", big"0.3")
                med = median_sum(E, big"0.0", ħ; stokes = S1, order = 1)
                @test abs(imag(med)) < big"1e-40"
                # leading-sector method: the average of the two lateral sums
                med0 = median_sum(sector(E, 0), ħ; order = 1)
                @test med ≈ med0 rtol = 1e-40
                # the two representatives agree: s₋ at σ − S₁/2 gives the same value
                alt = transseries_sum(E, -S1 / 2, ħ; side = :minus, order = 1)
                @test med ≈ alt rtol = 1e-30
            end
        end
    end

    @testset "Airy: |S₁| = 1 and the large-order ↔ alien-derivative relation" begin
        setprecision(BigFloat, 512) do
            Φa = FormalSeries(:airy, 60)
            fit = large_order_fit(Φa; order = 6)
            @test fit.A ≈ 4 / big"3.0" rtol = 1e-8
            @test abs(fit.b) < 1e-6                          # b = 0 for Airy
            @test stokes_constant(Φa; order = 6) ≈ 1 rtol = 1e-6   # |S₁| = 1

            # resurgence large-order relation with b = 0, A_c = −4/3:
            #   |a_n| (4/3)^n / Γ(n) → |S₁|/2π,  and the 1/n correction carries
            #   the first one-instanton coefficient c₁ of Δ_A Φ₀ = S₁ Φ₁:
            #   2π w_n = 1 + c₁·(−4/3)/(n−1) + O(1/n²)
            a = coefficients(Φa)
            w = [abs(a[n + 1]) * (big"4.0" / 3)^n / factorial(big(n - 1))
                 for n in 1:59]
            @test 2 * big(π) * richardson(w, 6) ≈ 1 rtol = 1e-6
            v = [(2 * big(π) * w[n] - 1) * (n - 1) for n in 2:59]
            # the API prediction of c₁: sector 0 of Δ_A(Airy transseries) ∝ airy_bi
            D = alien_derivative(Transseries(:airy, 6); stokes = im)
            c = coefficients(sector(D, 0))
            c1 = real(c[2] / c[1])                           # = 5/48
            @test c1 == 5//48
            @test richardson(v, 4) ≈ -c1 * (4 / big"3.0") rtol = 1e-3
        end
    end

    @testset "quartic instanton: A = 1/3, b = 1/2, S = √6/π^{3/2}" begin
        setprecision(BigFloat, 512) do
            # the 1/n corrections are large (Zinn-Justin: −95/72n), so this needs
            # more terms than Euler/Airy — 80 exact terms cost ~0.3 s
            Φq = FormalSeries(:quartic, 80)
            fit = large_order_fit(Φq; order = 8)
            @test fit.A ≈ 1 / big"3.0" rtol = 1e-6
            @test fit.b ≈ 1 / big"2.0" atol = 1e-4
            @test fit.S ≈ sqrt(big"6.0") / big(π)^(big"3.0" / 2) rtol = 1e-4
            # |S₁| = 2π S A^b = 2√(2/π)
            @test stokes_constant(Φq; order = 8) ≈ 2 * sqrt(2 / big(π)) rtol = 1e-4
        end
    end
end
