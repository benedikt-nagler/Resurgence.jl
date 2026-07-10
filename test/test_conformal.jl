@testset "conformal" begin
    @testset "map identities" begin
        m = conformal_map(big(-1 // 1))
        @test m(0) == 0
        # round trip w ∘ ζ on a grid off the cut
        for w in (big"0.3" + 0im, big"-0.4" + big"0.2" * im, big"0.7" * im)
            @test m(inverse(m, w)) ≈ w atol = 1e-70
        end
        for ζ in (big"0.5" + 0im, big"2.0" + big"1.0" * im, big"-0.9" + 0im)
            @test inverse(m, m(ζ)) ≈ ζ atol = 1e-68
        end
        # the branch point maps to w = 1; points on the cut to |w| = 1
        @test m(big(-1)) ≈ 1 atol = 1e-70
        @test abs(m(big"-2.5")) ≈ 1 atol = 1e-70
        # symmetric pair map
        mp = conformal_map(big(2im); pair = true)
        @test mp(0) == 0
        for ζ in (big"0.5" + 0im, big"1.0" + big"0.3" * im)
            @test inverse(mp, mp(ζ)) ≈ ζ atol = 1e-68
        end
        @test abs(mp(big(2im))) ≈ 1 atol = 1e-70
        @test abs(mp(big(-2im))) ≈ 1 atol = 1e-70
        @test_throws InvalidArgument conformal_map(0)
    end

    @testset "Euler: exact rational composition, [2/2] in w" begin
        B = borel(FormalSeries(:euler, 13))
        c = conformal_borel(B; zeta0 = -1 // 1)
        # −1/(1+ζ(w)) = −(1+w)²/(1−w)²: rational in w, caught exactly by the
        # reduced inner Padé
        @test c.inner isa PadeApproximant{Rational{BigInt}}
        @test numerator_degree(c.inner) == 2 && denominator_degree(c.inner) == 2
        # inner is exactly −(1+w)²/(1−w)² (exact rational composition)
        for w in (1 // 3, -1 // 5)
            @test c.inner(w) == -(1 + w)^2 / (1 - w)^2
        end
        # end-to-end c(ζ) = −1/(1+ζ): the map's sqrt needs a BigFloat argument to
        # keep full precision (rational ζ would route through a Float64 sqrt)
        setprecision(BigFloat, 256) do
            for ζ in (big"0.4", big"-0.3", big"1.5")
                @test c(ζ) ≈ -1 / (1 + ζ) atol = 1e-60
            end
        end
        # auto-detected ζ₀ agrees
        c2 = conformal_borel(B)
        @test c2.map.zeta0 ≈ -1 atol = 1e-15
        # poles: the branch point −1 is reported (and guards the cut direction)
        @test minimum(abs.(poles(c) .+ 1)) < 1e-15
        # integrating along the cut (θ = π, ħ < 0 so e^{-ζ/ħ} decays) hits the
        # branch point → PoleOnRay (Float64 ħ keeps the on-ray tolerance loose
        # enough to absorb the Float64-π angle)
        @test_throws PoleOnRay laplace_sum(B, c, -0.5; θ = Float64(π))
    end

    @testset "Euler summation through the seam is quadrature-exact" begin
        # (the Euler Borel function is rational, so plain Padé is already exact —
        # the conformal approximant must reproduce the same value; the genuine
        # digits-vs-order gain is tested on Airy below)
        setprecision(BigFloat, 256) do
            ħ = big"0.7"
            B = borel(FormalSeries(:euler, 13))
            exact = laplace_sum(B, ħ; order = 1)
            conf = laplace_sum(B, conformal_borel(B; zeta0 = -1 // 1), ħ)
            @test conf ≈ exact rtol = 1e-40
        end
    end

    @testset "digits vs order: conformal beats plain Borel–Padé (Airy)" begin
        setprecision(BigFloat, 256) do
            ħ = big"1.0"
            # high-order reference (validated to ~1e-25 against Ai(1) in
            # test_laplace.jl)
            Φref = FormalSeries(BigFloat.(coefficients(FormalSeries(:airy, 220))))
            ref = lateral_sum(borel(Φref), ħ; side = :plus, tilt = 3 // 10,
                              rtol = big"1e-40")
            prev = nothing
            for n in (12, 20, 30)
                Φ = FormalSeries(BigFloat.(coefficients(FormalSeries(:airy, n))))
                B = borel(Φ)
                plain = lateral_sum(B, pade(B), ħ; side = :plus, tilt = 3 // 10)
                conf = lateral_sum(B, conformal_borel(B; zeta0 = -4 // 3), ħ;
                                   side = :plus, tilt = 3 // 10)
                err_plain = abs(plain - ref)
                err_conf = abs(conf - ref)
                @test err_conf < err_plain
                prev !== nothing && @test err_conf < prev
                prev = err_conf
            end
        end
    end

    @testset "NoSingularityFound" begin
        # a polynomial Borel transform has no Padé poles anywhere
        B = borel(FormalSeries([0 // 1, 1, 0, 0, 0, 0, 0]; power_offset = 1))
        @test_throws NoSingularityFound conformal_borel(B)
    end
end
