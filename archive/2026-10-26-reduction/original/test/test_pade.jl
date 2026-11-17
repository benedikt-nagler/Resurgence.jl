@testset "pade" begin
    @testset "Euler oracle: [1/1] reproduces -1/(1+ζ) exactly" begin
        B = borel(FormalSeries(:euler, 8))
        r = pade(B.series.coeffs, 1, 1)
        @test r isa PadeApproximant{Rational{BigInt}}
        @test numerator_degree(r) == 1 && denominator_degree(r) == 1
        for ζ in (0//1, 1//2, 2//1, -1//2)
            @test r(ζ) == -1 // (1 + ζ)
        end
        # pole and residue at ζ = −1
        ζs = poles(r)
        @test length(ζs) == 1
        @test ζs[1] ≈ -1 atol = 1e-70
        ζ, res = only(residues(r))
        @test res ≈ -1 atol = 1e-70
    end

    @testset "exp oracle: [2/2] closed form, exact" begin
        Φ = FormalSeries(:exp, 6)
        r = pade(Φ, 2, 2)
        # known: exp ≈ (1 + x/2 + x²/12)/(1 − x/2 + x²/12)
        @test r.p == [1, 1//2, 1//12]
        @test r.q == [1, -1//2, 1//12]
    end

    @testset "degenerate Toeplitz block throws" begin
        # 1/(1−x) admits a [0/1] representation: the [1/2] system is singular
        c = fill(1//1, 5)
        @test_throws DegeneratePade pade(c, 1, 2)
        @test pade(c, 0, 1)(1//3) == 3//2            # while [0/1] is fine
    end

    @testset "guards" begin
        @test_throws InvalidArgument pade([1//1, 1//1], 1, 1)   # too few terms
        @test_throws InvalidArgument pade([1//1], -1, 0)
        @test pade([2//1, 3//1, 4//1], 2, 0).q == [1//1]         # M = 0: polynomial
    end

    @testset "numeric coefficients (RowMaximum path)" begin
        c = Float64[(-1.0)^(n + 1) for n in 0:6]
        r = pade(c, 1, 1)
        @test r(0.5) ≈ -1 / 1.5
        @test poles(r)[1] ≈ -1 atol = 1e-12
    end

    @testset "exact field elements (RowNonZero path)" begin
        F = AbstractAlgebra.GF(101)
        c = [F(1), F(1), F(1), F(1), F(1)]
        r = pade(c, 0, 1)
        @test r(F(3)) == F(1) / (F(1) - F(3))
    end

    @testset "degeneracy reduction (reduce = true)" begin
        # binding oracle: the [n/n] Padé of the exactly rational Euler Borel
        # transform −1/(1+ζ) reduces to the exact [1/1] for every n
        for n in 2:8
            r = pade(borel(FormalSeries(:euler, 2n + 1)); order = n, reduce = true)
            @test r.p == [-1//1, 0//1]
            @test r.q == [1//1, 1//1]
        end
        # the default still throws with the originally requested orders
        err = try
            pade(borel(FormalSeries(:euler, 9)); order = 4)
        catch e
            e
        end
        @test err isa DegeneratePade && err.L == 4 && err.M == 4
        # reduction of a non-degenerate request is a no-op
        c = Float64[(-1.0)^(k + 1) for k in 0:6]
        @test pade(c, 1, 1; reduce = true).q == pade(c, 1, 1).q
        # geometric series 1/(1−x): [1/2] reduces down to [0/1]
        r01 = pade(fill(1//1, 5), 1, 2; reduce = true)
        @test numerator_degree(r01) == 0 && denominator_degree(r01) == 1
        @test r01(1//3) == 3//2
        # a series that is exactly polynomial: the reduced denominator is trivial
        rpoly = pade([1//1, 2//1, 0//1, 0//1, 0//1], 2, 2; reduce = true)
        @test rpoly.q == [1//1, 0//1]
        @test rpoly(5//1) == 11
        # borel_pade_poles passes reduce through
        ζs = borel_pade_poles(FormalSeries(:euler, 11); reduce = true)
        @test only(ζs) ≈ -1 atol = 1e-70
    end

    @testset "pade(B; order) and borel_pade_poles" begin
        # near-diagonal default on a non-rational Borel function; the airy
        # constant term is split off, leaving 8 coefficients → [4/3]
        r = pade(borel(FormalSeries(:airy, 9)))
        @test numerator_degree(r) == 4 && denominator_degree(r) == 3
        Φ = FormalSeries(:euler, 9)
        B = borel(Φ)
        # the [4/4] Padé of the exactly rational -1/(1+ζ) is degenerate by design
        @test_throws DegeneratePade pade(B)
        r2 = pade(B; order = 1)
        @test numerator_degree(r2) == 1 && denominator_degree(r2) == 1
        ζs = borel_pade_poles(Φ; order = 1)
        @test ζs[1] ≈ -1 atol = 1e-70
        # Airy: leading Borel singularity at ζ = −4/3 (a branch point, so the
        # nearest Padé pole converges slowly)
        ζa = borel_pade_poles(FormalSeries(:airy, 40))
        @test ζa[1] ≈ -4 // 3 atol = 1e-2
    end
end
