# Constructed two-pole oracle: a_n = n!((-1)^n + s(-1/2)^n), s = 1/3, power
# offset 1 ⇒ reduced Borel function exactly 1/(1+ζ) + s/(1+ζ/2), poles -1 and -2,
# multi-saddle triples (A, b, S) = (-1, 0, -1) and (-2, 0, -2s).
@testset "singularities" begin
    s = 1 // 3
    two_pole = FormalSeries([Rational{BigInt}(factorial(big(n))) *
                             ((-1)^n + s * (-1 // 2)^n) for n in 0:24];
                            power_offset = 1)
    Bp = borel(two_pole)

    @testset "guards" begin
        @test_throws InvalidArgument subtract_singularity(Bp, 0)
        @test_throws InvalidArgument subtract_singularity(Bp, -1; beta = 0)
        @test_throws InvalidArgument large_order_fit(two_pole; saddles = 0)
        # a polynomial Borel transform has no singularity to find
        Bpoly = borel(FormalSeries([0 // 1, 1, 0, 0, 0, 0, 0]; power_offset = 1))
        @test_throws NoSingularityFound subtract_singularity(Bpoly)
        # more saddles than Borel–Padé poles
        @test_throws NoSingularityFound large_order_fit(FormalSeries(:euler, 21);
                                                        saddles = 3)
    end

    @testset "Darboux peeling: simple pole, residue route" begin
        out = subtract_singularity(Bp, -1 // 1)
        @test out.c ≈ 1 atol = 1e-25
        @test out.beta == 1
        # auto-detected ζ₀ agrees
        auto = subtract_singularity(Bp)
        @test auto.zeta0 ≈ -1 atol = 1e-25
        @test auto.c ≈ 1 atol = 1e-25
        # the numerically peeled series exposes the subleading pole at −2
        ζs = poles(pade(out.peeled; reduce = true))
        @test minimum(abs.(ζs .+ 2)) < 1e-20
    end

    @testset "Darboux peeling: exact amplitude peels exactly" begin
        out = subtract_singularity(Bp, -1 // 1; c = 1 // 1)
        # peeled reduced coefficients are exactly s(−1/2)ⁿ
        @test out.peeled.series.coeffs == [s * (-1 // 2)^n for n in 0:24]
        @test out.peeled.beta == Bp.beta
        @test only(poles(pade(out.peeled; reduce = true))) ≈ -2 atol = 1e-70
        # the subleading action is now visible to the large-order fit
        fit = large_order_fit(inverse_borel(out.peeled))
        @test fit.A ≈ 2 atol = 1e-20
    end

    @testset "Darboux peeling: branch point via the acceleration route" begin
        # b(ζ) = 2(1 − ζ/(−1))^{-1/2} + regular part Σ(−ζ/5)ⁿ, exact rationals
        poch = Rational{BigInt}[1]
        for k in 0:23
            push!(poch, poch[end] * (1 // 2 + k))
        end
        b = [2 * poch[n + 1] // factorial(big(n)) * (-1)^n + (-1 // 5)^n
             for n in 0:24]
        Φh = FormalSeries([b[n + 1] * factorial(big(n)) for n in 0:24];
                          power_offset = 1)
        out = subtract_singularity(borel(Φh), -1; beta = 1 // 2)
        @test out.c ≈ 2 atol = 1e-10
    end

    @testset "multi-saddle fit: two real actions, exact model" begin
        fit = large_order_fit(two_pole; saddles = 2)
        @test fit.converged
        @test length(fit.saddles) == 2
        s1, s2 = fit.saddles           # sorted by |A|
        @test s1.A ≈ -1 atol = 1e-15
        @test s1.b ≈ 0 atol = 1e-15
        @test s1.S ≈ -1 atol = 1e-15
        @test s2.A ≈ -2 atol = 1e-15
        @test s2.b ≈ 0 atol = 1e-15
        @test s2.S ≈ -2 * s atol = 1e-15
    end

    @testset "multi-saddle fit: oscillatory conjugate pair" begin
        # A = (3+4i)/5 on the unit circle, S = 1+2i:
        # a_n = n!·2Re[S A^{-(n+1)}] — the case ratio-based fits cannot see
        A = Complex{Rational{BigInt}}(3 // 5, 4 // 5)
        S = Complex{Rational{BigInt}}(1, 2)
        aosc = [factorial(big(n)) * 2 * real(S * conj(A)^(n + 1)) for n in 0:29]
        Φosc = FormalSeries(Rational{BigInt}.(aosc); power_offset = 1)
        fit = large_order_fit(Φosc; saddles = 2)
        @test fit.converged
        @test length(fit.saddles) == 2
        up = fit.saddles[findfirst(x -> imag(x.A) > 0, fit.saddles)]
        @test up.A ≈ (3 + 4im) // 5 atol = 1e-12
        @test up.b ≈ 0 atol = 1e-12
        @test up.S ≈ 1 + 2im atol = 1e-12
        # the two members are conjugate
        dn = fit.saddles[findfirst(x -> imag(x.A) < 0, fit.saddles)]
        @test dn.A ≈ conj(up.A) atol = 1e-30
        @test dn.S ≈ conj(up.S) atol = 1e-30
        # documented failure of the single-saddle ratio fit on oscillatory data
        single_ok = try
            fit1 = large_order_fit(Φosc)
            isfinite(fit1.A) && abs(abs(fit1.A) - 1) < 1e-2 && abs(fit1.b) < 1e-2
        catch e
            e isa ResurgenceError ? false : rethrow()
        end
        @test !single_ok
    end
end
