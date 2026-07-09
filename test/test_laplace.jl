# The [n/n] Padé of the *exactly rational* Euler/geometric Borel functions is
# degenerate for n ≥ 2 (opt-in automatic reduction: `reduce = true`), so these
# tests pass `order = 1` explicitly.

# A minimal hand-rolled approximant exercising the AbstractBorelApproximant seam:
# the exact reduced Euler Borel function −1/(1+ζ) with its pole declared manually.
struct EulerDummyApprox <: AbstractBorelApproximant end
(::EulerDummyApprox)(ζ) = -1 / (1 + ζ)
Resurgence.poles(::EulerDummyApprox; refine::Bool = true) = [Complex{BigFloat}(-1)]

@testset "laplace" begin
    @testset "_gamma cross-check vs SpecialFunctions" begin
        for x in (0.5, 1.0, 1.5, 2.0, 3.25, 7.0)
            @test Resurgence._gamma(x) ≈ SpecialFunctions.gamma(x) rtol = 1e-10
        end
        @test Resurgence._gamma(big"0.75") ≈ SpecialFunctions.gamma(0.75) rtol = 1e-10
        @test_throws InvalidArgument Resurgence._gamma(-1.5)
    end

    @testset "Euler oracle: borel_sum matches -e^{1/ħ}E₁(1/ħ)" begin
        Φ = FormalSeries(:euler, 12)
        E₁(x) = first(quadgk(t -> exp(-t) / t, x, typeof(x)(Inf);
                             rtol = eps(typeof(x))^(4 // 5)))
        setprecision(BigFloat, 256) do
            for ħ in (big"0.1", big"0.25", big"0.5")
                closed = -exp(1 / ħ) * E₁(1 / ħ)
                s = borel_sum(Φ, ħ; order = 1)
                @test abs(imag(s)) < 1e-50
                @test real(s) ≈ closed rtol = 1e-45
            end
        end
        # Float64 path, cross-checked against SpecialFunctions
        s64 = borel_sum(Φ, 0.2; order = 1)
        @test real(s64) ≈ -exp(5.0) * SpecialFunctions.expint(5.0) rtol = 1e-10
    end

    @testset "convergent sanity: Σ ħ^{n+1} → ħ/(1-ħ)" begin
        Φ = FormalSeries(fill(1//1, 41); power_offset = 1//1)
        # reduced Borel function is e^ζ; an odd-degree Padé denominator of exp has
        # one real (spurious, far-out) zero that trips the ray guard → even order.
        # [n/n] error ~ (ħ/(1-ħ))^{2n+2}/binomial(2n,n), worst at ħ = 1/2.
        for ħ in (0.1, 0.3, 0.5)
            @test real(borel_sum(Φ, ħ; order = 18)) ≈ ħ / (1 - ħ) rtol = 1e-8
        end
    end

    @testset "PoleOnRay guard and lateral sums" begin
        # Σ n! ħ^{n+1}: reduced Borel function 1/(1-ζ), pole on the θ = 0 ray
        Φ = FormalSeries([Rational{BigInt}(factorial(big(n))) for n in 0:11];
                         power_offset = 1//1)
        B = borel(Φ)
        ħ = big"0.1"
        @test_throws PoleOnRay laplace_sum(B, ħ; order = 1)
        up = lateral_sum(B, ħ; side = :plus, order = 1)
        dn = lateral_sum(B, ħ; side = :minus, order = 1)
        @test up ≈ conj(dn) rtol = 1e-30            # real series: Schwarz reflection
        # lateral sums are tilt-independent (contour deformation, no pole crossed)
        up2 = lateral_sum(B, ħ; side = :plus, order = 1, tilt = 1//25)
        @test up ≈ up2 rtol = 1e-30
        @test_throws InvalidArgument lateral_sum(B, ħ; side = :sideways, order = 1)
        # e^{-ζ/ħ} must decay along the ray
        @test_throws InvalidArgument laplace_sum(B, -ħ; θ = 0, order = 1,
                                                 check_poles = false)
    end

    @testset "Stokes discontinuity oracle: +2πi e^{-1/ħ}" begin
        Φ = FormalSeries([Rational{BigInt}(factorial(big(n))) for n in 0:11];
                         power_offset = 1//1)
        B = borel(Φ)
        setprecision(BigFloat, 256) do
            for ħ in (big"0.1", big"0.2")
                disc = stokes_discontinuity(B, ħ; order = 1)
                # sign convention pinned here (see docstring)
                @test disc ≈ 2im * big(π) * exp(-1 / ħ) rtol = 1e-30
            end
        end
    end

    @testset "approximant seam (AbstractBorelApproximant)" begin
        Φ = FormalSeries(:euler, 9)
        B = borel(Φ)
        r = pade(B; order = 1)
        ħ = big"0.3"
        # an explicit Padé approximant takes the identical code path
        @test laplace_sum(B, r, ħ) == laplace_sum(B, ħ; order = 1)
        # a hand-rolled approximant satisfying the interface reproduces Euler
        dummy = EulerDummyApprox()
        @test laplace_sum(B, dummy, ħ) ≈ laplace_sum(B, ħ; order = 1) rtol = 1e-40
        # PoleOnRay fires through the seam (pole at −1, ray θ = π, ħ < 0 decays)
        @test_throws PoleOnRay laplace_sum(B, dummy, -ħ; θ = Float64(π))
        # lateral sums and the Stokes discontinuity through the seam
        @test lateral_sum(B, dummy, -ħ; θ = Float64(π), side = :plus) ≈
              lateral_sum(B, -ħ; θ = Float64(π), side = :plus, order = 1) rtol = 1e-40
        @test stokes_discontinuity(B, dummy, -ħ; θ = Float64(π)) ≈
              stokes_discontinuity(B, -ħ; θ = Float64(π), order = 1) rtol = 1e-30
    end

    @testset "Airy end-to-end: prefactor · borel_sum ≈ Ai(1) to ~30 digits" begin
        setprecision(BigFloat, 512) do
            # reference Ai(1) from the Taylor series y'' = zy at 0, with
            # Ai(0) = 3^{-2/3}/Γ(2/3), Ai'(0) = -3^{-1/3}/Γ(1/3),
            # Γ(2/3) = 2π/(√3 Γ(1/3)) — Γ(1/3) from _gamma (validated above)
            f1 = zero(BigFloat); c = big"1.0"
            g1 = zero(BigFloat); d = big"1.0"
            for k in 0:60
                f1 += c; c /= (3k + 3) * (3k + 2)
                g1 += d; d /= (3k + 4) * (3k + 3)
            end
            Γ13 = Resurgence._gamma(big"1.0" / 3)
            Γ23 = 2 * big(π) / (sqrt(big"3.0") * Γ13)
            ai0 = big"3.0"^(-big"2.0" / 3) / Γ23
            dai0 = -big"3.0"^(-big"1.0" / 3) / Γ13
            ai1 = ai0 * f1 + dai0 * g1
            @test Float64(ai1) ≈ SpecialFunctions.airyai(1.0) rtol = 1e-12

            # Borel–Padé–Laplace of the asymptotic series at z = 1: the series is
            # in ħ = z^{-3/2} = 1, with ξ = 2/(3ħ) = 2/3. High-order Padé grows
            # spurious Froissart pole–zero doublets near the positive axis, so
            # integrate on a tilted ray (exact by contour deformation — the true
            # Borel function is analytic off the cut ζ ≤ -4/3).
            Φ = FormalSeries(:airy, 220)
            Φb = FormalSeries(BigFloat.(coefficients(Φ)))
            ħ = big"1.0"
            ξ = 2 / (3ħ)
            prefactor = exp(-ξ) / (2 * sqrt(big(π)))     # z^{1/4} = 1 at z = 1
            s = prefactor * lateral_sum(borel(Φb), ħ; side = :plus, tilt = 3//10,
                                        rtol = big"1e-40")
            @test abs(s - ai1) / abs(ai1) < big"1e-25"
        end
    end
end
