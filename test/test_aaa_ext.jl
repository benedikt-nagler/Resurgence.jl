@testset "BaryRational extension" begin
    # the pre-load stub is asserted in test_extension_stubs.jl, which runs first
    @eval using BaryRational

    @testset "Euler: pole hunter and Laplace seam" begin
        Φ = FormalSeries(:euler, 21)
        B = borel(Φ)
        a = aaa_borel(B)
        @test a isa AbstractBorelApproximant
        # AAA runs in Float64 — a pole/branch-cut hunter, accurate to ~1e-7, not a
        # high-precision summation tool (that is Padé/BigFloat's job). It emits
        # several spurious (Froissart) poles too, so locate the physical one at −1.
        pol = poles(a)
        @test minimum(abs.(pol .+ 1)) < 1e-7
        rs = residues(a)
        _, res = rs[argmin(abs(first(p) + 1) for p in rs)]
        @test res ≈ -1 atol = 1e-5           # residue of −1/(1+ζ) at −1
        # summation through the approximant seam
        ħ = 0.3
        @test laplace_sum(B, a, ħ) ≈ laplace_sum(B, ħ; order = 1) rtol = 1e-6
    end

    @testset "two-pole hunter feeds subtract_singularity" begin
        s = 1 // 3
        Φ2 = FormalSeries([Rational{BigInt}(factorial(big(n))) *
                           ((-1)^n + s * (-1 // 2)^n) for n in 0:24];
                          power_offset = 1)
        B2 = borel(Φ2)
        a2 = aaa_borel(B2)
        ζs = poles(a2)
        @test minimum(abs.(ζs .+ 1)) < 1e-8
        @test minimum(abs.(ζs .+ 2)) < 1e-6
        # the AAA pole works as the peeling input (residue route is exact once ζ₀
        # is supplied, so c is recovered far better than the AAA pole itself)
        out = subtract_singularity(B2, ζs[argmin(abs.(ζs .+ 1))])
        @test out.c ≈ 1 atol = 1e-6
    end

    @testset "aaa_approximant from explicit samples" begin
        zs = [ComplexF64(0.3 * cospi(2k / 40), 0.3 * sinpi(2k / 40)) for k in 0:39]
        fs = [-1 / (1 + z) for z in zs]
        a = aaa_approximant(zs, fs)
        @test a(0.1) ≈ -1 / 1.1 rtol = 1e-10
        @test poles(a)[1] ≈ -1 atol = 1e-10
        @test_throws InvalidArgument aaa_approximant(zs, fs[1:3])
    end
end
