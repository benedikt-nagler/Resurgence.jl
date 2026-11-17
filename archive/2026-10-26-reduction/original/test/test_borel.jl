@testset "borel" begin
    @testset "Euler oracle: exactly -1/(1+ζ)" begin
        Φ = FormalSeries(:euler, 8)
        B = borel(Φ)
        # β = 1: (1)_n = n!, b_n = (−1)^{n+1} n!/n! = (−1)^{n+1} - geometric −1/(1+ζ)
        @test B isa BorelSeries{Rational{BigInt}}
        @test B.beta == 1
        @test variable(B) == :ζ
        @test power_offset(B.series) == 0
        @test coefficients(B) == [(-1)^(n + 1) for n in 0:7]
        @test iszero(B.constant_term)
        @test B.source_var == :ħ
    end

    @testset "round trip is exact" begin
        # rationals, offset > 0
        Φ = FormalSeries(:euler, 10)
        @test inverse_borel(borel(Φ)) == Φ
        # offset 0: constant term split off and restored
        Ψ = FormalSeries(:airy, 10)
        B = borel(Ψ)
        @test B.constant_term == 1
        @test inverse_borel(B) == Ψ
        # rational power offset
        Θ = FormalSeries([3//1, 1//1, 4//1]; power_offset = 1//2)
        @test inverse_borel(borel(Θ)) == Θ
        # BigFloat coefficients
        G = FormalSeries([big"1.5", big"2.5"]; power_offset = 2//1)
        @test inverse_borel(borel(G)) == G
        # AbstractAlgebra field elements by duck typing
        F = AbstractAlgebra.QQ
        H = FormalSeries([F(1), F(1, 2), F(1, 6)]; power_offset = 1//1)
        @test inverse_borel(borel(H)) == H
    end

    @testset "type errors & guards" begin
        Φ = FormalSeries(:euler, 5)
        @test_throws InvalidArgument borel(borel(Φ))          # double Borel
        @test_throws InvalidArgument borel(FormalSeries([1//1]; power_offset = -1//1))
        @test_throws InvalidArgument borel(FormalSeries([1//1]))  # pure constant
    end

    @testset "accessors delegate" begin
        B = borel(FormalSeries(:euler, 6))
        @test n_terms(B) == 6
        @test is_exact(B)
        @test B[0] == -1 && B[1] == 1
        @test eltype(typeof(B)) == Rational{BigInt}
    end
end
