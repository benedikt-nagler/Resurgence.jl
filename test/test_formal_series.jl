@testset "formal_series" begin
    @testset "construction & accessors" begin
        Φ = FormalSeries([1//1, 2//1, 3//1])
        @test Φ isa FormalSeries{Rational{Int}} && Φ isa AbstractSeries
        @test n_terms(Φ) == 3
        @test coefficients(Φ) == [1//1, 2//1, 3//1]
        @test variable(Φ) == :ħ
        @test power_offset(Φ) == 0
        @test is_exact(Φ)
        @test !is_exact(FormalSeries([1.0, 2.0]))
        @test is_exact(FormalSeries([1 + 2im, 3 + 0im]))
        @test !is_exact(FormalSeries([big"1.0" + 0im]))
        Ψ = FormalSeries([1.0, 2.0], :x; power_offset = 1//2)
        @test variable(Ψ) == :x && power_offset(Ψ) == 1//2
        @test_throws InvalidArgument FormalSeries(Int[])
        # immutability: constructor copies its input
        v = [1, 2]
        Θ = FormalSeries(v)
        v[1] = 99
        @test coefficients(Θ) == [1, 2]
    end

    @testset "indexing" begin
        Φ = FormalSeries([5//1, 7//1])
        @test Φ[0] == 5 && Φ[1] == 7
        @test Φ[10] == 0                       # beyond truncation → zero
        @test_throws InvalidArgument Φ[-1]
    end

    @testset "arithmetic" begin
        Φ = FormalSeries([1//1, 1//1, 1//1])
        Ψ = FormalSeries([0//1, 1//1, 2//1])
        @test coefficients(Φ + Ψ) == [1//1, 2//1, 3//1]
        @test coefficients(Φ - Ψ) == [1//1, 0//1, -1//1]
        @test coefficients(-Φ) == [-1//1, -1//1, -1//1]
        @test coefficients(3 * Φ) == [3//1, 3//1, 3//1]
        @test Φ * Ψ == Ψ * Φ
        # Cauchy product: (1+x+x²)(x+2x²) = x + 3x² + 3x³ → truncated to 3 terms
        @test coefficients(Φ * Ψ) == [0//1, 1//1, 3//1]
        # offsets add under *, align under + (integer shift)
        A = FormalSeries([1//1, 1//1]; power_offset = 1//1)
        @test power_offset(A * A) == 2
        S = A + FormalSeries([1//1, 1//1, 1//1]; power_offset = 2//1)
        @test power_offset(S) == 1 && coefficients(S) == [1//1, 2//1]
        # incompatibilities are typed errors
        @test_throws IncompatibleSeries Φ + FormalSeries([1//1], :x)
        @test_throws IncompatibleSeries A + FormalSeries([1//1]; power_offset = 1//2)
    end

    @testset "truncate & evaluate" begin
        Φ = FormalSeries([1//1, 2//1, 3//1])
        @test coefficients(truncate(Φ, 2)) == [1//1, 2//1]
        @test n_terms(truncate(Φ, 3)) == 3
        @test_throws InvalidArgument truncate(Φ, 0)
        @test_throws InvalidArgument truncate(Φ, 4)
        @test evaluate(Φ, 1//2) == 1 + 2 * (1//2) + 3 * (1//4)
        G = FormalSeries([1.0, 1.0]; power_offset = 1//2)
        @test evaluate(G, 4.0) ≈ 2.0 * 5.0     # √4 (1 + 4)
    end

    @testset "generic coefficient types (duck typing)" begin
        # Complex{BigFloat}
        Φ = FormalSeries([big"1.0" + im * big"2.0", big"3.0" + 0im])
        @test (Φ + Φ)[0] == big"2.0" + im * big"4.0"
        # AbstractAlgebra ring elements work without Resurgence depending on AA
        F = AbstractAlgebra.GF(7)
        Ψ = FormalSeries([F(1), F(2), F(3)])
        @test is_exact(Ψ)
        @test coefficients(Ψ + Ψ) == [F(2), F(4), F(6)]
        @test coefficients(Ψ * Ψ) == [F(1), F(4), F(3)]   # (1+2x+3x²)² mod 7
    end
end
