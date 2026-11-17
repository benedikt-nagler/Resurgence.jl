@testset "large_order" begin
    @testset "richardson" begin
        # s_n = 1 + 1/n: one sweep eliminates the 1/n term exactly
        seq = [1.0 + 1 / n for n in 1:10]
        @test richardson(seq, 1) ≈ 1.0 atol = 1e-14
        # s_n = π + 2/n + 3/n²
        seq2 = [Float64(π) + 2 / n + 3 / n^2 for n in 1:20]
        @test richardson(seq2, 3) ≈ π atol = 1e-10
        @test richardson([5.0], 0) == 5.0
        @test_throws InvalidArgument richardson(seq, 10)
        @test_throws InvalidArgument richardson(seq, -1)
    end

    @testset "Euler oracle: (A, b, S) = (1, 0, 1)" begin
        Φ = FormalSeries(:euler, 30)
        fit = large_order_fit(Φ)
        @test fit.A isa BigFloat                    # precision follows coefficients
        @test fit.A ≈ 1 atol = 1e-25
        @test fit.b ≈ 0 atol = 1e-25
        @test fit.S ≈ 1 atol = 1e-25
    end

    @testset "Airy oracle: A = 4/3" begin
        fit = large_order_fit(FormalSeries(:airy, 60); order = 4)
        @test fit.A ≈ 4 / big"3.0" atol = 1e-7
    end

    @testset "guards" begin
        @test_throws InvalidArgument large_order_fit(FormalSeries(:euler, 4))
        @test_throws InvalidArgument large_order_fit(FormalSeries([1//1, 0//1, 1//1,
                                                                   1//1, 1//1, 1//1,
                                                                   1//1, 1//1]))
    end
end
