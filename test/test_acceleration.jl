@testset "acceleration" begin
    @testset "argument validation" begin
        @test_throws InvalidArgument accelerate([1.0])
        @test_throws InvalidArgument accelerate([1.0, 2.0, 3.0]; order = -1)
        @test_throws InvalidArgument accelerate([1.0, 2.0, 3.0]; method = :nope)
        @test_throws InvalidArgument accelerate([1.0, 2.0, 3.0];
                                                method = :shanks, order = 2)
        @test_throws InvalidArgument accelerate([1.0, 2.0, 3.0];
                                                method = :wynn, order = 2)
        @test_throws InvalidArgument accelerate([1.0, 2.0, 3.0];
                                                method = :levin_u, order = 2)
    end

    @testset "shanks exact on geometric approach" begin
        # s_n = 3 + 2 q^n is annihilated by a single Shanks step
        s = [3 + 2 * (big(1) // 2)^n for n in 0:10]
        @test abs(accelerate(s; method = :shanks, order = 1) - 3) < big"1e-70"
        # iterating further keeps the exact limit (zero-denominator guard)
        @test abs(accelerate(s; method = :shanks) - 3) < big"1e-70"
    end

    @testset "wynn column 2 is the Shanks transform" begin
        s = accumulate(+, [(-big(1))^(n + 1) / n for n in 1:12])   # ln 2 partial sums
        @test accelerate(s; method = :wynn, order = 1) ≈
              accelerate(s; method = :shanks, order = 1) rtol = big"1e-60"
    end

    @testset "richardson delegation" begin
        s = [1 + big(1) / n for n in 1:9]
        @test accelerate(s; method = :richardson, order = 3) == richardson(s, 3)
        @test accelerate(s; method = :richardson) == richardson(s, 4)
    end

    @testset "levin on convergent alternating series" begin
        # Leibniz π/4 and ln 2, 20 terms: Levin-u gains ≥ 10 digits over raw
        leib = accumulate(+, [(-big(1))^n / (2n + 1) for n in 0:19])
        ln2 = accumulate(+, [(-big(1))^(n + 1) / n for n in 1:20])
        for (s, limit) in ((leib, big(π) / 4), (ln2, log(big(2))))
            raw = abs(s[end] - limit)
            acc = abs(accelerate(s; method = :levin_u) - limit)
            @test acc < raw * big"1e-10"
        end
        # Wynn also converges (weaker bound)
        @test abs(accelerate(leib; method = :wynn) - big(π) / 4) < big"1e-10"
    end

    @testset "levin sums the divergent Euler series" begin
        Φ = FormalSeries(:euler, 25)
        s = partial_sums(Φ, big"0.1")
        reference = borel_sum(Φ, big"0.1"; order = 1)
        @test abs(accelerate(s; method = :levin_u) - reference) < 1e-6
        @test abs(accelerate(s; method = :levin_t) - reference) < 1e-6
    end

    @testset "partial_sums" begin
        Φ = FormalSeries([1//1, 2, 3], :x)
        @test partial_sums(Φ, 1//2) == [1//1, 2//1, 11//4]
        # power offset enters every term
        Ψ = FormalSeries([1//1, 1], :x; power_offset = 2)
        @test partial_sums(Ψ, 1//2) == [1//4, 3//8]
        # last partial sum agrees with evaluate
        Φe = FormalSeries(:euler, 8)
        @test partial_sums(Φe, big"0.1")[end] ≈ evaluate(Φe, big"0.1")
    end

    @testset "large_order_fit with accelerate methods" begin
        fit = large_order_fit(FormalSeries(:airy, 40); method = :wynn)
        @test fit.A ≈ 4 / 3 atol = 1e-4
        # Euler stays exact through the alternative extrapolators
        fitₑ = large_order_fit(FormalSeries(:euler, 25); method = :shanks, order = 4)
        @test fitₑ.A ≈ 1 atol = 1e-20
    end
end
