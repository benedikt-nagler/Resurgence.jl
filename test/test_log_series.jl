@testset "log_series" begin
    fs(c; kw...) = FormalSeries(c, :ħ; kw...)
    zero1 = fs([0//1])
    one1 = fs([1//1])
    logħ = LogSeries([zero1, one1])                     # log ħ
    Φ = fs([1//1, 2//1, 3//1])                          # 1 + 2ħ + 3ħ²

    @testset "construction and accessors" begin
        L = LogSeries([Φ, one1])
        @test L isa LogSeries{Rational{Int}}
        @test variable(L) == :ħ
        @test eltype(typeof(L)) == Rational{Int}
        @test is_exact(L)
        @test log_degree(L) == 1
        @test log_blocks(L) == [Φ, one1]
        @test log_block(L, 0) == Φ
        @test log_block(L, 1) == one1
        @test log_block(L, 7) == zero1                  # absent ⇒ zero, as Φ[n]
        @test_throws InvalidArgument log_block(L, -1)
        @test n_terms(L) == 3                           # the longest block
        # trailing zero blocks are dropped, so log_degree is exact
        @test log_degree(LogSeries([Φ, zero1, zero1])) == 0
        @test log_degree(LogSeries([zero1, zero1])) == 0
        @test LogSeries([Φ, zero1]) == LogSeries([Φ])
        @test_throws InvalidArgument LogSeries(FormalSeries{Rational{Int}}[])
        @test_throws IncompatibleSeries LogSeries([Φ, FormalSeries([1//1], :z)], :ħ)
    end

    @testset "round-trip with FormalSeries (mandatory)" begin
        @test FormalSeries(LogSeries(Φ)) == Φ
        @test log_degree(LogSeries(Φ)) == 0
        # a genuine log block has no FormalSeries representation
        @test_throws IncompatibleSeries FormalSeries(logħ)
        e = try; FormalSeries(logħ); catch err; err; end
        @test e.field == :log_degree && e.lhs == 1 && e.rhs == 0
    end

    @testset "algebra" begin
        L = LogSeries([Φ, one1])
        @test L + L == LogSeries([Φ + Φ, one1 + one1])
        @test Resurgence._is_zero_series(L - L)
        @test log_degree(L - L) == 0
        @test -L == LogSeries([-Φ, -one1])
        # FormalSeries `+` truncates to the shorter tail, so logħ's 1-term zero block 0
        # shortens Φ — the package-wide convention, inherited blockwise
        @test (L + logħ) == LogSeries([Φ + zero1, one1 + one1])
        @test n_terms(log_block(L + logħ, 0)) == 1
        # log degrees ADD under * (a log block is exact — unlike the truncated ħ direction)
        @test log_degree(logħ * logħ) == 2
        @test log_block(logħ * logħ, 2) == one1
        @test log_degree(L * L) == 2
        @test log_block(L * L, 0) == Φ * Φ
        @test log_block(L * L, 1) == Φ * one1 + one1 * Φ
        @test log_block(L * L, 2) == one1 * one1
        @test logħ * LogSeries(one1) == logħ            # multiplying by 1
        # scalar * promotes
        @test 2 * logħ == LogSeries([zero1, fs([2//1])])
        @test logħ * 2 == 2 * logħ
        @test (0.5 * logħ) isa LogSeries{Float64}
        @test_throws IncompatibleSeries logħ + LogSeries(FormalSeries([1//1], :z))
        @test_throws IncompatibleSeries logħ * LogSeries(FormalSeries([1//1], :z))
    end

    @testset "truncate in both indices" begin
        L = LogSeries([Φ, Φ, Φ])
        @test log_degree(truncate(L; degree = 1)) == 1
        @test log_degree(truncate(L; degree = 5)) == 2   # capped at the true degree
        @test n_terms(truncate(L; terms = 2)) == 2
        @test log_block(truncate(L; terms = 2), 1) == truncate(Φ, 2)
        @test log_degree(truncate(L; degree = 0, terms = 1)) == 0
        @test truncate(L; degree = nothing, terms = nothing) == L
        @test_throws InvalidArgument truncate(L; degree = -1)
    end

    @testset "derivative" begin
        # ∂(log ħ) = 1/ħ  — the defining property of the log basis
        @test derivative(logħ) == LogSeries([fs([1//1]; power_offset = -1)])
        @test log_degree(derivative(logħ)) == 0
        # ħ²∂(log ħ) = ħ: the cokernel identity the whole resonance layer rests on
        @test Resurgence._shift_power(derivative(logħ), 2) ==
              LogSeries([fs([1//1]; power_offset = 1)])
        # exact on rationals, and matching FormalSeries' own derivative on block 0
        @test derivative(LogSeries(Φ)) == LogSeries(derivative(Φ))
        @test derivative(Φ) == fs([0//1, 2//1, 6//1]; power_offset = -1)
        # Leibniz
        for (L, M) in ((logħ, LogSeries([Φ, one1])), (LogSeries([Φ]), logħ * logħ))
            @test derivative(L * M) == derivative(L) * M + L * derivative(M)
        end
        # a fractional power offset differentiates exactly through the Rational path
        H = fs([1//1]; power_offset = 1//2)
        @test derivative(H) == fs([1//2]; power_offset = -1//2)
    end

    @testset "evaluate and the log branch" begin
        ħ = 0.25
        @test evaluate(logħ, ħ) ≈ log(complex(ħ))
        @test evaluate(LogSeries(Φ), ħ) ≈ evaluate(Φ, ħ)
        L = LogSeries([Φ, one1])
        @test evaluate(L, ħ) ≈ evaluate(Φ, ħ) + log(complex(ħ)) * 1
        # log_ħ is a keyword: the caller owns the branch (universal cover of ħ = 0)
        @test evaluate(logħ, ħ; log_x = log(complex(ħ)) + 2π * im) ≈
              log(complex(ħ)) + 2π * im
        @test evaluate(logħ, ħ; log_x = 0) == 0
        # log-free blocks do not see the branch at all
        @test evaluate(LogSeries(Φ), ħ; log_x = 99.0) ≈ evaluate(Φ, ħ)
    end
end
