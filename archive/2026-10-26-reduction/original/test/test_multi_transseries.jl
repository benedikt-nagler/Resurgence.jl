@testset "multi_transseries" begin
    # a small rank-2 toy: factorized sectors Φ_{(i,j)} = φ_i · ψ_j over actions (1, 2)
    φ = [FormalSeries([1//1, k]) for k in 0:1]          # φ_0 = 1 + ħ, φ_1 = 1 + 2ħ
    ψ = [FormalSeries([2//1, k]) for k in 0:1]          # ψ_0 = 2 + ħ, ψ_1 = 2 + 3ħ ... (k index)
    toy(secs) = MultiTransseries((1//1, 2//1), secs)
    all_secs = Dict((i, j) => φ[i + 1] * ψ[j + 1] for i in 0:1, j in 0:1)
    F = toy(all_secs)

    @testset "construction and accessors" begin
        @test F isa MultiTransseries{Rational{Int},Rational{Int},2}
        @test actions(F) == (1//1, 2//1)
        @test n_actions(F) == 2
        @test n_sectors(F) == 4
        @test variable(F) == :ħ
        @test eltype(typeof(F)) == Rational{Int}
        @test is_exact(F)
        @test sector(F, (0, 0)) == φ[1] * ψ[1]
        @test sector(F, 1, 1) == φ[2] * ψ[2]            # varargs form
        @test all(iszero, coefficients(sector(F, (5, 5))))   # absent ⇒ zero
        @test weight(F, (1, 1)) == 3                    # 1·1 + 1·2
        @test weight(F, (2, 0)) == 2
        @test sector_indices(F) == [(0, 0), (0, 1), (1, 0), (1, 1)]
        @test !is_resonant(F)                           # weights {0,2,1,3} distinct

        @test_throws InvalidArgument sector(F, (-1, 0))
        @test_throws InvalidArgument sector(F, (0, 0, 0))    # wrong rank
        @test_throws InvalidArgument MultiTransseries((0//1, 2//1),
                                                      Dict((0, 0) => φ[1] * ψ[1]))
        @test_throws InvalidArgument MultiTransseries((1//1, 2//1),
            Dict{NTuple{2,Int},FormalSeries{Rational{Int}}}())
        Ψx = FormalSeries([1//1], :x)
        @test_throws IncompatibleSeries MultiTransseries((1//1, 2//1),
            Dict((0, 0) => φ[1] * ψ[1], (1, 0) => Ψx))
    end

    @testset "arithmetic" begin
        G = toy(Dict((0, 0) => φ[1] * ψ[1]))            # single-sector
        S = F + G
        @test sector(S, (0, 0)) == 2 * (φ[1] * ψ[1])    # overlapping adds
        @test sector(S, (1, 1)) == φ[2] * ψ[2]          # missing passes through (untruncated)
        @test n_terms(sector(S, (1, 1))) == n_terms(φ[2] * ψ[2])
        @test S - G == F
        @test -(-F) == F
        @test sector(3 * F, (1, 0)) == 3 * (φ[2] * ψ[1]) && F * 3 == 3 * F

        # weight-graded Cauchy product truncated to the min bounding box
        P = F * G                                       # G's box is (0,0) ⇒ only (0,0)
        @test sector_indices(P) == [(0, 0)]
        @test sector(P, (0, 0)) == (φ[1] * ψ[1]) * (φ[1] * ψ[1])

        P2 = F * F
        @test sector(P2, (0, 0)) == sector(F, (0, 0)) * sector(F, (0, 0))
        # (F·F)_{(1,1)} = Σ_{i+j=(1,1)} Φ_i Φ_j
        expected = sector(F, (0, 0)) * sector(F, (1, 1)) +
                   sector(F, (1, 0)) * sector(F, (0, 1)) +
                   sector(F, (0, 1)) * sector(F, (1, 0)) +
                   sector(F, (1, 1)) * sector(F, (0, 0))
        @test sector(P2, (1, 1)) == expected

        H = MultiTransseries((1//1, 3//1), Dict((0, 0) => φ[1] * ψ[1]))
        @test_throws IncompatibleSeries F + H           # action mismatch
        K1 = MultiTransseries((1//1,), Dict((0,) => FormalSeries([1//1])))
        @test_throws IncompatibleSeries F + K1          # rank mismatch
    end

    @testset "truncate" begin
        Ft = truncate(F; degree = (0, 1))
        @test sector_indices(Ft) == [(0, 0), (0, 1)]
        Ftd = truncate(F; degree = 0)
        @test sector_indices(Ftd) == [(0, 0)]
        Ftt = truncate(F; terms = 1)
        @test all(n_terms(sector(Ftt, i)) == 1 for i in sector_indices(Ftt))
    end

    @testset "k=1 embedding round-trip (mandatory)" begin
        E = Transseries(:euler, 12)
        M = MultiTransseries(E)
        @test M isa MultiTransseries{Rational{BigInt},Rational{Int},1}
        @test actions(M) == (1//1,)
        @test n_actions(M) == 1
        @test sector(M, (0,)) == sector(E, 0) && sector(M, (1,)) == sector(E, 1)
        @test Transseries(M) == E                       # round-trip

        # arithmetic commutes with the embedding
        E2 = Transseries(1//1, [sector(E, 0), 2 * sector(E, 1)])
        @test MultiTransseries(E) + MultiTransseries(E2) == MultiTransseries(E + E2)
        @test MultiTransseries(E) * MultiTransseries(E2) == MultiTransseries(E * E2)

        @test_throws InvalidArgument Transseries(F)     # rank 2 has no projection
    end

    @testset "summation reductions [rank-1 Euler embedding]" begin
        E = Transseries(:euler, 12)
        M = MultiTransseries(E)
        setprecision(BigFloat, 256) do
            ħ = big"0.2"
            B = borel(sector(E, 0))
            s0 = lateral_sum(B, ħ; side = :plus, order = 1)
            # scalar-σ convenience and tuple form agree, and reduce to the M2 sum
            @test transseries_sum(M, big"0.0", ħ; order = 1) ≈ s0 rtol = 1e-40
            σ = big"0.7"
            @test transseries_sum(M, (σ,), ħ; order = 1) ≈
                  transseries_sum(E, σ, ħ; order = 1) rtol = 1e-40
            @test_throws InvalidArgument transseries_sum(M, (σ, σ), ħ; order = 1)
        end
    end

    @testset "log sectors (M6b)" begin
        # The M6a suite above is the regression gate: the trailing sector-type parameter S
        # must be inert for non-resonant transseries. Here S = LogSeries.
        Λ = LogSeries([FormalSeries([0//1]), FormalSeries([1//1])])       # log ħ
        R = MultiTransseries((1//1, -1//1), Dict((0, 0) => Λ, (1, 1) => Λ))

        @testset "sector type" begin
            @test R isa MultiTransseries{Rational{Int},Rational{Int},2}   # UnionAll still
            @test sector(R, (0, 0)) isa LogSeries
            @test sector(R, (5, 5)) isa LogSeries                         # absent ⇒ zero
            @test log_degree(R) == 1
            @test log_degree(R, (0, 0)) == 1
            @test log_degree(R, 1, 1) == 1
            # a non-resonant transseries keeps FormalSeries sectors and log_degree 0
            @test sector(F, (0, 0)) isa FormalSeries
            @test log_degree(F) == 0
            @test log_degree(F, (1, 1)) == 0
        end

        @testset "mixed operands promote to LogSeries" begin
            P = MultiTransseries((1//1, -1//1), Dict((0, 0) => FormalSeries([1//1])))
            @test sector(P, (0, 0)) isa FormalSeries
            for S in (R + P, P + R, R * P, P * R)
                @test sector(S, (0, 0)) isa LogSeries
            end
            @test log_block(sector(R + P, (0, 0)), 0) == FormalSeries([1//1])
            @test log_block(sector(R + P, (0, 0)), 1) == FormalSeries([1//1])
            # FormalSeries-only operands stay FormalSeries - no silent widening
            @test sector(P + P, (0, 0)) isa FormalSeries
        end

        @testset "graded product: log degrees add, the n-box takes the min" begin
            Q = R * R
            @test log_degree(sector(Q, (0, 0))) == 2              # 1 + 1
            @test log_block(sector(Q, (0, 0)), 2) == FormalSeries([1//1])
            @test Resurgence._bounding_box(Q) == (1, 1)            # min of the boxes
        end

        @testset "truncate reaches both indices" begin
            L2 = LogSeries([FormalSeries([1//1, 2//1]), FormalSeries([3//1, 4//1])])
            T = MultiTransseries((1//1, -1//1), Dict((0, 0) => L2, (1, 1) => L2))
            @test n_terms(sector(truncate(T; terms = 1), (0, 0))) == 1
            @test log_degree(sector(truncate(T; terms = 1), (0, 0))) == 1
            @test n_sectors(truncate(T; degree = (0, 0))) == 1
        end

        @testset "projection to Transseries refuses genuine log blocks" begin
            R1 = MultiTransseries((1//1,), Dict((0,) => Λ))
            @test_throws IncompatibleSeries Transseries(R1)
            # a log-free LogSeries sector projects fine (round-trip preserved)
            R0 = MultiTransseries((1//1,), Dict((0,) => LogSeries(FormalSeries([2//1]))))
            @test sector(Transseries(R0), 0) == FormalSeries([2//1])
        end
    end
end
