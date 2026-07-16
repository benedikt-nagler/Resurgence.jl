@testset "resonance" begin
    fs(c; kw...) = FormalSeries(c, :ħ; kw...)
    zero1 = fs([0//1])
    ħ_series = fs([1//1]; power_offset = 1)             # the series "ħ"

    # ħ²∂ + c applied to a LogSeries — the operator resonant_solve inverts. At c = 0 the
    # `c*u` term is skipped rather than added: `0*u` is a 1-term zero, and FormalSeries `+`
    # truncates to the shorter tail, which would silently eat the ħ²∂u we are checking.
    channel(u, c) = iszero(c) ? Resurgence._shift_power(derivative(u), 2) :
                    Resurgence._shift_power(derivative(u), 2) + c * u
    _is_zero_series_test(x) = Resurgence._is_zero_series(x)

    # Coefficientwise equality of two LogSeries over a window of powers — robust to how a
    # zero block happens to be represented (offset/length), which plain `==` is not.
    function agree(L, M; powers = (-2//1):(1//2):(4//1))
        P = max(log_degree(L), log_degree(M))
        all(Resurgence._coeff_at(log_block(L, p), q) ==
            Resurgence._coeff_at(log_block(M, p), q) for p in 0:P, q in powers)
    end

    @testset "resonance_lattice" begin
        # A = (A, −A): the double-well / Painlevé I case, kernel (1,1) ≽ 0
        @test resonance_lattice((1//1, -1//1)) == [(1, 1)]
        @test resonance_lattice((3//1, -3//1)) == [(1, 1)]
        # commensurate but not opposite: the primitive generator, denominators cleared
        @test resonance_lattice((1//1, -2//1)) == [(2, 1)]
        @test resonance_lattice((2//1, -3//1)) == [(3, 2)]
        # a nonzero kernel need not be ≽ 0: (2,0) and (0,1) collide at weight 2
        @test resonance_lattice((1//1, 2//1)) == [(2, -1)]
        @test resonance_lattice((1//1, 1//1)) == [(1, -1)]     # duplicated actions
        # non-resonant: K = 1 always (A₁ ≠ 0), and ℚ-independent complex actions
        @test isempty(resonance_lattice((1//1,)))
        @test isempty(resonance_lattice((5//2,)))
        @test isempty(resonance_lattice((Complex(1//1, 0//1), Complex(0//1, 1//1))))
        # rank-3 with a 2-dimensional kernel (integral pivots ⇒ provably saturated)
        @test resonance_lattice((1//1, 1//1, -2//1)) == [(1, -1, 0), (2, 0, 1)]
        # sign canonicalization: first nonzero entry positive
        @test resonance_lattice((-1//1, 1//1)) == [(1, 1)]

        @testset "is_resonant on the lattice" begin
            @test is_resonant((1//1, -1//1))
            @test is_resonant((1//1, 2//1))             # collision without a ≽ 0 walk
            @test !is_resonant((1//1,))
            @test !is_resonant((Complex(1//1, 0//1), Complex(0//1, 1//1)))
        end

        @testset "inexact actions refuse to guess" begin
            # Painlevé I's A = 8√6/5 is exactly this case: a ℚ-linear relation between
            # floats is not decidable, so we throw rather than fake a tolerance
            @test_throws InvalidArgument resonance_lattice((1.0, -1.0))
            @test_throws InvalidArgument resonance_lattice((8 * sqrt(6) / 5, -8 * sqrt(6) / 5))
            # …but the caller can supply the kernel
            @test resonance_lattice((1.0, -1.0); lattice = [(1, 1)]) == [(1, 1)]
            @test_throws InvalidArgument resonance_lattice((1.0, -1.0); lattice = [(1, 1, 1)])
            @test_throws InvalidArgument resonance_lattice((1.0, -1.0); lattice = [(0, 0)])
        end
    end

    @testset "resonance_depth" begin
        mk(A, secs) = MultiTransseries(A, secs)
        F = mk((1//1, -1//1), Dict((i, j) => fs([1//1]) for i in 0:3, j in 0:3))
        # k = (1,1) ⇒ depth = min(n₁, n₂), the familiar double-well bound
        for n in ((0, 0), (1, 0), (0, 2), (2, 1), (3, 3), (2, 3))
            @test resonance_depth(F, n) == min(n[1], n[2])
        end
        @test resonance_depth(F, 2, 1) == 1                     # varargs form
        # a commensurate lattice walks in steps of k = (2,1)
        G = mk((1//1, -2//1), Dict((i, j) => fs([1//1]) for i in 0:4, j in 0:4))
        @test resonance_depth(G, (4, 1)) == 1
        @test resonance_depth(G, (4, 2)) == 2
        @test resonance_depth(G, (1, 4)) == 0
        # non-resonant ⇒ depth 0 everywhere
        H = mk((1//1, im * 1//1), Dict((i, j) => fs([1//1]) for i in 0:2, j in 0:2))
        @test resonance_depth(H, (2, 2)) == 0
        # no ≽ 0 generator ⇒ the depth is model-dependent, and we say so
        D = mk((1//1, 2//1), Dict((i, j) => fs([1//1]) for i in 0:2, j in 0:2))
        @test_throws InvalidArgument resonance_depth(D, (2, 2))
        E = mk((1//1, 1//1), Dict((i, j) => fs([1//1]) for i in 0:2, j in 0:2))
        @test_throws InvalidArgument resonance_depth(E, (1, 1))
    end

    @testset "charges: the fiber over a singularity" begin
        F = MultiTransseries((1//1, -1//1),
                             Dict((i, j) => fs([1//1]) for i in 0:3, j in 0:3))
        # resonant: the fiber over ω = A is the whole tower (1,0), (2,1), (3,2)
        @test charges(F, 1//1) == [(1, 0), (2, 1), (3, 2)]
        @test charges(F, 0//1) == [(1, 1), (2, 2), (3, 3)]      # the kernel monomials
        @test charges(F, 2//1) == [(2, 0), (3, 1)]
        @test isempty(charges(F, 17//1))
        # non-resonant: a single charge per location
        G = MultiTransseries((1//1, 3//1), Dict((i, j) => fs([1//1]) for i in 0:1, j in 0:1))
        @test charges(G, 3//1) == [(0, 1)]
        @test charges(G, 1//1) == [(1, 0)]
    end

    @testset "the mechanism: resonant_solve" begin
        @testset "detuned ⇒ unique and log-free" begin
            # (ħ²∂ + c)u = g has a unique power-series solution for c ≠ 0
            for c in (2//1, -1//1, 3//2), g in (fs([1//1]), fs([1//1, 2//1, 3//1]),
                                                fs([0//1, 1//1]))
                u = resonant_solve(g, c)
                @test log_degree(u) == 0                        # NO logs — the theorem
                @test agree(channel(u, c), LogSeries(g))        # round-trip, exact
            end
            # log degree is preserved, not grown, when the source already has blocks
            L = LogSeries([fs([1//1]), fs([2//1])])
            u = resonant_solve(L, 2//1)
            @test log_degree(u) == log_degree(L) == 1
            @test agree(channel(u, 2//1), L)
        end

        @testset "resonant ⇒ log degree grows by one" begin
            # the cokernel identity: ħ²∂ misses ħ¹ exactly, and log ħ is its preimage
            u = resonant_solve(ħ_series, 0//1)
            @test u == LogSeries([zero1, fs([1//1])])           # u = log ħ
            @test log_degree(u) == 1
            @test agree(channel(u, 0//1), LogSeries(ħ_series))
            # a source with no ħ¹ coefficient stays log-free even at zero detuning
            v = resonant_solve(fs([0//1, 0//1, 1//1]), 0//1)    # g = ħ², not ħ
            @test log_degree(v) == 0
            @test agree(channel(v, 0//1), LogSeries(fs([0//1, 0//1, 1//1])))
            @test log_degree(resonant_solve(zero1, 0//1)) == 0
        end

        @testset "the top block is known in closed form" begin
            # u^{[P+1]} = ([ħ¹] g^{[P]}) / (P+1)
            for (P, a) in ((0, 1//1), (0, 5//1), (1, 3//1), (2, 7//2))
                blocks = [zero1 for _ in 0:P]
                blocks[P + 1] = a * ħ_series                    # [ħ¹]g^{[P]} = a
                g = LogSeries(blocks)
                u = resonant_solve(g, 0//1)
                @test log_degree(u) == P + 1
                @test log_block(u, P + 1) == fs([a / (P + 1)])
                @test agree(channel(u, 0//1), g)
            end
        end

        @testset "fractional power offsets never hit the cokernel" begin
            # zero detuning is necessary but NOT sufficient: if power_offset ∉ ℤ the power
            # grid misses ħ¹, so ħ²∂ stays invertible and no log is forced
            g = fs([1//1, 2//1]; power_offset = 1//2)
            u = resonant_solve(g, 0//1)
            @test log_degree(u) == 0
            @test agree(channel(u, 0//1), LogSeries(g))
            @test power_offset(log_block(u, 0)) == -1//2
        end

        @testset "the integration constant is fixed to zero" begin
            # ħ²∂ annihilates constants, so u^{[0]}'s ħ⁰ coefficient is free; we pin it
            u = resonant_solve(fs([0//1, 0//1, 1//1]), 0//1)
            @test log_block(u, 0)[0] == 0
        end
    end

    @testset "forced-resonance oracle: Y = σ₁e^{-A/ħ}·ħ^{σ₁σ₂}" begin
        # Y solves (ħ²∂_ħ − A)Y = ħ·(σ₁σ₂)·Y exactly, with A = (A,−A) and σ₁σ₂ the
        # weight-0 kernel monomial. Grading by (n,m) gives the sector recursion
        #     (ħ²∂_ħ + (n−m−1)A) Φ_{(n,m)} = ħ · Φ_{(n−1,m−1)}                     (†)
        # and, seeded with Φ_{(1,0)} = 1, ħ^{σ₁σ₂} = exp(σ₁σ₂ log ħ) predicts
        #     Φ_{(n+1,n)} = log^n(ħ)/n!,   every other sector zero.
        # This is the M6 milestone oracle: forced resonance, log coefficient in closed form.
        A = 1//1
        # The summation checks below resum Σ_n xⁿ/n! with x = σ₁σ₂·log ħ, so the tail is
        # |x|^{NMAX+1}/(NMAX+1)! — controlled, but it has to clear the rtol. The monodromy
        # case is the binding one (|x| ≈ 0.95 once 2πi joins log ħ): NMAX = 12 puts its
        # tail near 1e-10, comfortably under the 1e-8 asserted. Blocks stay single-term, so
        # the extra orders are cheap.
        NMAX = 12
        Φ = Dict{Tuple{Int,Int},LogSeries{Rational{BigInt}}}()
        Φ[(1, 0)] = LogSeries(FormalSeries([one(Rational{BigInt})], :ħ))
        for tot in 0:(2NMAX + 1), n in 0:tot            # (NMAX+1, NMAX) has total 2NMAX+1
            m = tot - n
            (n, m) == (1, 0) && continue
            src = get(Φ, (n - 1, m - 1), nothing)
            g = src === nothing ? LogSeries(FormalSeries([zero(Rational{BigInt})], :ħ)) :
                Resurgence._shift_power(src, 1)                 # the ħ· on (†)'s RHS
            Φ[(n, m)] = resonant_solve(g, (n - m - 1) * A)
        end

        for k in 0:NMAX
            u = Φ[(k + 1, k)]
            @test log_degree(u) == k
            @test log_block(u, k) == FormalSeries([1 // factorial(big(k))], :ħ)
            # lower blocks vanish: Φ_{(k+1,k)} is exactly log^k(ħ)/k!
            @test all(_is_zero_series_test(log_block(u, p)) for p in 0:(k - 1))
        end

        nonzero = sort([k for (k, v) in Φ if !_is_zero_series_test(v)])
        @test nonzero == [(k + 1, k) for k in 0:NMAX]           # only the m = n−1 diagonal
        @test _is_zero_series_test(Φ[(2, 0)])                   # detuned, empty source
        @test _is_zero_series_test(Φ[(0, 0)])

        @testset "the log-degree bound is saturated" begin
            # log_degree(Φ_n) ≤ resonance_depth(n), and (†) attains it on the diagonal
            F = MultiTransseries((A, -A), Dict(k => v for (k, v) in Φ))
            @test is_resonant(actions(F))
            @test is_resonant(F)                                # support-level too
            for k in 0:NMAX
                n = (k + 1, k)
                @test resonance_depth(F, n) == k
                @test log_degree(F, n) == resonance_depth(F, n)
            end
            @test log_degree(F) == NMAX
        end

        @testset "summation reproduces σ₁e^{-A/ħ}ħ^{σ₁σ₂}" begin
            F = MultiTransseries((A, -A), Dict(k => v for (k, v) in Φ))
            for (σ1, σ2, ħ) in ((0.3, 0.5, 0.4), (1.0, 0.25, 0.2), (0.5, -0.4, 0.3))
                got = transseries_sum(F, (σ1, σ2), ħ)
                want = σ1 * exp(-A / ħ) * complex(ħ)^(σ1 * σ2)
                @test got ≈ want rtol = 1e-8
            end
            # monodromy: log ħ ↦ log ħ + 2πi multiplies Y by exactly e^{2πiσ₁σ₂}.
            # Resonance ⇔ nontrivial monodromy at ħ = 0 — the branch keyword is physics.
            σ1, σ2, ħ = 0.3, 0.5, 0.4
            base = transseries_sum(F, (σ1, σ2), ħ)
            turned = transseries_sum(F, (σ1, σ2), ħ; log_ħ = log(complex(ħ)) + 2π * im)
            @test turned ≈ base * exp(2π * im * σ1 * σ2) rtol = 1e-8
            @test !(turned ≈ base)                              # the branch really matters
        end
    end
end
