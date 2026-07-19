@testset "painleve" begin
    R = Rational{BigInt}
    cf(Φ, p) = Resurgence._coeff_at(Φ, p // 1)

    @testset "exact perturbative seed (GIKM eq. 2.1)" begin
        Φ = FormalSeries(:painleve1, 6)
        a = coefficients(Φ)
        @test eltype(a) == R
        @test a[1] == 1                       # a₀
        @test a[2] == -1 // 48                # a₁
        @test a[3] == -49 // 4608             # a₂
        # the recursion a_{n+1} = ((25n²−1)/48)aₙ − ½Σ_{m=1}^{n} aₘ a_{n+1−m} holds
        for n in 0:3
            lhs = a[n + 2]
            conv = sum((a[m + 1] * a[n + 1 - m + 1] for m in 1:n); init = zero(R))
            rhs = R(25 * n^2 - 1, 48) * a[n + 1] - conv // 2
            @test lhs == rhs
        end
        @test_throws InvalidArgument FormalSeries(:painleve1, 0)
    end

    @testset "exact action, both conventions" begin
        @test painleve1_action_squared() == 192 // 25               # GIKM
        @test painleve1_action_squared(convention = :gikm) == 192 // 25
        @test painleve1_action_squared(convention = :string) == 384 // 25
        @test painleve1_action() ≈ 8 * sqrt(big(3)) / 5
        @test painleve1_action(convention = :string) ≈ 8 * sqrt(big(6)) / 5
        @test painleve1_action()^2 ≈ 192 / 25
        @test_throws InvalidArgument painleve1_action_squared(convention = :bogus)
    end

    @testset "ODE residual vanishes (u₀ solves −u″/6 + u² = z)" begin
        for conv in (:gikm, :string)
            u0 = Resurgence._painleve1_u0(7; T = R, convention = conv)
            res = Resurgence._painleve1_residual(u0; convention = conv)
            # exact rational: every determined coefficient is zero
            @test all(iszero, coefficients(res)[1:40])
        end
    end

    @testset "sector-recursion solver" begin
        F = painleve1(6; sectors = 2)
        @test action(F) ≈ painleve1_action()
        @test n_sectors(F) == 2
        # sector 0 is the perturbative series cast to BigFloat
        Φ0 = sector(F, 0)
        @test cf(Φ0, -2) ≈ 1
        @test cf(Φ0, 8) ≈ -1 / 48            # a₁ at power −2 + 10
        @test cf(Φ0, 18) ≈ -49 / 4608        # a₂

        A = painleve1_action()
        u0 = Resurgence._painleve1_u0(6; T = BigFloat, pad = 75 - (10 * 5 + 1))

        @testset "leading balance fixes A² = 192/25" begin
            # the k=1 transmonomial is balanced iff the s^{-2} diagonal of L₁ vanishes
            β1 = Resurgence._painleve1_beta(1)
            unitβ = FormalSeries([i == 1 ? one(BigFloat) : zero(BigFloat) for i in 1:80],
                                 :s; power_offset = β1)
            Lu = Resurgence._painleve1_apply_L(unitβ, 1, A, u0, 1 // 6)
            @test abs(cf(Lu, β1 - 2)) < 1e-40         # s^{-2} coefficient ≈ 0
            # a wrong action leaves it nonzero
            Lbad = Resurgence._painleve1_apply_L(unitβ, 1, A + 1, u0, 1 // 6)
            @test abs(cf(Lbad, β1 - 2)) > 1e-3
        end

        @testset "one-instanton sector is self-consistent (L₁ψ₁ ≈ 0)" begin
            ψ1 = sector(F, 1)
            @test cf(ψ1, Resurgence._painleve1_beta(1)) ≈ 1     # normalization u_{0,1}=1
            L1 = Resurgence._painleve1_apply_L(ψ1, 1, A, u0, 1 // 6)
            # homogeneous: every determined coefficient of L₁ψ₁ vanishes
            for p in (17 // 2, 27 // 2, 37 // 2, 47 // 2)
                @test abs(cf(L1, p)) < 1e-30
            end
        end
    end

    @testset "large order recovers |A|² = 192/25 (Gevrey-2)" begin
        # u₀'s coefficients grow like Γ(2n)/A^{2n}: aₙ/aₙ₋₁ → 4n²/A², so
        # rₙ = 4n² aₙ₋₁/aₙ → A². Richardson-extrapolate rₙ in 1/n.
        a = coefficients(FormalSeries(:painleve1, 45))
        r = [4 * big(n)^2 * a[n] / a[n + 1] for n in 2:44]   # rₙ, n = 2 … 44
        A2 = richardson(r, 6)
        @test A2 ≈ 192 / 25 rtol = 1e-2
    end

    @testset "resonant (A, −A) structure needs lattice = [(1,1)]" begin
        A = painleve1_action()
        acts = (A, -A)
        # inexact action ⇒ resonance layer refuses to guess, must be told the kernel
        @test_throws InvalidArgument resonance_lattice(acts)
        @test resonance_lattice(acts; lattice = [(1, 1)]) == [(1, 1)]
        @test is_resonant(acts; lattice = [(1, 1)])
        # the double-well / Painlevé I kernel: depth = min(n₁, n₂) on the diagonal
        fs1(x) = FormalSeries([x], :ħ)
        mt = MultiTransseries(acts, Dict((i, j) => fs1(1.0 + 0im) for i in 0:3, j in 0:3))
        for n in ((1, 1), (2, 1), (3, 3), (2, 3))
            @test resonance_depth(mt, n; lattice = [(1, 1)]) == min(n[1], n[2])
        end
        # the diagonal n₁ = n₂ collapses onto the perturbative (weight-0) ray - the
        # resonance that forces logs. (`charges` over an inexact action is skipped: it
        # tests weight equality bit-exactly, which 3A−2A is not; the exact-lattice
        # instanton tower is covered in test_resonance.jl.)
        @test weight(mt, (2, 2)) ≈ 0 atol = 1e-30
        @test weight(mt, (1, 0)) ≈ A
    end

    @testset "TransseriesSolveError" begin
        e = TransseriesSolveError("boom")
        @test e isa ResurgenceError
        @test occursin("boom", sprint(showerror, e))
        # a sector solve that returns the wrong variable is rejected
        bad = (sofar, k) -> FormalSeries([1.0], :wrongvar)
        seed = FormalSeries([1.0], :s)
        @test_throws TransseriesSolveError transseries_solve(seed, 1.0, bad; sectors = 2)
    end
end
