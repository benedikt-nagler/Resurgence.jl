# Multi-parameter alien calculus oracles. Sign conventions inherit from alien.jl
# (Euler-pinned S₁ = −2πi, s₋(F)(σ) = s₊(F)(σ + S₁)).

# one-parameter Euler transseries with Borel pole at ζ = a (a > 0 on the ray θ = 0):
# Φ₀ = Σ_k (k!/a^k) ħ^{k+1} has reduced Borel Σ (ζ/a)^k = 1/(1 − ζ/a), pole at ζ = a.
function _euler_at(a, n)
    c = [factorial(big(k)) // big(a)^k for k in 0:(n - 1)]
    Φ0 = FormalSeries(c; power_offset = 1//1)
    Φ1 = FormalSeries([one(Rational{BigInt})])
    Transseries(Rational{BigInt}(a), [Φ0, Φ1])
end

@testset "multi_alien" begin
    @testset "k=1 embedding: bridge equation projects to the M2 result" begin
        E = Transseries(:euler, 12)
        ME = MultiTransseries(E)
        setprecision(BigFloat, 256) do
            S1 = stokes_constant(borel(sector(E, 0)); order = 1)
            @test S1 ≈ -2im * big(π) rtol = 1e-40
            D = alien_derivative(ME, (1,); stokes = S1)
            @test sector(D, (0,)) == S1 * sector(E, 1)            # Δ_A Φ₀ = S₁ Φ₁
            @test Transseries(D) == alien_derivative(E; stokes = S1)   # projects to M2
            # default charge e₁ agrees with the explicit (1,)
            @test alien_derivative(ME; stokes = S1) == D
            # Δ_A annihilates the one-instanton tail
            D2 = alien_derivative(D, (1,); stokes = S1)
            @test all(all(iszero, coefficients(Φ)) for Φ in values(sectors(D2)))
        end
    end

    @testset "one-parameter consistency: operator σ-shift = M2 connection [Euler]" begin
        E = Transseries(:euler, 12)
        ME = MultiTransseries(E)
        setprecision(BigFloat, 256) do
            S1 = stokes_constant(borel(sector(E, 0)); order = 1)
            for (σ, ħ) in ((big"0.0", big"0.2"), (big"0.7", big"0.1"),
                           (big"-1.3", big"0.3"))
                lhs = transseries_sum(ME, σ, ħ; side = :minus, order = 1)
                σ′ = stokes_automorphism(ME, σ; stokes = S1)          # scalar K=1 form
                @test σ′ == σ + S1
                rhs = transseries_sum(ME, σ′, ħ; side = :plus, order = 1)
                @test lhs ≈ rhs rtol = 1e-30
            end
        end
    end

    @testset "two-parameter toy: exactly-solvable alien algebra" begin
        # factorized sectors Φ_{(i,j)} = φ_i · ψ_j over the (weight-distinct) lattice (1,2)
        φ = [FormalSeries([1//1, k + 1]) for k in 0:1]
        ψ = [FormalSeries([1//1, k + 2]) for k in 0:1]
        G = MultiTransseries((1//1, 2//1),
                             Dict((i, j) => φ[i + 1] * ψ[j + 1] for i in 0:1, j in 0:1))
        Sd = Dict((1, 0) => 2//1, (0, 1) => 3//1)

        # closed-form bridge Δ_{ℓ·A} Φ_n = S_ℓ (n_j + c) Φ_{n+ℓ}
        D10 = alien_derivative(G, (1, 0); stokes = Sd)
        @test sector(D10, (0, 0)) == 2 * sector(G, (1, 0))       # μ = 0 + 1
        @test sector(D10, (0, 1)) == 2 * sector(G, (1, 1))
        @test sector_indices(D10) == [(0, 0), (0, 1)]
        D01 = alien_derivative(G, (0, 1); stokes = Sd)
        @test sector(D01, (0, 0)) == 3 * sector(G, (0, 1))

        # pointed derivatives for distinct fundamental directions commute (exact)
        P12 = pointed_alien_derivative(pointed_alien_derivative(G, (0, 1); stokes = Sd),
                                       (1, 0); stokes = Sd)
        P21 = pointed_alien_derivative(pointed_alien_derivative(G, (1, 0); stokes = Sd),
                                       (0, 1); stokes = Sd)
        @test P12 == P21
        @test sector(P12, (1, 1)) == 6 * sector(G, (1, 1))       # 2·3·μ·μ = 6

        # mixed charge needs an explicit multiplicity
        @test_throws InvalidArgument alien_derivative(G, (1, 1); stokes = Sd)
        Dmix = alien_derivative(G, (1, 1); stokes = Dict((1, 1) => 6//1),
                                multiplicity = (n, ℓ) -> (n[1] + 1) * (n[2] + 1))
        @test sector(Dmix, (0, 0)) == 6 * sector(G, (1, 1))
        # missing Stokes constant for the requested charge
        @test_throws InvalidArgument alien_derivative(G, (1, 0); stokes = Dict((0, 1) => 3//1))
    end

    @testset "Stokes automorphism as an operator (nilpotent, exact)" begin
        # rank-1 constant-sector toy so 𝒩² = 0: 𝔖F = F + S·Δ_{e₁}F
        F = MultiTransseries((1//1,), Dict((0,) => FormalSeries([5//1]),
                                           (1,) => FormalSeries([1//1])))
        SF = stokes_automorphism(F; stokes = 2//1)
        @test sector(SF, (0,)) == FormalSeries([7//1])           # 5 + 2·1·1
        @test sector(SF, (1,)) == FormalSeries([1//1])
        @test stokes_automorphism(F; stokes = 2//1, order = 1) == SF
    end

    @testset "two-parameter connection formula, per direction" begin
        e1 = _euler_at(1, 12)
        e2 = _euler_at(2, 12)
        Fa = MultiTransseries((1//1, 2//1),
                              Dict((0, 0) => sector(e1, 0), (1, 0) => sector(e1, 1)))
        Fb = MultiTransseries((1//1, 2//1),
                              Dict((0, 0) => sector(e2, 0), (0, 1) => sector(e2, 1)))
        setprecision(BigFloat, 256) do
            S1a = stokes_constant(borel(sector(e1, 0)); order = 1, location = 1)
            Sa = Dict((1, 0) => S1a, (0, 1) => zero(S1a))        # only direction 1 shifts
            for (σ1, σ2, ħ) in ((big"0.0", big"0.4", big"0.2"),
                                (big"0.7", big"-1.1", big"0.15"))
                lhs = transseries_sum(Fa, (σ1, σ2), ħ; side = :minus, order = 1)
                σ′ = stokes_automorphism(Fa, (σ1, σ2); stokes = Sa)
                @test σ′[1] ≈ σ1 + S1a rtol = 1e-40
                @test σ′[2] ≈ σ2 rtol = 1e-40
                rhs = transseries_sum(Fa, σ′, ħ; side = :plus, order = 1)
                @test lhs ≈ rhs rtol = 1e-30
            end

            S2b = stokes_constant(borel(sector(e2, 0)); order = 1, location = 2)
            Sb = Dict((0, 1) => S2b, (1, 0) => zero(S2b))        # only direction 2 shifts
            for (σ1, σ2, ħ) in ((big"0.3", big"0.0", big"0.2"),
                                (big"-0.9", big"0.6", big"0.15"))
                lhs = transseries_sum(Fb, (σ1, σ2), ħ; side = :minus, order = 1)
                σ′ = stokes_automorphism(Fb, (σ1, σ2); stokes = Sb)
                @test σ′[1] ≈ σ1 rtol = 1e-40
                @test σ′[2] ≈ σ2 + S2b rtol = 1e-40
                rhs = transseries_sum(Fb, σ′, ħ; side = :plus, order = 1)
                @test lhs ≈ rhs rtol = 1e-30
            end

            # median summation is real on the embedded Euler (mirrors the M2 oracle)
            E = Transseries(:euler, 12)
            ME = MultiTransseries(E)
            S1 = stokes_constant(borel(sector(E, 0)); order = 1)
            med = median_sum(ME, big"0.0", big"0.1"; stokes = S1, order = 1)
            @test abs(imag(med)) < big"1e-40"
        end
    end
end
