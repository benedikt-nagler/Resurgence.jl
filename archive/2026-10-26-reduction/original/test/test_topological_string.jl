R = Rational{BigInt}

# 1/Φ for a unit-grid series with nonzero constant term: d_m = -(1/c_0) Σ_{i≥1} c_i d_{m-i}.
function _inv_series(c::AbstractVector{R})
    d = Vector{R}(undef, length(c))
    d[1] = 1 // c[1]
    for k in 2:length(c)
        d[k] = -sum(c[i + 1] * d[k - i] for i in 1:(k - 1)) // c[1]
    end
    d
end

# (sinh(z/2)/(z/2))^2 and (sin(x/2)/(x/2))^2 on the unit grid, built from the elementary
# Taylor series and squared with the package's own FormalSeries product - no Bernoulli
# numbers anywhere, which is what makes them independent oracles below.
function _half_sinc_squared(n::Integer; alternating::Bool)
    c = zeros(R, n)
    for k in 0:((n - 1) ÷ 2)
        s = alternating && isodd(k) ? -1 : 1
        c[2k + 1] = R(s, big(4)^k * factorial(big(2k + 1)))
    end
    Φ = FormalSeries(c)
    coefficients(Φ * Φ)
end

@testset "topological_string.jl" begin
    @testset "Bernoulli numbers" begin
        B = Resurgence._bernoulli(12)
        @test B[1] == 1
        @test B[2] == -1//2
        @test B[3] == 1//6
        @test B[5] == -1//30
        @test B[7] == 1//42
        @test B[13] == -691//2730
        @test all(iszero, B[4:2:12])                 # B_3 = B_5 = ⋯ = 0
        @test_throws InvalidArgument Resurgence._bernoulli(-1)
    end

    @testset "negative-index polylogarithms" begin
        # Li_0 = Q/(1-Q), Li_{-1} = Q/(1-Q)^2, Li_{-2} = Q(1+Q)/(1-Q)^3
        Q = R(1, 4)
        @test Resurgence._polylog_negative(0, Q) == Q / (1 - Q)
        @test Resurgence._polylog_negative(1, Q) == Q / (1 - Q)^2
        @test Resurgence._polylog_negative(2, Q) == Q * (1 + Q) / (1 - Q)^3
        # Li_{-n}(Q) = Σ_{d≥1} d^n Q^d
        for n in 0:6
            direct = sum(big(d)^n * Q^d for d in 1:400)
            @test Resurgence._polylog_negative(n, Q) ≈ direct rtol = 1e-30
        end
        # Eulerian polynomials sum to n!: this is the exact conifold-gap limit
        # lim_{Q→1} (1-Q)^{n+1} Li_{-n}(Q) = n!
        for n in 0:8
            @test sum(Resurgence._eulerian_polynomial(n)) == factorial(big(n))
        end
        @test_throws InvalidArgument Resurgence._polylog_negative(-1, Q)
    end

    @testset "conifold gap series" begin
        F = conifold_gap_series(4)
        @test power_offset(F) == 2//1
        @test n_terms(F) == 7
        @test is_exact(F)
        # F_g = B_{2g}/(2g(2g-2)) at t = 1, odd powers exactly zero
        @test coefficients(F) == R[-1//240, 0, 1//1008, 0, -1//1440, 0, 1//1056]
        # t enters only as t^{2-2g}
        Ft = conifold_gap_series(3; t = 1//2)
        @test coefficients(Ft) == R[-1//240 * 4, 0, 1//1008 * 16, 0, -1//1440 * 64]
        @test_throws InvalidArgument conifold_gap_series(0)
        @test_throws InvalidArgument conifold_gap_series(3; t = 0//1)
    end

    @testset "exact closed-form Borel transform" begin
        # Reconstruct  g(z) = (1 - 1/s₂ - z²/12)/z³,  s₂ = (sinh(z/2)/(z/2))²,
        # from the sinh Taylor series alone, and demand that it equal the Borel
        # transform of the Bernoulli-built gap series coefficient by coefficient.
        # This pins both the closed form and the λ-grid normalization.
        N = 26
        s2 = _half_sinc_squared(N; alternating = false)
        u = _inv_series(s2)                                    # 1/s₂
        w = -u                                                 # 1 - u - z²/12
        w[1] += 1
        w[3] -= 1//12
        @test all(iszero, w[1:4])                              # starts at z⁴
        @test w[5] == -1//240
        # w_k is the coefficient of z^{k-1}, and g = w/z³, so [z^m]g = w[m+4];
        # B̂(ζ) = g(ζ/t)/t contributes [ζ^m]B̂ = w[m+4]·t^{-m-1}.
        for t in (R(1), R(1, 2), R(3, 2))
            B = borel(conifold_gap_series(9; t))
            b = coefficients(B)                                # true B̂ = ζ · Σ bₖ ζᵏ
            for g in 2:9
                m = 2g - 3                                     # power of ζ (and of z)
                @test b[m] == w[m + 4] * t^(2 - 2g)            # t^{-m-1} = t^{2-2g}
            end
        end

        # the small-|z| Bernoulli branch agrees with the closed form, evaluated in
        # BigFloat where the cancellation is harmless
        closed(z) = -1 / (4z * sinh(z / 2)^2) + 1 / z^3 - 1 / (12z)
        for z in (big(0.99), big(0.5), big(0.1), big(0.7) * im)
            @test conifold_gap_borel(Float64(real(z)) + Float64(imag(z)) * im) ≈
                  ComplexF64(closed(z)) rtol = 1e-13
        end
        @test conifold_gap_borel(0.5; t = 2) ≈ conifold_gap_borel(0.25) / 2 rtol = 1e-14
    end

    @testset "Borel-Padé poles sit at 2πi m t (the gap tower)" begin
        # Padé renders each double pole as a coalescing pair of simple poles, so the
        # oracle is the pair's midpoint. This splitting is itself the signature of the
        # double pole - the Euler and Airy oracles show nothing of the kind.
        for t in (1, 2)
            p = borel_pade_poles(conifold_gap_series(20; t = R(t)); reduce = true)
            top = ComplexF64.(filter(z -> imag(z) > 0, p))
            sort!(top; by = abs)
            @test abs(real(top[1])) < 1e-6
            @test (imag(top[1]) + imag(top[2])) / 2 ≈ 2π * t rtol = 1e-8
            @test abs(top[1] - top[2]) < 1e-5 * abs(top[1])     # coalescing pair
            @test (imag(top[3]) + imag(top[4])) / 2 ≈ 4π * t rtol = 1e-3
        end
    end

    @testset "large order sees |A| = 2πt" begin
        # Gevrey-1 in λ but on the step-2 grid, so the ratio test runs on the genus
        # subsequence: F_g/F_{g+1} · (2g-2)(2g-1) → A².
        for t in (1, 2)
            F = conifold_gap_series(24; t = R(t))
            c = coefficients(F)[1:2:end]                       # F_2, F_3, …
            seq = [Float64((2g - 2) * (2g - 1) * c[g - 1] / c[g]) for g in 2:23]
            A2 = richardson(seq, 6)
            @test sqrt(abs(A2)) ≈ 2π * t rtol = 1e-6
            @test A2 < 0                                       # A² = (2πi t)² < 0
        end
    end

    @testset "Borel-Padé-Laplace reproduces the exact Laplace integral" begin
        B = borel(conifold_gap_series(20))
        r = pade(B; reduce = true)
        for λ in (0.3, 0.6, 1.0)
            approx = laplace_sum(B, r, λ)
            exact, _ = quadgk(ζ -> exp(-ζ / λ) * conifold_gap_borel(ζ), 0.0, Inf;
                              rtol = 1e-12)
            @test approx ≈ exact rtol = 1e-12
        end
    end

    @testset "Stokes discontinuity matches the exact Laurent data" begin
        # Sign pin: with disc = s₊ - s₋ (laplace.jl) a Borel pole c₂/(ζ-ζ_m)² + c₁/(ζ-ζ_m)
        # contributes [2πi c₂/λ - 2πi c₁] e^{-ζ_m/λ}. Every singularity on θ = π/2 has
        # m > 0, so the prediction is a sum over the whole positive tower.
        t = 1
        B = borel(conifold_gap_series(20; t = R(t)))
        r = pade(B; reduce = true)
        for ℓ in (0.5, 0.8, 1.2)
            λ = im * ℓ
            d = stokes_discontinuity(B, r, λ; θ = π / 2, tilt = 0.2)
            pred = sum(1:8) do m
                sd = conifold_gap_stokes_data(m; t)
                (2π * im * sd.double / λ - 2π * im * sd.residue) * exp(-sd.action / λ)
            end
            @test d ≈ pred rtol = 1e-11
            @test abs(real(d)) < 1e-12 * abs(d)                # purely imaginary on the ray
        end
        sd = conifold_gap_stokes_data(2; t = 3)
        @test sd.action ≈ 12π * im
        @test sd.double ≈ -3 / (4π * im)
        @test sd.residue ≈ -1 / (16π^2)
        @test_throws InvalidArgument conifold_gap_stokes_data(0)
    end

    @testset "action bookkeeping" begin
        a = conifold_gap_actions(; t = 2, mmax = 2)
        @test length(a) == 4
        @test a ≈ [4π * im, -4π * im, 8π * im, -8π * im]
        @test issorted(abs.(a))
        @test_throws InvalidArgument conifold_gap_actions(; mmax = 0)

        l = resolved_conifold_actions(; t = 1, mmax = 2, nmax = 1)
        @test length(l) == 12
        @test issorted(abs.(l))
        @test l[1] ≈ 2π                                        # the leading D2 tower
        @test minimum(abs.(l .- 2π * (1 + 2π * im))) < 1e-10    # a D2-D0 lattice point
        @test_throws InvalidArgument resolved_conifold_actions(; nmax = -1)
    end

    @testset "resolved conifold free energies" begin
        Q = R(1, 4)
        F = resolved_conifold_series(4; Q)
        @test power_offset(F) == 2//1
        @test n_terms(F) == 7
        @test is_exact(F)
        @test all(iszero, coefficients(F)[2:2:end])
        # F_g = (-1)^{g-1} B_{2g}/(2g(2g-2)!) Li_{3-2g}(Q), checked against the
        # Gopakumar-Vafa resummation built from the *sine* series (no Bernoulli numbers):
        # 1/(2 sin(x/2))² = Σ_g κ_g x^{2g-2} with κ_g = [x^{2g}] (sin(x/2)/(x/2))^{-2}.
        v = _inv_series(_half_sinc_squared(24; alternating = true))
        for g in 2:5
            κ = v[2g + 1]
            li = sum(big(d)^(2g - 3) * Q^d for d in 1:400)
            @test coefficients(F)[2g - 3] ≈ κ * li rtol = 1e-25
        end
        @test_throws InvalidArgument resolved_conifold_series(0; Q)
        @test_throws InvalidArgument resolved_conifold_series(3; Q = 1//1)
    end

    @testset "the conifold gap is the Q → 1 limit of the resolved conifold" begin
        # (1-Q)^{2g-2} F_g^{RC}(Q) → (-1)^{g-1} B_{2g}/(2g(2g-2)) = the gap coefficient
        # at t ↦ i t, since i^{2-2g} = (-1)^{g-1}.
        Q = R(10^12 - 1, 10^12)
        RC = coefficients(resolved_conifold_series(5; Q))
        gap = coefficients(conifold_gap_series(5))
        for g in 2:6
            scaled = RC[2g - 3] * (1 - Q)^(2g - 2)
            @test Float64(scaled) ≈ (isodd(g) ? 1 : -1) * Float64(gap[2g - 3]) rtol = 1e-9
        end
    end

    @testset "resolved conifold Borel poles land on the derived lattice" begin
        # t = 1, so the n = 0 tower ζ = 2π m t sits on the real axis - the Stokes ray
        # θ = 0 whose jumps carry the Donaldson-Thomas invariants.
        t = 1.0
        p = borel_pade_poles(resolved_conifold_series(20; Q = exp(-big(1.0))); reduce = true)
        right = ComplexF64.(filter(z -> real(z) > 0, p))
        sort!(right; by = abs)
        @test abs(imag(right[1])) < 1e-6
        @test (real(right[1]) + real(right[2])) / 2 ≈ 2π * t rtol = 1e-8
        @test (real(right[3]) + real(right[4])) / 2 ≈ 4π * t rtol = 1e-3
        lat = resolved_conifold_actions(; t, mmax = 2, nmax = 0)
        for z in right[1:4]
            @test minimum(abs.(lat .- z)) < 1e-2 * abs(z)
        end
    end
end
