# Named course-canon series via exact rational recursions, mirroring the
# `Quiver(:A, 3)` named-constructor style of ClusterAlgebras.jl.

"""
    FormalSeries(name::Symbol, n::Integer; var = :ħ) -> FormalSeries{Rational{BigInt}}

Named course-canon series with `n` exact rational coefficients:

- `:euler` — the Euler-equation series ``\\sum_{k≥0} k!\\,(-ħ)^{k+1}
  = -ħ + ħ^2 - 2ħ^3 + …`` (power offset 1); its Borel transform is exactly
  ``-1/(1+ζ)``.
- `:airy` — the Airy asymptotic series in the WKB variable ``ħ = z^{-3/2}``:
  ``\\sum_{k≥0} (-1)^k u_k (\\tfrac32 ħ)^k`` with ``u_0 = 1``,
  ``u_{k+1} = u_k (6k+5)(6k+1) / (72(k+1))``, the series ``Φ(ħ)`` in
  ``\\mathrm{Ai}(z) ∼ e^{-ξ} Φ(ħ) / (2\\sqrt{π} z^{1/4})``,
  ``ξ = \\tfrac23 z^{3/2} = \\tfrac{2}{3ħ}``. Its Borel singularity sits at
  ``ζ = -4/3``, the Airy action.
- `:airy_bi` — the second Airy sector ``\\sum_{k≥0} u_k (\\tfrac32 ħ)^k`` (the same
  ``u_k`` without the alternating sign), the series ``Ψ(ħ)`` in
  ``\\mathrm{Bi}(z) ∼ e^{ξ} Ψ(ħ) / (\\sqrt{π} z^{1/4})``; Borel singularity at
  ``ζ = +4/3``. This is the one-instanton sector of the Airy transseries:
  ``Δ_A Φ = S_1 Ψ`` (see [`alien_derivative`](@ref)).
- `:quartic` — the Bender–Wu perturbative ground-state energy of the quartic
  anharmonic oscillator ``H = p^2/2 + x^2/2 + g\\,x^4``:
  ``E_0(g) = \\tfrac12 + \\tfrac34 g - \\tfrac{21}{8} g^2 + \\tfrac{333}{16} g^3 - …``
  (exact rationals from the wavefunction recursion; large order
  ``E_n \\sim -(\\sqrt{6}/π^{3/2})\\,(-3)^n\\,Γ(n+\\tfrac12)``, instanton action 1/3).
- `:phi4` — the zero-dimensional φ⁴ partition function
  ``Z(g) = \\tfrac{1}{\\sqrt{2π}} ∫_{-∞}^{∞} e^{-x^2/2 - g x^4}\\,dx
  = \\sum_{n≥0} \\tfrac{(-1)^n (4n-1)!!}{n!}\\, g^n``; the ``n``-th coefficient counts
  the Wick contractions (vacuum Feynman diagrams) at order ``n``. Closed form
  ``a_n = (-4)^n\\,Γ(2n+\\tfrac12)/(\\sqrt{π}\\, n!)``, large order
  ``a_n \\sim (-16)^n\\,Γ(n+\\tfrac12)/\\sqrt{π}``; Borel singularity at
  ``ζ = -1/16``, the action of the nontrivial saddle of ``x^2/2 + g x^4``.
- `:exp` — the exponential ``\\sum_{k≥0} ħ^k / k!`` (convergent sanity example).
- `:painleve1` — the perturbative (0-instanton) series of Painlevé I in the
  Garoufalidis–Its–Kapaev–Mariño normalization ``-\\tfrac16 u'' + u^2 = z``
  (arXiv:1002.3634, eq. 1.2): ``u_0(z) = \\sqrt{z}\\,\\sum_{n≥0} a_n z^{-5n/2}`` with
  ``a_0 = 1`` and the exact rational recursion (eq. 2.1)
  ``a_{n+1} = \\tfrac{25n^2-1}{48}\\,a_n - \\tfrac12 \\sum_{m=1}^{n} a_m a_{n+1-m}``,
  giving ``a_1 = -1/48``, ``a_2 = -49/4608``, …. This returns the bare fluctuation
  ``\\sum_n a_n x^n`` (``x = z^{-5/2}``, power offset 0); the ``\\sqrt{z}`` prefactor and
  the transmonomial ``e^{-kA z^{5/4}}`` (action ``A = 8\\sqrt3/5``) are supplied by
  [`painleve1`](@ref). Its Borel singularities sit at the instanton actions ``±A``.
"""
function FormalSeries(name::Symbol, n::Integer; var::Symbol = :ħ)
    n ≥ 1 || throw(InvalidArgument("need n ≥ 1 coefficients, got $n"))
    if name === :euler
        # Σ_{k≥0} k! (−ħ)^{k+1} = ħ Σ_{k≥0} k! (−1)^{k+1} ħ^k
        c = Vector{Rational{BigInt}}(undef, n)
        c[1] = -1
        for k in 1:(n - 1)
            c[k + 1] = -k * c[k]
        end
        FormalSeries(c, var; power_offset = 1//1)
    elseif name === :airy || name === :airy_bi
        # :airy — coefficients (−1)^k u_k (3/2)^k: c_{k+1} = −c_k (6k+5)(6k+1)/(48(k+1));
        # :airy_bi — the same without the sign.
        s = name === :airy ? -1 : 1
        c = Vector{Rational{BigInt}}(undef, n)
        c[1] = 1
        for k in 0:(n - 2)
            c[k + 2] = s * c[k + 1] * Rational{BigInt}((6k + 5) * (6k + 1), 48 * (k + 1))
        end
        FormalSeries(c, var)
    elseif name === :quartic
        FormalSeries(_bender_wu(n), var)
    elseif name === :phi4
        # a_n = (−1)^n (4n−1)!!/n!: a_{k+1} = −a_k (4k+1)(4k+3)/(k+1)
        c = Vector{Rational{BigInt}}(undef, n)
        c[1] = 1
        for k in 0:(n - 2)
            c[k + 2] = -c[k + 1] * Rational{BigInt}((4k + 1) * (4k + 3), k + 1)
        end
        FormalSeries(c, var)
    elseif name === :exp
        c = Vector{Rational{BigInt}}(undef, n)
        c[1] = 1
        for k in 1:(n - 1)
            c[k + 1] = c[k] // k
        end
        FormalSeries(c, var)
    elseif name === :painleve1
        FormalSeries(_painleve1(n), var)
    else
        throw(InvalidArgument("unknown named series :$name (expected :euler, :airy, " *
                              ":airy_bi, :quartic, :phi4, :exp, or :painleve1)"))
    end
end

# Painlevé I perturbative recursion for −c·u″ + u² = z, u₀ = √z Σ aₙ z^{−5n/2}: matching
# powers gives, for k ≥ 1,
#   aₖ = (c/8)(25k² − 50k + 24) a_{k−1} − (1/2) Σ_{i=1}^{k−1} aᵢ a_{k−i},   a₀ = 1.
# At c = 1/6 (GIKM) this is a_{n+1} = ((25n²−1)/48)aₙ − ½Σ, so a₁ = −1/48, a₂ = −49/4608;
# at c = 1/12 (:string) a₁ = −1/96. Every aₙ is an exact rational.
function _painleve1(n::Integer, c::Rational{Int} = 1 // 6)
    a = Vector{Rational{BigInt}}(undef, n)
    a[1] = 1
    for k in 1:(n - 1)                       # fill a_k = a[k+1]
        conv = sum((a[i + 1] * a[k - i + 1] for i in 1:(k - 1)); init = zero(Rational{BigInt}))
        a[k + 1] = (c / 8) * (25 * k^2 - 50 * k + 24) * a[k] - conv // 2
    end
    a
end

# Bender–Wu recursion for H = p²/2 + x²/2 + g x⁴, ground state E₀(g) = Σ_m E_m g^m.
# With ψ = e^{-x²/2} Σ_m g^m P_m(x), P_0 = 1, P_m = Σ_{j=1}^{2m} c_{m,j} x^{2j},
# the Schrödinger equation gives, per power x^{2j},
#   2j c_{m,j} = (j+1)(2j+1) c_{m,j+1} − c_{m-1,j-2} + Σ_{k=1}^{m-1} E_k c_{m-k,j},
# solved downward from c_{m,2m} = −c_{m-1,2m-2}/(4m); the j = 0 row gives
# E_m = −c_{m,1}.
function _bender_wu(n::Integer)
    E = Vector{Rational{BigInt}}(undef, n)
    E[1] = 1//2
    c = Vector{Vector{Rational{BigInt}}}(undef, n)
    c[1] = Rational{BigInt}[]
    # c_{m,j} with the conventions c_{m,0} = δ_{m,0} and c_{m,j} = 0 outside 1 ≤ j ≤ 2m
    at(m, j) = j == 0 ? Rational{BigInt}(m == 0 ? 1 : 0) :
               (1 ≤ j ≤ 2m ? c[m + 1][j] : zero(Rational{BigInt}))
    for m in 1:(n - 1)
        row = zeros(Rational{BigInt}, 2m)
        row[2m] = -at(m - 1, 2m - 2) // (4m)
        for j in (2m - 1):-1:1
            acc = (j + 1) * (2j + 1) * (j + 1 ≤ 2m ? row[j + 1] : zero(Rational{BigInt}))
            acc -= at(m - 1, j - 2)
            for k in 1:(m - 1)
                acc += E[k + 1] * at(m - k, j)
            end
            row[j] = acc // (2j)
        end
        c[m + 1] = row
        E[m + 1] = -row[1]
    end
    E
end

"""
    Transseries(name::Symbol, n::Integer; var = :ħ) -> Transseries{Rational{BigInt}}

Named course-canon transseries, each with `n` coefficients per sector:

- `:euler` — action ``A = 1``, sectors ``Φ_0 = \\sum_k k!\\,ħ^{k+1}`` (the
  ``ħ → -ħ`` mirror of `FormalSeries(:euler)`, placing the Borel pole at
  ``ζ = +1 = A`` on the Stokes ray ``θ = 0``) and ``Φ_1 = 1``; all higher sectors
  vanish identically. Stokes constant ``S_1 = -2πi`` (see
  [`stokes_constant`](@ref)).
- `:airy` — action ``A = -4/3`` (the Borel singularity of the Ai sector; the
  relative weight ``e^{-A/ħ} = e^{4/(3ħ)} = e^{2ξ}`` is the Bi/Ai ratio), sectors
  ``Φ_0 = `` `FormalSeries(:airy, n)` and ``Φ_1 = `` `FormalSeries(:airy_bi, n)`.
"""
function Transseries(name::Symbol, n::Integer; var::Symbol = :ħ)
    n ≥ 1 || throw(InvalidArgument("need n ≥ 1 coefficients, got $n"))
    if name === :euler
        c = Vector{Rational{BigInt}}(undef, n)
        c[1] = 1
        for k in 1:(n - 1)
            c[k + 1] = k * c[k]
        end
        Φ0 = FormalSeries(c, var; power_offset = 1//1)
        Φ1 = FormalSeries([one(Rational{BigInt})], var)
        Transseries(1//1, [Φ0, Φ1], var)
    elseif name === :airy
        Transseries(-4//3, [FormalSeries(:airy, n; var), FormalSeries(:airy_bi, n; var)],
                    var)
    else
        throw(InvalidArgument("unknown named transseries :$name (expected :euler or :airy)"))
    end
end
