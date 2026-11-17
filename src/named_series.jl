# Exact coefficient reference cases, not general model solvers.
"""
    FormalSeries(name::Symbol, n; var = :ħ)

Construct `n` exact rational coefficients of a retained reference series.
Names: `:euler`, `:airy`, `:airy_bi`, `:quartic`, `:phi4`, `:exp`.
Euler uses power offset 1 and Borel function `-1/(1+zeta)`.
Airy uses hbar = z^(-3/2). Quartic uses H = p^2/2 + x^2/2 + g*x^4.
Phi4 is the formal expansion of the normalized integral of exp(-x^2/2-g*x^4).
These constructors supply coefficients, not transseries or physical Stokes data.
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
        # :airy - coefficients (−1)^k u_k (3/2)^k: c_{k+1} = −c_k (6k+5)(6k+1)/(48(k+1));
        # :airy_bi - the same without the sign.
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
    else
        throw(InvalidArgument("unknown named series :$name (expected :euler, :airy, " *
                              ":airy_bi, :quartic, :phi4, :exp)"))
    end
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
