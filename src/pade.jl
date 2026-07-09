# Padé approximants over any field: one generic Toeplitz solve, pivoting chosen by an
# exactness trait (partial pivoting for numeric types, first-nonzero for exact rings).
# Polynomials.jl's Pade has a known correctness bug and other packages are
# numeric-only, hence self-owned (~50 lines of core).

"""
    AbstractBorelApproximant

Abstract supertype of approximants of the *reduced* Borel function ``Σ b_n ζ^n``
(the ``ζ^{β-1}/Γ(β)`` factor stays with [`laplace_sum`](@ref)). The interface a
subtype must provide:

- callable at complex ``ζ``: `a(ζ)` evaluates the approximant,
- `poles(a; refine = true)` returning the physical-sheet singularities as a
  `Vector{Complex{BigFloat}}` sorted by absolute value.

Any such approximant plugs into the Laplace pipeline via
`laplace_sum(B, a, ħ)` / `lateral_sum(B, a, ħ)`. Subtypes:
[`PadeApproximant`](@ref), [`ConformalPade`](@ref), and the `AAAApproximant` of
the BaryRational package extension.
"""
abstract type AbstractBorelApproximant end

"""
    PadeApproximant{T}

The rational function ``[L/M](z) = p(z)/q(z)`` with ``\\deg p = L``, ``\\deg q = M``,
``q(0) = 1``, matching a power series through order ``z^{L+M}``. Callable:
`r(z)` evaluates ``p(z)/q(z)`` by Horner. Fields `p`, `q` hold ascending
coefficients; `var` names the variable.
"""
struct PadeApproximant{T} <: AbstractBorelApproximant
    p::Vector{T}
    q::Vector{T}
    var::Symbol
end

numerator_degree(r::PadeApproximant) = length(r.p) - 1
denominator_degree(r::PadeApproximant) = length(r.q) - 1
variable(r::PadeApproximant) = r.var

_horner(c::Vector, z) = foldr((cᵢ, acc) -> cᵢ + z * acc, c; init = zero(c[1]) * zero(z))

(r::PadeApproximant)(z) = _horner(r.p, z) / _horner(r.q, z)

# Dense Gaussian elimination, in place; pivoting by the exactness trait — max-|·|
# for numeric types, first nonzero for exact rings (where `abs`/`one(::Type)` need
# not exist, so LinearAlgebra's generic lu is unusable). Returns nothing on a
# singular matrix.
function _solve_dense!(A::Matrix{T}, b::Vector{T}) where {T}
    n = length(b)
    exact = _is_exact_type(T)
    for k in 1:n
        piv = 0
        if exact
            for i in k:n
                iszero(A[i, k]) && continue
                piv = i
                break
            end
        else
            piv = k
            for i in (k + 1):n
                abs(A[i, k]) > abs(A[piv, k]) && (piv = i)
            end
            iszero(A[piv, k]) && (piv = 0)
        end
        piv == 0 && return nothing
        if piv != k
            for j in k:n
                A[k, j], A[piv, j] = A[piv, j], A[k, j]
            end
            b[k], b[piv] = b[piv], b[k]
        end
        for i in (k + 1):n
            f = A[i, k] / A[k, k]
            for j in k:n
                A[i, j] -= f * A[k, j]
            end
            b[i] -= f * b[k]
        end
    end
    x = b
    for k in n:-1:1
        acc = b[k]
        for j in (k + 1):n
            acc -= A[k, j] * x[j]
        end
        x[k] = acc / A[k, k]
    end
    x
end

"""
    pade(c::AbstractVector, L, M; var = :ζ, reduce = false) -> PadeApproximant
    pade(Φ::FormalSeries, L, M; reduce = false) -> PadeApproximant

The ``[L/M]`` Padé approximant of the series with coefficients `c` (needs
`L + M + 1` of them; on a `FormalSeries` the polynomial part ``Σ a_n x^n`` is used,
ignoring the power offset). Works over any field; throws
[`DegeneratePade`](@ref) when the Toeplitz system is singular. With
`reduce = true` a singular system is instead retried at ``[L-1/M-1]`` (walking
down the degeneracy block, by Baker's block theorem, until the system is regular
or ``M = 0``) — the reduced approximant still matches the series through the
originally requested order when the degeneracy is exact.
"""
function pade(c::AbstractVector{T}, L::Integer, M::Integer; var::Symbol = :ζ,
              reduce::Bool = false) where {T}
    (L ≥ 0 && M ≥ 0) || throw(InvalidArgument("Padé orders must be ≥ 0, got [$L/$M]"))
    length(c) ≥ L + M + 1 ||
        throw(InvalidArgument("[$L/$M] Padé needs $(L + M + 1) coefficients, got $(length(c))"))
    at(m) = m < 0 ? zero(c[1]) : c[m + 1]
    if M == 0
        return PadeApproximant(collect(c[1:(L + 1)]), [one(c[1])], var)
    end
    # q₀ = 1; solve Σ_{j=1}^{M} q_j c_{L+i-j} = -c_{L+i},  i = 1..M
    A = [at(L + i - j) for i in 1:M, j in 1:M]
    rhs = [-at(L + i) for i in 1:M]
    qtail = _solve_dense!(A, rhs)
    if qtail === nothing
        reduce || throw(DegeneratePade(Int(L), Int(M)))
        return pade(c, max(L - 1, 0), M - 1; var, reduce = true)
    end
    q = vcat(one(c[1]), qtail)
    p = [sum(q[j + 1] * at(k - j) for j in 0:min(k, M)) for k in 0:L]
    PadeApproximant(p, q, var)
end

pade(Φ::FormalSeries, L::Integer, M::Integer; reduce::Bool = false) =
    pade(Φ.coeffs, L, M; var = Φ.var, reduce)

"""
    pade(B::BorelSeries; order = nothing, reduce = false) -> PadeApproximant

Padé approximant of the *reduced* Borel function ``Σ b_n ζ^n`` (the stored
coefficients; the overall ``ζ^{β-1}/Γ(β)`` factor does not move the poles and is
re-applied by [`laplace_sum`](@ref)). With `order = m` builds the diagonal
``[m/m]`` approximant (needs ``2m+1`` terms); by default uses all available terms,
near-diagonally: ``M = ⌊N/2⌋``, ``L = N - M`` for ``N = n\\_terms - 1``.
`reduce = true` degrades a degenerate request gracefully (see [`pade`](@ref)).
"""
function pade(B::BorelSeries; order::Union{Nothing,Integer} = nothing,
              reduce::Bool = false)
    N = n_terms(B) - 1
    if order === nothing
        M = N ÷ 2
        L = N - M
    else
        L = M = Int(order)
    end
    pade(B.series.coeffs, L, M; var = variable(B), reduce)
end

"""
    poles(r::PadeApproximant; refine = true) -> Vector{Complex}

Zeros of the denominator ``q``: located via the companion-matrix eigenvalues in
`ComplexF64`, then (with `refine = true`) polished by Newton iteration in
`Complex{BigFloat}`. Sorted by absolute value — for a Borel–Padé approximant the
smallest is the leading Borel singularity.
"""
function poles(r::PadeApproximant; refine::Bool = true)
    M = denominator_degree(r)
    M == 0 && return Complex{BigFloat}[]
    q64 = ComplexF64.(Complex{BigFloat}.(r.q))
    # companion matrix of q normalized to monic
    C = zeros(ComplexF64, M, M)
    for i in 2:M
        C[i, i - 1] = 1
    end
    C[:, M] .= -q64[1:M] ./ q64[M + 1]
    ζs = LinearAlgebra.eigvals(C)
    out = Complex{BigFloat}[Complex{BigFloat}(z) for z in ζs]
    if refine
        qb = Complex{BigFloat}.(r.q)
        dqb = [k * qb[k + 1] for k in 1:M]
        for (i, z) in pairs(out)
            for _ in 1:8
                dq = _horner(dqb, z)
                iszero(dq) && break
                step = _horner(qb, z) / dq
                z -= step
                abs(step) ≤ eps(BigFloat) * (1 + abs(z)) && break
            end
            out[i] = z
        end
    end
    sort!(out; by = abs)
end

"""
    residues(r::PadeApproximant; refine = true) -> Vector{Pair}

Pairs `ζ₀ => p(ζ₀)/q′(ζ₀)` for each simple pole of `r`, in the order of
[`poles`](@ref). This is the data alien derivatives consume (M2).
"""
function residues(r::PadeApproximant; refine::Bool = true)
    ζs = poles(r; refine)
    pb = Complex{BigFloat}.(r.p)
    qb = Complex{BigFloat}.(r.q)
    dqb = [k * qb[k + 1] for k in 1:denominator_degree(r)]
    [ζ => _horner(pb, ζ) / _horner(dqb, ζ) for ζ in ζs]
end

"""
    borel_pade_poles(Φ::FormalSeries; beta = power_offset(Φ), order = nothing,
                     refine = true, reduce = false) -> Vector{Complex}

The one-liner: Borel transform, Padé approximation, pole extraction — the Padé
estimate of the Borel-plane singularities of `Φ`.
"""
function borel_pade_poles(Φ::FormalSeries;
                          beta::Union{Rational{Int},Integer} = power_offset(Φ),
                          order::Union{Nothing,Integer} = nothing,
                          refine::Bool = true, reduce::Bool = false)
    poles(pade(borel(Φ; beta); order, reduce); refine)
end
