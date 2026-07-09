# Sequence-acceleration toolbox: Shanks, Wynn-ε, and Levin transforms behind one
# common interface (`accelerate`), alongside the M1 `richardson`. Used both to sum
# (accelerated) sequences of partial sums and to sharpen `large_order_fit`.

"""
    accelerate(seq::AbstractVector; method = :wynn, order = nothing) -> Number

Accelerate the convergence of the sequence `seq` (element `seq[n]` is ``s_n``) and
return the extrapolated limit. Available `method`s:

- `:richardson` - delegates to [`richardson`](@ref) (`order` sweeps, default
  `min(4, length - 1)`); best for sequences with ``1/n`` power-law error.
- `:shanks` - `order` iterated Shanks transformations
  ``S(s)_n = (s_{n+1} s_{n-1} - s_n^2)/(s_{n+1} + s_{n-1} - 2 s_n)`` (default: as
  many as the length allows); exact on ``s_n = s + α q^n``.
- `:wynn` - Wynn's ε-algorithm, descending to even column ``ε_{2\\,order}``
  (default: the deepest reachable even column); an efficient equivalent of iterated
  Shanks, robust for alternating and divergent-series partial sums.
- `:levin_u`, `:levin_t` - Levin transforms with ``β = 1`` and remainder estimates
  ``ω_n = (n+1)\\,Δs_{n-1}`` (u) or ``ω_n = Δs_{n-1}`` (t), using `order + 2`
  sequence elements (default: all); typically the fastest accelerator for both
  convergent and alternating divergent series.

Precision follows the element type (`Rational` input is promoted to `BigFloat`);
complex sequences work. See [`partial_sums`](@ref) for building the input from a
[`FormalSeries`](@ref).
"""
function accelerate(seq::AbstractVector; method::Symbol = :wynn,
                    order::Union{Nothing,Integer} = nothing)
    N = length(seq)
    N ≥ 2 || throw(InvalidArgument("accelerate needs ≥ 2 sequence elements, got $N"))
    order === nothing || order ≥ 0 ||
        throw(InvalidArgument("acceleration order must be ≥ 0, got $order"))
    s = float.(collect(seq))
    if method === :richardson
        richardson(s, order === nothing ? min(4, N - 1) : Int(order))
    elseif method === :shanks
        k = order === nothing ? (N - 1) ÷ 2 : Int(order)
        2k + 1 ≤ N ||
            throw(InvalidArgument("shanks order $k needs ≥ $(2k + 1) elements, got $N"))
        _shanks(s, k)
    elseif method === :wynn
        kmax = (N - 1) ÷ 2
        k = order === nothing ? kmax : Int(order)
        k ≤ kmax ||
            throw(InvalidArgument("wynn order $k needs ≥ $(2k + 1) elements, got $N"))
        _wynn_epsilon(s, k)
    elseif method === :levin_u || method === :levin_t
        k = order === nothing ? N - 2 : Int(order)
        k ≤ N - 2 ||
            throw(InvalidArgument("levin order $k needs ≥ $(k + 2) elements, got $N"))
        _levin(s, k, method === :levin_u)
    else
        throw(InvalidArgument("unknown acceleration method :$method (expected " *
                              ":richardson, :shanks, :wynn, :levin_u, or :levin_t)"))
    end
end

# Iterated Shanks transformation; a zero denominator means the local triple already
# converged exactly, so the exact limit s_{n+1} is passed through.
function _shanks(s::Vector{T}, order::Integer) where {T}
    for _ in 1:order
        t = Vector{T}(undef, length(s) - 2)
        for n in 2:(length(s) - 1)
            d = s[n + 1] + s[n - 1] - 2 * s[n]
            t[n - 1] = iszero(d) ? s[n + 1] : (s[n + 1] * s[n - 1] - s[n]^2) / d
        end
        s = t
    end
    s[end]
end

# Wynn ε-algorithm: ε_{-1}^{(n)} = 0, ε_0^{(n)} = s_n,
# ε_{k+1}^{(n)} = ε_{k-1}^{(n+1)} + 1/(ε_k^{(n+1)} - ε_k^{(n)}). Only even columns
# estimate the limit (ε_2^{(n)} is the Shanks transform); descend to column 2·kmax,
# stopping early (last finite even-column value) on a zero denominator.
function _wynn_epsilon(s::Vector{T}, kmax::Integer) where {T}
    prev = zeros(T, length(s))
    curr = s
    best = curr[end]
    col = 0
    while col < 2 * kmax
        next = Vector{T}(undef, length(curr) - 1)
        stop = false
        for i in 1:(length(curr) - 1)
            d = curr[i + 1] - curr[i]
            iszero(d) && (stop = true; break)
            next[i] = prev[i + 1] + 1 / d
        end
        stop && break
        prev, curr = curr, next
        col += 1
        iseven(col) && (best = curr[end])
    end
    best
end

# Levin transform with β = 1, base index n₀ = 1 (0-based), depth k:
# L_k = Σ_j (-1)^j C(k,j) c_{jk} s_{n₀+j}/ω_{n₀+j} / Σ_j (-1)^j C(k,j) c_{jk}/ω_{n₀+j}
# with c_{jk} = ((1+n₀+j)/(1+n₀+k))^{k-1} and ω_n = Δs_{n-1} (t) or (n+1)Δs_{n-1} (u).
# A zero difference means the sequence already converged exactly.
function _levin(s::Vector{T}, k::Integer, u::Bool) where {T}
    R = real(T)
    num = zero(T)
    den = zero(T)
    C = one(R)                              # binomial C(k, j), built iteratively
    for j in 0:k
        n = 1 + j                           # 0-based sequence index used
        a = s[n + 1] - s[n]                 # Δs_{n-1} = s_n - s_{n-1}
        iszero(a) && return s[end]
        ω = u ? (n + 1) * a : a
        c = C * (R(2 + j) / R(2 + k))^(k - 1)
        t = (isodd(j) ? -c : c) / ω
        num += t * s[n + 1]
        den += t
        C *= R(k - j) / R(j + 1)
    end
    num / den
end

"""
    partial_sums(Φ::FormalSeries, x) -> Vector

The partial sums ``s_k = \\sum_{n=0}^{k} a_n x^{p_n}`` of the series evaluated at
`x`, where ``p_n = n + \\mathrm{power\\_offset}`` - the natural input for
[`accelerate`](@ref), both for convergent series and for divergent ones (where the
Levin/Wynn transforms act as summation methods).
"""
function partial_sums(Φ::FormalSeries, x)
    xn = iszero(Φ.power_offset) ? one(x) : x^Φ.power_offset
    acc = zero(Φ.coeffs[1] * xn)
    out = Vector{typeof(acc)}(undef, n_terms(Φ))
    for n in 1:n_terms(Φ)
        acc += Φ.coeffs[n] * xn
        out[n] = acc
        xn *= x
    end
    out
end
