# Large-order analysis: Richardson extrapolation and the factorial-growth fit
# a_n ~ S · Γ(n+b) / A^n that reads off the leading Borel singularity (|A| = its
# distance) from the coefficients alone.

"""
    richardson(seq::AbstractVector, order::Integer) -> Number

Richardson extrapolation of a sequence ``s_n = s_∞ + c_1/n + c_2/n^2 + …`` (the
element `seq[n]` is taken as ``s_n``): performs `order` elimination sweeps
``s_n ← ((n+k)\\,s_{n+1} - n\\,s_n)/k`` and returns the last surviving element,
accurate to ``O(n^{-(order+1)})``.
"""
function richardson(seq::AbstractVector, order::Integer)
    N = length(seq)
    0 ≤ order < N ||
        throw(InvalidArgument("richardson order must be in 0:$(N - 1), got $order"))
    s = float.(collect(seq))
    for k in 1:order
        for n in 1:(length(s) - 1)
            s[n] = ((n + k) * s[n + 1] - n * s[n]) / k
        end
        pop!(s)
    end
    s[end]
end

"""
    large_order_fit(Φ::FormalSeries; order = 4, method = :richardson)
        -> NamedTuple{(:A, :b, :S)}
    large_order_fit(Φ::FormalSeries; saddles = J, order = nothing,
                    window = nothing, maxiter = 50)
        -> NamedTuple{(:saddles, :converged)}

Fit the factorial growth ``|a_n| \\sim S\\, Γ(p_n + b)\\, A^{-p_n}`` of the
coefficients, where ``p_n = n + \\mathrm{power\\_offset}`` is the actual ħ-power of
each term. Returns `(A = …, b = …, S = …)`; ``A`` is the modulus of the action —
the distance of the leading Borel-plane singularity — extracted by
Richardson-extrapolating (`order` sweeps) the ratio sequences

- ``ρ_{n+1} - ρ_n → 1/A`` with ``ρ_n = |a_{n+1}/a_n|``,
- ``A ρ_n - p_n → b``,
- ``|a_n|\\, A^{p_n} / Γ(p_n + b) → S``.

Oracles: the Euler series gives ``(1, 0, 1)`` exactly; the Airy series gives
``A = 4/3``. Precision follows the coefficient type (`Rational{BigInt}` coefficients
are fitted in `BigFloat`). `method` selects the extrapolation used on the three
sequences — any [`accelerate`](@ref) method (`order` is passed through as its depth).

With `saddles = J > 1` the multi-saddle model
``a_n ≈ \\sum_{j=1}^{J} S_j\\, Γ(p_n + b_j)\\, A_j^{-p_n}`` is fitted instead
(complex actions ``A_j`` — the Borel singularity positions themselves, which also
initialize the fit; real ``b_j``): amplitudes by linear least squares over the
fit `window` (0-based coefficient indices; default the last half), actions and
exponents by damped Gauss–Newton (`maxiter` iterations; here `order` is passed to
the initializing [`pade`](@ref)). For real coefficients, complex-conjugate action
pairs are detected and fitted as one oscillatory group
``2\\,\\Re[S\\,Γ(p_n+b)\\,A^{-p_n}]`` — the equal-modulus ``A, \\bar A`` case where
the single-saddle ratio fit fails. Returns
`(saddles = [(A, b, S), …] sorted by |A|, converged::Bool)` without throwing on
non-convergence.
"""
function large_order_fit(Φ::FormalSeries; order = nothing,
                         method::Symbol = :richardson, saddles::Integer = 1,
                         window::Union{Nothing,AbstractRange} = nothing,
                         maxiter::Integer = 50)
    saddles ≥ 1 || throw(InvalidArgument("saddles must be ≥ 1, got $saddles"))
    if saddles > 1
        return _multi_saddle_fit(Φ, Int(saddles); order, window, maxiter)
    end
    order = order === nothing ? 4 : Int(order)
    N = n_terms(Φ)
    N ≥ order + 3 ||
        throw(InvalidArgument("large_order_fit needs ≥ $(order + 3) terms, got $N"))
    m = [float(abs(c)) for c in Φ.coeffs]
    any(iszero, m) &&
        throw(InvalidArgument("large_order_fit requires nonzero coefficients"))
    T = typeof(m[1])
    p0 = T(power_offset(Φ))
    ρ = [m[n + 1] / m[n] for n in 1:(N - 1)]                 # ρ_n ≈ (p_n + b)/A
    dρ = [ρ[n + 1] - ρ[n] for n in 1:(N - 2)]                # → 1/A
    A = 1 / accelerate(dρ; method, order)
    bseq = [A * ρ[n] - ((n - 1) + p0) for n in 1:(N - 1)]    # → b  (p_n, 0-based n)
    b = accelerate(bseq; method, order)
    # S via one Γ evaluation plus upward recurrence Γ(x+1) = xΓ(x)
    x1 = p0 + b                                              # p_0 + b
    g = _gamma(x1 > 0 ? x1 : x1 + ceil(-x1) + 1)
    if x1 ≤ 0                                                # undo the shift
        z = x1 + ceil(-x1) + 1
        while z > x1 + 1//2
            z -= 1
            g /= z
        end
    end
    Sseq = Vector{T}(undef, N)
    Apow = A^p0
    for n in 0:(N - 1)
        Sseq[n + 1] = m[n + 1] * Apow / g
        g *= (n + p0 + b)                                    # Γ(p_{n+1} + b)
        Apow *= A
    end
    S = accelerate(Sseq; method, order)
    (A = A, b = b, S = S)
end
