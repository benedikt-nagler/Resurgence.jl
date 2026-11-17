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
    large_order_fit(Φ::FormalSeries; order = nothing, method = :richardson)

Estimate `(A, b, S)` from one dominant factorial growth law of nonzero coefficients.
`A` is the action modulus, not a complex singularity location.
Assumes a single asymptotic law; oscillatory and sparse input is unsupported.
Extrapolation is diagnostic. A small residual does not certify the fitted model.
"""
function large_order_fit(Φ::FormalSeries; order = nothing, method::Symbol = :richardson)
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
    isfinite(A) && A > 0 || throw(InvalidArgument(
        "single-law action estimate is not finite and positive; the model is unresolved"))
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
