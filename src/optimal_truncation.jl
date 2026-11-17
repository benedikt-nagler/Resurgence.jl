# Stored coefficients beyond the truncation are unknown, not zero.
function _term_values(Φ::FormalSeries, x)
    xn = iszero(power_offset(Φ)) ? one(x) : x^power_offset(Φ)
    out = [c*xn*x^(i-1) for (i,c) in enumerate(coefficients(Φ))]
    all(isfinite,out) || throw(InvalidArgument("term magnitudes overflowed; use higher precision"))
    for (c,v) in zip(coefficients(Φ),out)
        !iszero(c) && !iszero(x) && iszero(v) && throw(InvalidArgument(
            "a nonzero term underflowed; use higher precision"))
    end
    out
end

"""
    optimal_truncation(Φ, ħ; rule = :smallest_term, A = nothing)

Select a truncation from nonzero stored terms and return `(N, value, error)`.
`N` counts retained coefficient slots, including structural zeros.
`error` is the first omitted nonzero term's magnitude, not a remainder bound.
If the stored terms are exhausted without a resolved minimum, keep all terms
and return `error = nothing`. The same applies to an all-zero stored series.
`:action` uses a supplied action modulus or a single-law large-order estimate.
Throws [`InvalidArgument`](@ref) on numerical overflow or underflow.
"""
function optimal_truncation(Φ::FormalSeries, ħ::Number; rule::Symbol = :smallest_term,
                            A::Union{Nothing,Number} = nothing)
    terms = _term_values(Φ, ħ)
    indices = findall(!iszero, coefficients(Φ))
    N = n_terms(Φ)
    omitted = nothing
    if rule === :smallest_term
        if !isempty(indices)
            k = argmin(abs.(terms[indices]))
            # Without a later nonzero term, the stored range does not resolve a turn.
            k < length(indices) && (omitted = indices[k])
        end
    elseif rule === :action
        iszero(ħ) && throw(InvalidArgument("action truncation requires nonzero hbar"))
        Amod = abs(A === nothing ? large_order_fit(Φ).A : A)
        isfinite(Amod) && Amod > 0 || throw(InvalidArgument("action modulus must be finite and positive"))
        estimated = Amod/abs(ħ)-float(power_offset(Φ))
        target = estimated >= N ? N : estimated <= 0 ? 0 : round(Int,estimated)
        k = findfirst(i -> i > target, indices)
        k === nothing || (omitted = indices[k])
    else
        throw(InvalidArgument("rule must be :smallest_term or :action"))
    end
    kept = omitted === nothing ? N : omitted-1
    value = sum(terms[1:kept];init=zero(first(terms)))
    (N=kept,value=value,error=omitted===nothing ? nothing : abs(terms[omitted]))
end
