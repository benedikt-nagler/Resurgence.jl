# Text and LaTeX display for the three public types. Shows the first few terms plus
# a truncation-order tag; never dumps hundreds of coefficients.

const _SHOW_TERMS = 5

function _term_string(c, var::Symbol, p::Rational{Int}; latex::Bool = false)
    pow = if iszero(p)
        ""
    elseif isone(p)
        string(var)
    elseif denominator(p) == 1
        latex ? "$var^{$(numerator(p))}" : "$var^$(numerator(p))"
    else
        latex ? "$var^{$(numerator(p))/$(denominator(p))}" : "$var^($p)"
    end
    cs = latex && c isa Rational ?
         (denominator(c) == 1 ? string(numerator(c)) :
          "\\frac{$(numerator(c))}{$(denominator(c))}") : string(c)
    sep = latex ? "\\," : "*"
    isempty(pow) ? cs : (isone(c) ? pow : (c == -1 ? "-$pow" : "$cs$sep$pow"))
end

function _series_string(coeffs, var::Symbol, offset::Rational{Int}; latex::Bool = false)
    parts = String[]
    shown = 0
    for (k, c) in pairs(coeffs)
        shown ≥ _SHOW_TERMS && break
        iszero(c) && continue
        push!(parts, _term_string(c, var, offset + (k - 1); latex))
        shown += 1
    end
    body = isempty(parts) ? "0" : join(parts, " + ")
    body = replace(body, "+ -" => "- ")
    tail_pow = offset + length(coeffs)
    order = latex ? "O($var^{$(numerator(tail_pow))" *
                    (denominator(tail_pow) == 1 ? "}" : "/$(denominator(tail_pow))}") * ")" :
                    "O($var^" * (denominator(tail_pow) == 1 ?
                                 string(numerator(tail_pow)) : "($tail_pow)") * ")"
    "$body + $order"
end

function Base.show(io::IO, Φ::FormalSeries{T}) where {T}
    print(io, "FormalSeries{$T}: ",
          _series_string(Φ.coeffs, Φ.var, Φ.power_offset))
end

function Base.show(io::IO, ::MIME"text/latex", Φ::FormalSeries)
    print(io, "\$", _series_string(Φ.coeffs, Φ.var, Φ.power_offset; latex = true), "\$")
end

function Base.show(io::IO, B::BorelSeries{T}) where {T}
    print(io, "BorelSeries{$T} (β = $(B.beta), from $(B.source_var)): ")
    iszero(B.constant_term) || print(io, "[const $(B.constant_term)] + ")
    print(io, "1/Γ(β) * (",
          _series_string(B.series.coeffs, variable(B), B.series.power_offset), ")")
end

function Base.show(io::IO, ::MIME"text/latex", B::BorelSeries)
    print(io, "\$\\frac{1}{\\Gamma($(B.beta))}\\left(",
          _series_string(B.series.coeffs, variable(B), B.series.power_offset;
                         latex = true), "\\right)\$")
end

function Base.show(io::IO, r::PadeApproximant{T}) where {T}
    print(io, "PadeApproximant{$T}: [$(numerator_degree(r))/",
          "$(denominator_degree(r))] in $(r.var)")
end

function Base.show(io::IO, m::ConformalMap{T}) where {T}
    print(io, "ConformalMap{$T}: cut plane → unit disk, branch point ζ₀ = ",
          m.zeta0, m.pair ? " (symmetric pair ±ζ₀)" : "")
end

function Base.show(io::IO, c::ConformalPade{T}) where {T}
    print(io, "ConformalPade{$T}: [$(numerator_degree(c.inner))/",
          "$(denominator_degree(c.inner))] in w(", c.var, "), ζ₀ = ", c.map.zeta0,
          c.map.pair ? " (pair)" : "")
end

function Base.show(io::IO, F::Transseries{T}) where {T}
    print(io, "Transseries{$T} (A = $(F.action), $(n_sectors(F)) sectors in $(F.var)): ",
          "Φ₀ = ", _series_string(F.sectors[1].coeffs, F.var, F.sectors[1].power_offset))
end

function Base.show(io::IO, ::MIME"text/latex", F::Transseries)
    A = F.action
    As = A isa Rational ?
         (denominator(A) == 1 ? string(numerator(A)) :
          "\\frac{$(numerator(A))}{$(denominator(A))}") : string(A)
    print(io, "\$\\sum_{n=0}^{$(n_sectors(F) - 1)} \\sigma^n ",
          "e^{-n \\left($As\\right)/$(F.var)}\\,\\Phi_n, \\quad \\Phi_0 = ",
          _series_string(F.sectors[1].coeffs, F.var, F.sectors[1].power_offset;
                         latex = true), "\$")
end
