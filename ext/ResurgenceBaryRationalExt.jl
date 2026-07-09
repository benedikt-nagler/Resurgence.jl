# AAA rational approximation of the reduced Borel function (BaryRational.jl),
# behind a package extension: an alternative pole/branch-cut hunter that plugs into
# the Laplace pipeline through the AbstractBorelApproximant seam.
module ResurgenceBaryRationalExt

using BaryRational: BaryRational, aaa, prz
using Resurgence: Resurgence, BorelSeries, AbstractBorelApproximant,
                  NoSingularityFound, InvalidArgument, pade, poles

"""
    AAAApproximant{A} <: AbstractBorelApproximant

An AAA barycentric rational approximant of the reduced Borel function (wrapping
`BaryRational.AAAapprox`). Callable at complex ζ; `poles`/`residues` via
`BaryRational.prz`. Built by `aaa_approximant` (from samples) or `aaa_borel`
(from a `BorelSeries`).
"""
struct AAAApproximant{A} <: AbstractBorelApproximant
    inner::A
    var::Symbol
end

(a::AAAApproximant)(ζ) = a.inner(ζ)

Resurgence.variable(a::AAAApproximant) = a.var

function Resurgence.aaa_approximant(zs::AbstractVector, fs::AbstractVector;
                                    tol::Real = 1e-13, mmax::Integer = 100,
                                    var::Symbol = :ζ)
    length(zs) == length(fs) ||
        throw(InvalidArgument("aaa_approximant needs equally many sample points " *
                              "and values, got $(length(zs)) and $(length(fs))"))
    AAAApproximant(aaa(collect(zs), collect(fs); tol = Float64(tol),
                       mmax = Int(mmax)), var)
end

function Resurgence.aaa_borel(B::BorelSeries; samples::Integer = 128,
                              radius::Union{Nothing,Real} = nothing,
                              tol::Real = 1e-13, mmax::Integer = 100)
    b = B.series.coeffs
    N = length(b)
    ρ = radius
    if ρ === nothing
        ζs = poles(pade(B; reduce = true))
        R = if !isempty(ζs)
            Float64(abs(ζs[1]))
        else
            # crude radius-of-convergence estimate from the last coefficient ratio
            r = abs(Float64(float(real(b[N - 1]))) + im * Float64(float(imag(b[N - 1])))) /
                abs(Float64(float(real(b[N]))) + im * Float64(float(imag(b[N]))))
            (isfinite(r) && r > 0) ||
                throw(NoSingularityFound("aaa_borel could not estimate a sampling " *
                                         "radius; pass radius explicitly"))
            r
        end
        # keep the series-truncation error at the sample points below tol:
        # (ρ/R)^N ≤ tol
        ρ = R * min(Float64(tol)^(1 / N), 0.9)
    end
    # sample the truncated reduced Borel sum on two concentric circles
    c64 = [ComplexF64(ComplexF64(float(real(x)), float(imag(x)))) for x in b]
    horner(z) = foldr((cᵢ, acc) -> cᵢ + z * acc, c64; init = zero(z))
    m = max(Int(samples) ÷ 2, 8)
    zs = ComplexF64[]
    for r in (ρ, ρ / 2)
        append!(zs, [r * cispi(2k / m) for k in 0:(m - 1)])
    end
    AAAApproximant(aaa(zs, horner.(zs); tol = Float64(tol), mmax = Int(mmax)),
                   Resurgence.variable(B))
end

function Resurgence.poles(a::AAAApproximant; refine::Bool = true)
    pol, _, _ = prz(a.inner)
    sort!(Complex{BigFloat}.(pol); by = abs)
end

function Resurgence.residues(a::AAAApproximant; refine::Bool = true)
    pol, res, _ = prz(a.inner)
    perm = sortperm(pol; by = abs)
    [Complex{BigFloat}(pol[i]) => Complex{BigFloat}(res[i]) for i in perm]
end

end # module ResurgenceBaryRationalExt
