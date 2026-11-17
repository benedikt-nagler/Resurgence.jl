# Borel–Padé–Laplace summation: QuadGK ray integrals through the Padé approximant of
# the reduced Borel function. The 1/Γ(β) dropped by the rising-factorial
# normalization (borel.jl) is re-applied here, numerically, in the caller's
# precision. Precision flows from the argument types - no global configuration.

# Γ(x) for real x > 0 in the precision of x: exact for integers, else reduced to
# [1, 2) by the recurrence Γ(x+1) = xΓ(x) and integrated there with QuadGK.
# Keeps SpecialFunctions out of [deps]; the tests cross-check against it.
function _gamma(x::Real)
    x > 0 || throw(InvalidArgument("_gamma requires x > 0, got $x"))
    if isinteger(x)
        return oftype(float(x), factorial(big(Int(x) - 1)))
    end
    T = typeof(float(x))
    frac = float(x) - floor(float(x))          # in (0, 1)
    y = frac + 1                                # in (1, 2): integrand is smooth
    g, _ = quadgk(t -> t^(y - 1) * exp(-t), zero(T), T(Inf);
                  rtol = eps(T)^(4//5))
    # upward recurrence Γ(x) = (x-1)(x-2)⋯(y)·Γ(y), downward for x < 1
    z = y
    while z < float(x) - 1//2
        g *= z
        z += 1
    end
    while z > float(x) + 1//2
        z -= 1
        g /= z
    end
    g
end

_real_float_type(ħ::Number) = typeof(float(real(ħ)))

# poles that lie on the ray arg(ζ) = θ (within an angular tolerance)
function _pole_on_ray(ζs, θ, tol)
    for ζ in ζs
        iszero(ζ) && continue
        Δ = mod(angle(ζ) - θ + π, 2π) - π
        abs(Δ) < tol && return ζ
    end
    nothing
end

"""
    laplace_sum(B::BorelSeries, ħ; θ = 0, order = nothing, rtol = eps^(3/4),
                check_poles = true) -> Complex

Laplace integral of the Borel transform along the ray ``\\arg ζ = θ``:

```math
φ(ħ) = a_0 + \\frac{1}{Γ(β)} \\int_0^{∞ e^{iθ}} e^{-ζ/ħ}\\, ζ^{β-1}\\, r(ζ)\\, dζ ,
```

where ``r`` is the Padé approximant of the reduced Borel function (`order` is passed
to [`pade`](@ref)) and ``a_0`` the split-off constant term. The result's precision
follows the type of `ħ` (pass a `BigFloat` for high precision). Throws
[`PoleOnRay`](@ref) when a Padé pole lies on the integration ray - a Stokes ray; use
[`lateral_sum`](@ref) there. Requires ``\\Re(e^{iθ}/ħ) > 0`` for convergence.

    laplace_sum(B::BorelSeries, r::AbstractBorelApproximant, ħ; θ = 0,
                rtol = eps^(3/4), check_poles = true) -> Complex

The same integral with an explicit (possibly prebuilt or non-Padé) approximant `r`
of the reduced Borel function - the seam every [`AbstractBorelApproximant`](@ref)
plugs into, e.g. `laplace_sum(B, conformal_borel(B), ħ)`.
"""
function laplace_sum(B::BorelSeries, ħ::Number; θ::Real = 0,
                     order::Union{Nothing,Integer} = nothing,
                     rtol::Union{Nothing,Real} = nothing,
                     check_poles::Bool = true)
    laplace_sum(B, pade(B; order), ħ; θ, rtol, check_poles)
end

function laplace_sum(B::BorelSeries, r::AbstractBorelApproximant, ħ::Number;
                     θ::Real = 0, rtol::Union{Nothing,Real} = nothing,
                     check_poles::Bool = true)
    T = _real_float_type(ħ)
    if check_poles
        bad = _pole_on_ray(poles(r), θ, sqrt(eps(T)))
        bad === nothing || throw(PoleOnRay(ComplexF64(bad), Float64(θ)))
    end
    _laplace_integral(B, r, ħ, T(θ), rtol === nothing ? eps(T)^(3//4) : rtol)
end

function _laplace_integral(B::BorelSeries, r::AbstractBorelApproximant, ħ::Number,
                           θ::Real, rtol::Real)
    T = _real_float_type(ħ)
    real(cis(θ) / ħ) > 0 ||
        throw(InvalidArgument("e^{-ζ/ħ} does not decay along θ = $θ for ħ = $ħ"))
    eiθ = cis(T(θ))
    β = T(B.beta)
    # ζ = t e^{iθ}:  ∫ e^{-t e^{iθ}/ħ} (t e^{iθ})^{β-1} r(t e^{iθ}) e^{iθ} dt
    integrand = t -> exp(-t * eiθ / ħ) * (t * eiθ)^(β - 1) * r(t * eiθ) * eiθ
    val, _ = quadgk(integrand, zero(T), T(Inf); rtol)
    Complex{T}(B.constant_term) + val / _gamma(β)
end

"""
    lateral_sum(B::BorelSeries, ħ; θ = 0, side = :plus, tilt = 1//100,
                order = nothing, rtol = eps^(3/4)) -> Complex

Lateral Borel sum on a Stokes ray: integrates along the tilted ray
``θ ± \\mathrm{tilt}`` (`side = :plus` above, `:minus` below). Because the Padé
approximant is meromorphic, the value is independent of `tilt` as long as no pole
lies between the tilted ray and ``θ`` - it *is* the lateral sum, not an
approximation to it. The `lateral_sum(B, r::AbstractBorelApproximant, ħ; …)`
method takes an explicit approximant, like [`laplace_sum`](@ref).
"""
function lateral_sum(B::BorelSeries, ħ::Number; θ::Real = 0, side::Symbol = :plus,
                     tilt::Real = 1//100, order::Union{Nothing,Integer} = nothing,
                     rtol::Union{Nothing,Real} = nothing)
    lateral_sum(B, pade(B; order), ħ; θ, side, tilt, rtol)
end

function lateral_sum(B::BorelSeries, r::AbstractBorelApproximant, ħ::Number;
                     θ::Real = 0, side::Symbol = :plus, tilt::Real = 1//100,
                     rtol::Union{Nothing,Real} = nothing)
    side in (:plus, :minus) ||
        throw(InvalidArgument("side must be :plus or :minus, got :$side"))
    tilt > 0 || throw(InvalidArgument("tilt must be > 0, got $tilt"))
    T = _real_float_type(ħ)
    θtilted = T(θ) + (side === :plus ? T(tilt) : -T(tilt))
    _laplace_integral(B, r, ħ, θtilted, rtol === nothing ? eps(T)^(3//4) : rtol)
end

"""
    stokes_discontinuity(B::BorelSeries, ħ; θ = 0, kwargs...) -> Complex

The Stokes discontinuity across the ray ``\\arg ζ = θ``:
``\\mathrm{disc}_θ = s_θ^+(ħ) - s_θ^-(ħ)``, the difference of the two
[`lateral_sum`](@ref)s. Sign convention (pinned by the Euler oracle test): for
``Φ = Σ_n n!\\,ħ^{n+1}`` (Borel pole at ``ζ = 1``) and real ``ħ > 0``,
`stokes_discontinuity(borel(Φ), ħ) ≈ +2πi\\, e^{-1/ħ}`.
"""
function stokes_discontinuity(B::BorelSeries, ħ::Number; θ::Real = 0, kwargs...)
    lateral_sum(B, ħ; θ, side = :plus, kwargs...) -
    lateral_sum(B, ħ; θ, side = :minus, kwargs...)
end

function stokes_discontinuity(B::BorelSeries, r::AbstractBorelApproximant,
                              ħ::Number; θ::Real = 0, kwargs...)
    lateral_sum(B, r, ħ; θ, side = :plus, kwargs...) -
    lateral_sum(B, r, ħ; θ, side = :minus, kwargs...)
end

"""
    borel_sum(Φ::FormalSeries, ħ; beta = power_offset(Φ), θ = 0, order = nothing,
              rtol = eps^(3/4)) -> Complex

Full Borel–Padé–Laplace pipeline convenience:
`laplace_sum(borel(Φ; beta), ħ; θ, order, rtol)`.
"""
function borel_sum(Φ::FormalSeries, ħ::Number;
                   beta::Union{Rational{Int},Integer} = power_offset(Φ),
                   kwargs...)
    laplace_sum(borel(Φ; beta), ħ; kwargs...)
end
