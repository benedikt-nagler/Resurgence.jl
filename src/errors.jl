# Typed error hierarchy. Every invalid operation throws one of these; messages start
# with the error-type name so failures are self-identifying in test output.

"""
    ResurgenceError

Abstract supertype of all typed errors thrown by Resurgence.jl.
"""
abstract type ResurgenceError <: Exception end

"""
    InvalidArgument(msg)

Generic escape hatch for invalid inputs not covered by a more specific error type.
"""
struct InvalidArgument <: ResurgenceError
    msg::String
end

Base.showerror(io::IO, e::InvalidArgument) = print(io, "InvalidArgument: ", e.msg)

"""
    IncompatibleSeries(field, lhs, rhs)

Thrown by series arithmetic when the two operands disagree in `field`
(e.g. `:var`, or a `:power_offset` mismatch that is not an integer shift).
"""
struct IncompatibleSeries <: ResurgenceError
    field::Symbol
    lhs::Any
    rhs::Any
end

function Base.showerror(io::IO, e::IncompatibleSeries)
    print(io, "IncompatibleSeries: operands differ in $(e.field): ",
          "$(e.lhs) vs $(e.rhs)")
end

"""
    DegeneratePade(L, M)

Thrown when the Toeplitz system of the `[L/M]` Padé approximant is singular — the
underlying function admits a lower-degree rational representation. Request a lower
order, or pass `reduce = true` to [`pade`](@ref) for automatic degree reduction.
"""
struct DegeneratePade <: ResurgenceError
    L::Int
    M::Int
end

function Base.showerror(io::IO, e::DegeneratePade)
    print(io, "DegeneratePade: the [$(e.L)/$(e.M)] Padé system is singular; ",
          "the series is matched by a lower-degree rational function — ",
          "request a lower order or pass `reduce = true` for automatic reduction")
end

"""
    NoSingularityFound(msg)

Thrown when a Borel-plane singularity is required but none could be detected —
e.g. [`conformal_borel`](@ref) or [`subtract_singularity`](@ref) on a series whose
Padé approximant has no poles. Pass the singularity position explicitly.
"""
struct NoSingularityFound <: ResurgenceError
    msg::String
end

Base.showerror(io::IO, e::NoSingularityFound) =
    print(io, "NoSingularityFound: ", e.msg)

"""
    PoleOnRay(pole, theta)

Thrown by Laplace summation when a Borel-plane pole lies (within tolerance) on the
integration ray at angle `theta` — a Stokes ray. Use [`lateral_sum`](@ref) with
`side = :plus` or `:minus` to integrate around it.
"""
struct PoleOnRay <: ResurgenceError
    pole::Number
    theta::Float64
end

function Base.showerror(io::IO, e::PoleOnRay)
    print(io, "PoleOnRay: Borel-plane pole at ζ ≈ $(e.pole) lies on the Laplace ray ",
          "θ = $(e.theta) (a Stokes ray); use `lateral_sum(B, ħ; side=:plus)` or ",
          "`side=:minus` to integrate around it")
end
