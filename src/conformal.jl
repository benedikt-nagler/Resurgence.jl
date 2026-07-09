# Conformal maps of the Borel plane: map the cut ζ-plane to the unit disk, re-expand
# the reduced Borel series in the map variable w, and approximate there — the
# standard resummation accelerator when Padé alone converges slowly. Composable with
# the Laplace pipeline through the AbstractBorelApproximant seam.

"""
    ConformalMap{T}

The conformal map of the ζ-plane cut along ``\\{t\\,ζ_0 : t ≥ 1\\}`` (and, with
`pair = true`, also along ``\\{-t\\,ζ_0 : t ≥ 1\\}``) onto the unit w-disk, with
``0 ↦ 0`` and the branch point(s) mapped to ``|w| = 1``:

- single cut: ``w(ζ) = \\dfrac{1 - \\sqrt{1 - ζ/ζ_0}}{1 + \\sqrt{1 - ζ/ζ_0}}``,
  inverse ``ζ(w) = 4ζ_0\\, w/(1+w)^2``;
- symmetric pair (``±ζ_0``): ``w(ζ) = (1 - \\sqrt{1 - ζ^2/ζ_0^2})\\,ζ_0/ζ``,
  inverse ``ζ(w) = 2ζ_0\\, w/(1+w^2)``.

Non-symmetric singularity pairs at arbitrary angles are out of scope — use the AAA
extension or the multi-saddle [`large_order_fit`](@ref) in that regime. Callable:
`m(ζ)` evaluates the forward map (in complex arithmetic); [`inverse`](@ref) maps
back. Construct with [`conformal_map`](@ref) or via [`conformal_borel`](@ref).
"""
struct ConformalMap{T}
    zeta0::T
    pair::Bool
end

"""
    conformal_map(zeta0; pair = false) -> ConformalMap

The [`ConformalMap`](@ref) placing the branch point at `zeta0` (and, with
`pair = true`, its mirror at `-zeta0`).
"""
function conformal_map(zeta0::Number; pair::Bool = false)
    iszero(zeta0) && throw(InvalidArgument("the conformal-map branch point must be ≠ 0"))
    ConformalMap(zeta0, pair)
end

function (m::ConformalMap)(ζ::Number)
    r = ζ / m.zeta0
    if m.pair
        u = sqrt(complex(1 - r^2))
        iszero(r) ? zero(u) : (1 - u) / r
    else
        u = sqrt(complex(1 - r))
        (1 - u) / (1 + u)
    end
end

"""
    inverse(m::ConformalMap, w) -> Number

The inverse map ``w ↦ ζ``: ``4ζ_0 w/(1+w)^2`` (single cut) or ``2ζ_0 w/(1+w^2)``
(symmetric pair).
"""
inverse(m::ConformalMap, w::Number) =
    m.pair ? 2 * m.zeta0 * w / (1 + w^2) : 4 * m.zeta0 * w / (1 + w)^2

# Taylor coefficients of ζ(w) to order N (vector of length N + 1, zero constant
# term): single cut 4ζ₀ Σ_{k≥1} (-1)^{k+1} k wᵏ; pair 2ζ₀ Σ_{k≥0} (-1)^k w^{2k+1}.
# Exact when ζ₀ is exact.
function _map_series(m::ConformalMap, N::Integer)
    z = [zero(m.zeta0) for _ in 1:(N + 1)]
    if m.pair
        for k in 0:((N - 1) ÷ 2)
            z[2k + 2] = (isodd(k) ? -2 : 2) * m.zeta0
        end
    else
        for k in 1:N
            z[k + 1] = (isodd(k) ? 4k : -4k) * m.zeta0
        end
    end
    z
end

# c = a·b truncated to N coefficients (plain Cauchy product on raw vectors).
function _mul_trunc(a::Vector, b::Vector, N::Integer)
    c = [zero(a[1] * b[1]) for _ in 1:N]
    for i in 1:min(length(a), N)
        iszero(a[i]) && continue
        for j in 1:min(length(b), N - i + 1)
            c[i + j - 1] += a[i] * b[j]
        end
    end
    c
end

# Composition c(w) = b(z(w)) truncated to length(b) coefficients, by Horner in
# series space; requires z[1] == 0 (which _map_series guarantees).
function _compose_series(b::Vector, z::Vector)
    N = length(b)
    acc = [zero(b[1] * z[1]) for _ in 1:N]      # promoted element type
    acc[1] += b[N]
    for n in (N - 1):-1:1
        acc = _mul_trunc(acc, z, N)
        acc[1] += b[n]
    end
    acc
end

"""
    ConformalPade{T} <: AbstractBorelApproximant

A conformally mapped Borel approximant: the reduced Borel series re-expanded in the
map variable ``w`` of a [`ConformalMap`](@ref) and approximated there by an inner
[`PadeApproximant`](@ref). Callable at ζ (`c(ζ) = inner(w(ζ))`); [`poles`](@ref)
reports the branch point(s) plus the physical-sheet (``|w| < 1``) inner poles
pulled back to the ζ-plane. Construct with [`conformal_borel`](@ref); plugs into
[`laplace_sum`](@ref) through the approximant seam.
"""
struct ConformalPade{T} <: AbstractBorelApproximant
    map::ConformalMap
    inner::PadeApproximant{T}
    var::Symbol
end

variable(c::ConformalPade) = c.var

(c::ConformalPade)(ζ::Number) = c.inner(c.map(ζ))

"""
    conformal_borel(B::BorelSeries; zeta0 = nothing, pair = false, order = nothing,
                    inner_order = nothing) -> ConformalPade

Build the conformally mapped approximant of the reduced Borel function of `B`:
place the branch point at `zeta0` (default: the smallest-modulus Borel–Padé pole,
computed with `pade(B; order, reduce = true)`; throws [`NoSingularityFound`](@ref)
when there is none), re-expand the series in the map variable ``w`` (exact for
exact `zeta0` and coefficients), and approximate in ``w``: near-diagonally by
default, `inner_order = m` for a diagonal ``[m/m]``, and `inner_order = 0` for the
pure re-expanded polynomial (already convergent on the disk). Sum with
`laplace_sum(B, conformal_borel(B), ħ)`.
"""
function conformal_borel(B::BorelSeries; zeta0::Union{Nothing,Number} = nothing,
                         pair::Bool = false,
                         order::Union{Nothing,Integer} = nothing,
                         inner_order::Union{Nothing,Integer} = nothing)
    ζ0 = zeta0
    if ζ0 === nothing
        ζs = poles(pade(B; order, reduce = true))
        isempty(ζs) &&
            throw(NoSingularityFound("conformal_borel: the Borel–Padé approximant " *
                                     "of the series has no poles; pass zeta0"))
        ζ0 = ζs[1]
    end
    m = conformal_map(ζ0; pair)
    b = B.series.coeffs
    N = length(b)
    c = _compose_series(b, _map_series(m, N - 1))
    inner = if inner_order === nothing
        M = (N - 1) ÷ 2
        pade(c, N - 1 - M, M; var = variable(B), reduce = true)
    elseif iszero(inner_order)
        pade(c, N - 1, 0; var = variable(B))
    else
        pade(c, Int(inner_order), Int(inner_order); var = variable(B), reduce = true)
    end
    ConformalPade(m, inner, variable(B))
end

function poles(c::ConformalPade; refine::Bool = true)
    ζ0 = Complex{BigFloat}(c.map.zeta0)
    out = c.map.pair ? [ζ0, -ζ0] : [ζ0]
    for w in poles(c.inner; refine)
        abs(w) < 1 || continue                    # |w| ≥ 1 is the second sheet
        push!(out, Complex{BigFloat}(inverse(c.map, w)))
    end
    sort!(out; by = abs)
end
