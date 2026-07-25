```@meta
CurrentModule = Resurgence
```

# Resurgence.jl

Borel transforms, Padé approximants, Borel–Padé–Laplace summation, large-order analysis,
transseries and alien calculus.

Perturbative expansions in quantum mechanics and quantum field theory usually have zero
radius of convergence. The coefficients grow factorially, so the series is not a function —
but it still determines one. Resurgence theory, due to Écalle, recovers it through the Borel
transform. For ``\Phi(\hbar) = \sum_n a_n \hbar^{n+\beta}``,

```math
\hat{B}(\zeta) = \sum_{n} \frac{a_n}{\Gamma(n+\beta)}\, \zeta^{\,n+\beta-1}
```

converges on a disk, and the Laplace integral
``\int_0^\infty e^{-\zeta/\hbar}\hat{B}(\zeta)\,d\zeta`` resums the original series.

The point of the construction is what sits in the Borel plane. The singularities of
``\hat B`` are at the instanton actions ``A`` of the problem: they fix the growth
``a_n \sim S\,\Gamma(n+b)/A^n`` of the coefficients, and they correspond to the exponentially
small effects ``e^{-A/\hbar}`` that the perturbative series cannot see. When a singularity
lies on the integration ray, the two lateral sums differ by such a term — the Stokes
phenomenon. So the divergence of the series is not a defect; it encodes the nonperturbative
content.

Series with exact coefficients stay exact: `Rational` input runs through the whole pipeline,
and precision follows the caller's argument type, so the same code runs in `Float64` or in
`BigFloat` at whatever `setprecision` is active.

## Scope

This is an analysis engine for asymptotic series *however you produced them* — by hand, from
a recursion, from your own field-theory or string computation, or read in from a file. It is
not a coefficient generator for any particular model: the handful of built-in
[`FormalSeries`](@ref)`(:name, n)` recursions exist as test oracles and worked examples, not
as a library of physical expansions.

## Installation

```julia
pkg> add Resurgence
```

Requires Julia 1.10 or later.

## Quick start

```julia
using Resurgence

# The Euler series, exact rationals; divergent for every ħ ≠ 0
Φ = FormalSeries(:euler, 40)

# Borel transform, exactly −1/(1+ζ)
B = borel(Φ)
poles(pade(B; order = 1))         # ζ = −1, the Borel singularity

# Borel–Padé–Laplace summation; precision follows the argument type
borel_sum(Φ, 0.1; order = 1)      # ≈ −0.0915633, i.e. −e^{10}E₁(10)
borel_sum(Φ, big"0.1"; order = 1) # same pipeline in BigFloat

# Large-order analysis: aₙ ~ S·Γ(n+b)/Aⁿ
large_order_fit(Φ)                # (A, b, S) = (1, 0, 1) for Euler

# Stokes phenomenon: lateral sums, discontinuity, median summation
Ψ = FormalSeries(:airy_bi, 40)    # Borel pole at +4/3, so ℝ₊ is a Stokes ray
lateral_sum(borel(Ψ), 0.3; side = :plus)   # complex
median_sum(Ψ, 0.3)                         # real

# Transseries and alien calculus
F = Transseries(:euler, 20)
alien_derivative(F; stokes = -2π * im)
```

## A note on name collisions

Three exports — [`coefficients`](@ref), [`derivative`](@ref), [`evaluate`](@ref) — are also
exported by AbstractAlgebra.jl and by several other computer-algebra packages. They are the
right names for what they do, so they are kept; if you load both packages into one namespace,
import explicitly rather than relying on `using`:

```julia
using Resurgence: FormalSeries, borel, coefficients
import AbstractAlgebra as AA
```

## Related packages

[ExactWKB.jl](https://github.com/benedikt-nagler/ExactWKB.jl) uses this package to Borel sum
the Voros series of one-dimensional Schrödinger problems, and connects them to
[ClusterAlgebras.jl](https://github.com/benedikt-nagler/ClusterAlgebras.jl) through the
Iwaki–Nakanishi dictionary. Resurgence.jl itself has no dependency on either.

```@docs
Resurgence
```
