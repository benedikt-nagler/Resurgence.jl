```@meta
CurrentModule = Resurgence
```

# Resurgence.jl

Divergent and asymptotic series: Borel transforms, Padé approximants,
Borel–Padé–Laplace summation, sequence acceleration, large-order analysis, transseries and
alien calculus.

A series whose coefficients grow factorially has zero radius of convergence, so it is not a
function, but it still determines one, and it carries more information than a convergent
series would. Such expansions are the normal output of any problem where a small parameter
enters singularly: singular perturbation and boundary-layer problems, semiclassical and WKB
expansions, nonlinear ODEs near an irregular singular point, the coefficient asymptotics of
generating functions, saddle-point expansions of integrals, and perturbation theory in
quantum mechanics and quantum field theory.

Écalle's theory of resurgence is what this package implements. Its basic tool is the Borel
transform: for ``\Phi(\hbar) = \sum_n a_n \hbar^{n+\beta}``,

```math
\hat{B}(\zeta) = \sum_{n} \frac{a_n}{\Gamma(n+\beta)}\, \zeta^{\,n+\beta-1}
```

converges on a disk, and the Laplace integral
``\int_0^\infty e^{-\zeta/\hbar}\hat{B}(\zeta)\,d\zeta`` resums the original series to an
actual number.

The content of the theory is the geometry of the Borel plane. Each singularity of ``\hat B``
sits at some ``A``, and that one location fixes three things: the coefficient growth
``a_n \sim S\,\Gamma(n+b)/A^n``, an exponentially small scale ``e^{-A/\hbar}`` that no order
of the series can see, and, when the singularity lies on the integration ray, a difference of
exactly that size between the two lateral sums, which is the Stokes phenomenon. The growth
rate of the coefficients you already have therefore tells you about behaviour beyond all
orders. In a physics problem the ``A`` are instanton actions, in an ODE they are the exponents
of the subdominant solutions, and in coefficient asymptotics they are the dominant
singularities. The machinery is the same in each case.

Two things you can do with a list of coefficients:

- **get a number**: [`borel_sum`](@ref) / [`median_sum`](@ref) for a resummed value,
  [`accelerate`](@ref) for sequence acceleration, [`hyper_sum`](@ref) when optimal truncation
  is not accurate enough;
- **get structure**: [`large_order_fit`](@ref) and [`borel_pade_poles`](@ref) recover the
  singularity locations, exponents and Stokes constants from the coefficients alone.

Series with exact coefficients stay exact: `Rational` input runs through the whole pipeline,
and precision follows the caller's argument type, so the same code runs in `Float64` or in
`BigFloat` at whatever `setprecision` is active.

## Scope

This is an analysis engine for asymptotic series *however you produced them*: by hand, from
a recursion, from a symbolic computation, or read in from a file. It is not a coefficient
generator for any particular model. The handful of built-in
[`FormalSeries`](@ref)`(:name, n)` recursions exist as test oracles and worked examples, not
as a library of expansions.

## Installation

From the Julia REPL:

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

Three exports, [`coefficients`](@ref), [`derivative`](@ref) and [`evaluate`](@ref), are also
exported by AbstractAlgebra.jl and by several other computer-algebra packages. They are the
right names for what they do, so they are kept. If you load both packages into one namespace,
import explicitly instead of relying on `using`:

```julia
using Resurgence: FormalSeries, borel, coefficients
import AbstractAlgebra as AA
```

## Related packages

This package is self-contained and depends on nothing below. It is also one foundation of a
family of Julia packages for exact and asymptotic methods, in which a discrete, exactly
computable structure controls a continuous, only-asymptotically defined one.
[ClusterAlgebras.jl](https://github.com/benedikt-nagler/ClusterAlgebras.jl) is the discrete
side, and [ExactWKB.jl](https://github.com/benedikt-nagler/ExactWKB.jl) is the bridge: it
Borel sums the WKB series of a Schrödinger-type ODE with this package and matches the result
against cluster combinatorics via the Iwaki–Nakanishi dictionary.

```@docs
Resurgence
```
