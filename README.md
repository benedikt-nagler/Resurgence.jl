# Resurgence.jl

[![Version](https://juliahub.com/docs/General/Resurgence/stable/version.svg)](https://juliahub.com/ui/Packages/General/Resurgence)
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://benedikt-nagler.github.io/Resurgence.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://benedikt-nagler.github.io/Resurgence.jl/dev)
[![CI](https://github.com/benedikt-nagler/Resurgence.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/benedikt-nagler/Resurgence.jl/actions/workflows/CI.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

Divergent and asymptotic series in Julia: Borel transforms, Padé approximants,
Borel–Padé–Laplace summation, large-order analysis, transseries and alien calculus.

A series with factorially growing coefficients has zero radius of convergence, but it still
determines a function. The Borel transform of $\Phi(\hbar) = \sum_n a_n \hbar^{n+\beta}$,

$$\hat{B}(\zeta) = \sum_{n} \frac{a_n}{\Gamma(n+\beta)}\, \zeta^{\,n+\beta-1},$$

converges on a disk, and the Laplace integral $\int_0^\infty e^{-\zeta/\hbar}\hat{B}\,d\zeta$
resums the series to a number. Each singularity of $\hat B$ at $A$ fixes the coefficient
growth $a_n \sim S\,\Gamma(n+b)/A^n$ and marks a scale $e^{-A/\hbar}$ that no order of the
series can see.

Given a list of coefficients you can get a number (`borel_sum`, `median_sum`, `accelerate`,
`hyper_sum`) or get structure (`large_order_fit`, `borel_pade_poles`). `Rational` input runs
through the whole pipeline exactly, and precision follows the caller's argument type.

This is an analysis engine for series you produced elsewhere, not a coefficient generator.
The built-in `FormalSeries(:name, n)` recursions are test oracles and worked examples.

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
large_order_fit(FormalSeries(:quartic, 80); order = 8)
                                  # A ≈ 1/3, the quartic instanton action

# Stokes phenomenon: lateral sums, discontinuity, median summation
Ψ = FormalSeries(:airy_bi, 40)    # Borel pole at +4/3, so ℝ₊ is a Stokes ray
lateral_sum(borel(Ψ), 0.3; side = :plus)   # complex
median_sum(Ψ, 0.3)                         # real

# Transseries and alien calculus
F = Transseries(:euler, 20)
alien_derivative(F; stokes = -2π * im)
```

## Features

- **Formal series.** `FormalSeries` over any coefficient type, with arithmetic, truncation,
  differentiation and fractional power offsets. `LogSeries` adds a polynomial in
  $\log\hbar$ for resonant problems.
- **Borel machinery.** Exact `borel` and `inverse_borel` in rising-factorial normalization.
- **Padé approximants.** Over any field, with pole and residue extraction.
  `borel_pade_poles` locates the Borel singularities. AAA rational approximation is
  available through an extension.
- **Summation.** `borel_sum`, `laplace_sum`, `lateral_sum` on tilted rays,
  `stokes_discontinuity`, `median_sum`. A guard catches poles on the integration ray.
- **Acceleration and conformal maps.** `accelerate` (Shanks, Wynn-ε, Levin-u/t, Richardson)
  works on any sequence. `conformal_borel` handles a cut instead of isolated poles.
- **Large-order analysis.** `large_order_fit` recovers $(A, b, S)$, including multi-saddle
  fits with complex-conjugate action pairs. `subtract_singularity` exposes the subleading
  singularity.
- **Hyperasymptotics.** `optimal_truncation`, `dingle_terminant`, and level-one Berry–Howls
  summation in `hyper_sum`.
- **Transseries.** One-parameter `Transseries` $\sum_n \sigma^n e^{-nA/\hbar}\Phi_n$, and
  `MultiTransseries` for several actions, stored sparsely over the action lattice.
- **Alien calculus.** `alien_derivative` via the bridge equation, `stokes_automorphism` as
  both a $\sigma$-shift and an operator, and `stokes_constant` from Borel residues or
  large-order growth.
- **Resonance.** `resonance_lattice` computes integer relations among the actions exactly
  over $\mathbb{Z}$, and `resonant_solve` generates the log tower.
- **ODEs and Painlevé.** `transseries_solve` builds a solution sector by sector from a
  user-supplied linear solve, and `multi_transseries_solve` handles several actions and
  resonant sectors. `painleve1`, `painleve1_two_parameter` and `painleve1_stokes_constant`
  apply this to Painlevé I.
- **Named examples.** `:euler`, `:airy`, `:airy_bi`, `:quartic` (Bender–Wu), `:phi4`,
  `:exp`, `:painleve1`, all from exact recursions.
- **Typed errors.** Invalid input throws `PoleOnRay`, `DegeneratePade` and friends instead
  of returning a wrong number.

## Extensions

The core depends only on QuadGK, PrecompileTools and LinearAlgebra.

```julia
using CairoMakie                   # or GLMakie
plot_borel_plane(borel(Φ))         # pole scatter with Laplace rays
plot_large_order(Φ)
plot_optimal_truncation(Φ, 0.1)

using BaryRational
aaa_borel(borel(Φ))                # AAA rational approximation
```

## Name collisions

`coefficients`, `derivative` and `evaluate` are also exported by AbstractAlgebra.jl and
other computer-algebra packages. If you load both into one namespace, import explicitly:

```julia
using Resurgence: FormalSeries, borel, coefficients
import AbstractAlgebra as AA
```

## Documentation

The [documentation](https://benedikt-nagler.github.io/Resurgence.jl/stable) has a manual page
per layer. The notebook
[divergent perturbation theory](examples/divergent_perturbation_theory.ipynb) walks through
the toolbox on the Euler, quartic-oscillator and Airy series.

## Related packages

This package depends on nothing below it. It is one foundation of a family of Julia packages
for exact and asymptotic methods:

- [ClusterAlgebras.jl](https://github.com/benedikt-nagler/ClusterAlgebras.jl): quivers, seed
  mutation, exchange graphs, green sequences.
- [ExactWKB.jl](https://github.com/benedikt-nagler/ExactWKB.jl): Borel sums the WKB series of
  a Schrödinger-type ODE with this package, and matches the result against cluster
  combinatorics.

## License

MIT
