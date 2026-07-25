# Resurgence.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://benedikt-nagler.github.io/Resurgence.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://benedikt-nagler.github.io/Resurgence.jl/dev)
[![CI](https://github.com/benedikt-nagler/Resurgence.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/benedikt-nagler/Resurgence.jl/actions/workflows/CI.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

Resurgence in Julia - divergent and asymptotic series made usable: Borel transforms, Padé
approximants, Borel–Padé–Laplace summation, sequence acceleration, large-order analysis,
transseries and alien calculus.

A series whose coefficients grow factorially has zero radius of convergence, so it is not a
function - but it still determines one.
[Resurgence](https://en.wikipedia.org/wiki/Resurgent_function) is Écalle's theory of how, and
its basic tool is the **Borel transform**: for $\Phi(\hbar) = \sum_n a_n \hbar^{n+\beta}$,

$$\hat{B}(\zeta) = \sum_{n} \frac{a_n}{\Gamma(n+\beta)}\, \zeta^{\,n+\beta-1}$$

converges on a disk, and the Laplace integral
$\int_0^\infty e^{-\zeta/\hbar}\hat{B}(\zeta)\,d\zeta$ resums the original series to a number.
Each singularity of $\hat B$ at some $A$ then fixes the coefficient growth
$a_n \sim S\,\Gamma(n+b)/A^n$, marks an exponentially small scale $e^{-A/\hbar}$ invisible to
every order of the series, and - on the integration ray - produces the **Stokes phenomenon**.
Divergence is diagnostic: the coefficients you already have tell you about behaviour beyond
all orders.

So from a list of coefficients you can **get a number** (`borel_sum`, `median_sum`,
`accelerate`, `hyper_sum`) or **get structure** (`large_order_fit`, `borel_pade_poles` recover
actions, exponents and Stokes constants). Exact input stays exact - `Rational` runs through
the whole pipeline, and precision follows the caller's argument type.

**Scope.** An analysis engine for asymptotic series *however you produced them*, not a
coefficient generator: the built-in `FormalSeries(:name, n)` recursions are test oracles and
worked examples, not a library of expansions.

## Installation

Not yet registered in General. From the Julia REPL:

```julia
pkg> add https://github.com/benedikt-nagler/Resurgence.jl
```

Requires Julia 1.10 or later. Full API documentation is at
[benedikt-nagler.github.io/Resurgence.jl](https://benedikt-nagler.github.io/Resurgence.jl/stable).

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

The notebook [resurgence in perturbation theory](examples/divergent_perturbation_theory.ipynb)
walks through the toolbox on the Euler, quartic-oscillator and Airy series.

## Functionality

**Formal series.** `FormalSeries` over any coefficient type, with arithmetic, truncation,
differentiation and fractional power offsets; `LogSeries` adds a polynomial in $\log\hbar$ for
resonant problems.

**Borel machinery.** Exact `borel` and `inverse_borel` in rising-factorial normalization.

**Padé approximants.** Over any field, with pole and residue extraction. `borel_pade_poles`
locates the Borel singularities - the practical route from coefficients to the actions $A$.
AAA rational approximation is available as an alternative through an extension.

**Summation.** `borel_sum` and `laplace_sum`, `lateral_sum` on tilted rays,
`stokes_discontinuity` for the jump across a Stokes line, `median_sum` for the real average. A
guard catches poles sitting on the integration ray instead of returning nonsense.

**Acceleration and conformal maps.** `accelerate` (Shanks, Wynn-ε, Levin-u/t, Richardson) on
any sequence, useful on its own; `conformal_borel` for series whose singularities form a cut
rather than isolated poles, where Padé struggles.

**Large-order analysis.** `large_order_fit` recovers $(A, b, S)$ from
$a_n \sim S\,\Gamma(n+b)/A^n$, including multi-saddle fits with complex-conjugate action
pairs; `subtract_singularity` peels off a known singularity to expose the subleading one.

**Hyperasymptotics.** Optimal truncation, Dingle terminants, and level-one Berry–Howls
hyperasymptotic summation (`optimal_truncation`, `dingle_terminant`, `hyper_sum`).

**Transseries.** One-parameter `Transseries` $\sum_n \sigma^n e^{-nA/\hbar}\Phi_n$, and
`MultiTransseries` for several actions at once, stored sparsely over the action lattice.

**Alien calculus.** `alien_derivative` via the bridge equation (by lattice charge or by
singularity location), `stokes_automorphism` as both a $\sigma$-shift and an operator, and
`stokes_constant` from Borel-plane residues or large-order growth.

**Resonance.** Integer relations among the actions make the transseries develop logarithms;
`resonance_lattice` computes them exactly over $\mathbb{Z}$ and `resonant_solve` generates the
log tower.

**ODEs and Painlevé.** `transseries_solve` builds a transseries solution sector by sector from
a user-supplied linear solve; `painleve1` applies it to Painlevé I, in either normalization
convention for the instanton action.

**Named examples.** `:euler`, `:airy`, `:airy_bi`, `:quartic` (Bender–Wu), `:phi4`, `:exp`,
`:painleve1`, all from exact recursions. They double as test oracles - the suite checks
mathematical identities rather than stored numbers.

**Numerical hygiene.** Every operation returns a new object, invalid input throws a typed
error (`PoleOnRay`, `DegeneratePade`, …) rather than a silently wrong number, and precision
comes from the argument type rather than a global setting.

## Extensions

The core depends only on QuadGK, PrecompileTools and the LinearAlgebra stdlib; the rest is
loaded on demand.

```julia
using CairoMakie                   # or GLMakie
plot_borel_plane(borel(Φ))         # pole scatter with Laplace rays
plot_large_order(Φ)
plot_optimal_truncation(Φ, 0.1)

using BaryRational
aaa_borel(borel(Φ))                # AAA rational approximation of the Borel function
```

## A note on name collisions

Three exports - `coefficients`, `derivative`, `evaluate` - are also exported by
AbstractAlgebra.jl and by several other computer-algebra packages. They are the right names
for what they do, so they are kept; if you load both packages into one namespace, import
explicitly rather than relying on `using`:

```julia
using Resurgence: FormalSeries, borel, coefficients
import AbstractAlgebra as AA
```

## Related packages

This package is self-contained and depends on nothing below. It is also one foundation of a
family of Julia packages for **exact and asymptotic methods** - where a discrete, exactly
computable structure controls a continuous, only-asymptotically-defined one:

- [ClusterAlgebras.jl](https://github.com/benedikt-nagler/ClusterAlgebras.jl) - the discrete
  side: quivers, seed mutation, exchange graphs, green sequences.
- [ExactWKB.jl](https://github.com/benedikt-nagler/ExactWKB.jl) - the bridge, which Borel
  sums the WKB series of a Schrödinger-type ODE with this package and matches the result
  against cluster combinatorics via the Iwaki–Nakanishi dictionary.

## License

MIT
