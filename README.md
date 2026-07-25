# Resurgence.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://benedikt-nagler.github.io/Resurgence.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://benedikt-nagler.github.io/Resurgence.jl/dev)
[![CI](https://github.com/benedikt-nagler/Resurgence.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/benedikt-nagler/Resurgence.jl/actions/workflows/CI.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

Resurgence in Julia — divergent and asymptotic series made usable: Borel transforms, Padé
approximants, Borel–Padé–Laplace summation, sequence acceleration, large-order analysis,
transseries and alien calculus.

A series whose coefficients grow factorially has zero radius of convergence, so it is not a
function — but it still determines one, and it carries more information than a convergent
series would. Such expansions are the normal output of any problem where a small parameter
enters singularly: singular perturbation and boundary-layer problems, semiclassical and WKB
expansions, nonlinear ODEs near an irregular singular point (Painlevé and friends), the
coefficient asymptotics of generating functions, saddle-point expansions of integrals, and
perturbation theory in quantum mechanics and quantum field theory.

[Resurgence](https://en.wikipedia.org/wiki/Resurgent_function) is the theory that makes such
a series usable — Écalle's framework, and the one this package implements. Its basic tool is
the **Borel transform**: for $\Phi(\hbar) = \sum_n a_n \hbar^{n+\beta}$,

$$\hat{B}(\zeta) = \sum_{n} \frac{a_n}{\Gamma(n+\beta)}\, \zeta^{\,n+\beta-1}$$

converges on a disk, and the Laplace integral $\int_0^\infty e^{-\zeta/\hbar}\hat{B}(\zeta)\,d\zeta$
resums the original series to an actual number.

What makes this more than a summation trick is the geometry of the Borel plane. Each
singularity of $\hat B$ sits at some $A$, and that one location does three things at once: it
fixes the growth $a_n \sim S\,\Gamma(n+b)/A^n$ of the coefficients, it marks an exponentially
small scale $e^{-A/\hbar}$ that no order of the series can see, and — when it lies on the
integration ray — it makes the two lateral sums differ by exactly such a term, the **Stokes
phenomenon**. Divergence is therefore diagnostic: the growth rate of coefficients you already
have tells you about behaviour beyond all orders. (In a physics problem the $A$ are instanton
actions; in an ODE they are the exponents of the subdominant solutions; in coefficient
asymptotics they are the dominant singularities. The machinery does not care which.)

Two things you can do with a list of coefficients:

- **get a number** — `borel_sum` / `median_sum` for a resummed value, `accelerate` for
  sequence acceleration (Shanks, Wynn-ε, Levin, Richardson), `hyper_sum` when optimal
  truncation is not accurate enough;
- **get structure** — `large_order_fit` and `borel_pade_poles` recover the singularity
  locations, exponents and Stokes constants from the coefficients alone.

The notebook [resurgence in perturbation theory](examples/divergent_perturbation_theory.ipynb)
walks through the toolbox on the Euler, quartic-oscillator and Airy series.

Series with exact coefficients stay exact — `Rational` input runs through the whole pipeline
and precision follows the caller's argument type, so the same code runs in `Float64` or in
`BigFloat` at whatever `setprecision` is active.

**Scope.** This is an analysis engine for asymptotic series *however you produced them* — by
hand, from a recursion, from a symbolic computation, or read in from a file. It is not a
coefficient generator for any particular model: the handful of built-in
`FormalSeries(:name, n)` recursions exist as test oracles and worked examples, not as a
library of expansions.

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

## Functionality

**Formal series.** `FormalSeries` over any coefficient type, with arithmetic, truncation,
differentiation and fractional power offsets. `LogSeries` adds a polynomial in $\log\hbar$,
which is what resonant problems produce.

**Borel machinery.** Exact `borel` and `inverse_borel` in rising-factorial normalization.

**Padé approximants.** Over any field, with pole and residue extraction. `borel_pade_poles`
locates the Borel singularities of a series — the practical route from coefficients to the
actions $A$. AAA rational approximation is available as an alternative through an extension.

**Summation.** `borel_sum` and `laplace_sum` for Borel–Padé–Laplace sums, `lateral_sum` on
tilted rays, `stokes_discontinuity` for the jump across a Stokes line, and `median_sum` for
the real average. A guard detects poles sitting on the integration ray instead of silently
returning nonsense.

**Acceleration and conformal maps.** `accelerate` offers Shanks, Wynn-ε, Levin-u/t and
Richardson on any sequence (`partial_sums` turns a series into one), useful on its own as a
convergence accelerator; `conformal_borel` maps the Borel plane for series whose
singularities form a cut rather than isolated poles, where Padé struggles.

**Large-order analysis.** `large_order_fit` recovers $(A, b, S)$ from
$a_n \sim S\,\Gamma(n+b)/A^n$, including multi-saddle fits with complex-conjugate action
pairs — the numerical-analysis side of Darboux/singularity analysis, run directly on the
coefficients. `subtract_singularity` peels off a known singularity to expose the subleading
one.

**Hyperasymptotics.** Optimal truncation, Dingle terminants, and level-one Berry–Howls
hyperasymptotic summation (`optimal_truncation`, `dingle_terminant`, `hyper_sum`).

**Transseries.** One-parameter `Transseries` $\sum_n \sigma^n e^{-nA/\hbar}\Phi_n$, and
`MultiTransseries` for several actions at once, stored sparsely over the action lattice.

**Alien calculus.** `alien_derivative` via the bridge equation (single- and multi-parameter,
by lattice charge or by singularity location), `stokes_automorphism` as both a $\sigma$-shift
and an operator, and `stokes_constant` from either Borel-plane residues or large-order
growth.

**Resonance.** When the actions satisfy integer relations the transseries develops logarithms.
`resonance_lattice` computes the relations exactly over $\mathbb{Z}$, and `resonant_solve` is
the primitive that generates the log tower.

**ODEs and Painlevé.** `transseries_solve` is a model-agnostic driver that builds a
transseries solution sector by sector from a user-supplied linear solve. `painleve1` applies
it to Painlevé I, with both normalization conventions of the instanton action available
(`:gikm` by default, `:string` for the string-equation convention).

**Named examples.** `FormalSeries(:euler, n)` plus `:airy`, `:airy_bi`, `:quartic`
(Bender–Wu ground-state energy), `:phi4`, `:exp`, `:painleve1`, all from exact recursions.
They double as test oracles — the suite checks mathematical identities rather than stored
numbers.

**Numerical hygiene.** Every operation returns a new object, invalid input throws a typed
error (`PoleOnRay`, `DegeneratePade`, …) rather than a silently wrong number, and precision
is taken from the argument type rather than a global setting.

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

Three exports — `coefficients`, `derivative`, `evaluate` — are also exported by
AbstractAlgebra.jl and by several other computer-algebra packages. They are the right names
for what they do, so they are kept; if you load both packages into one namespace, import
explicitly rather than relying on `using`:

```julia
using Resurgence: FormalSeries, borel, coefficients
import AbstractAlgebra as AA
```

## Related packages

This package is self-contained and depends on nothing below. It is also one foundation of a
family of Julia packages for **exact and asymptotic methods** — where a discrete, exactly
computable structure controls a continuous, only-asymptotically-defined one:

- [ClusterAlgebras.jl](https://github.com/benedikt-nagler/ClusterAlgebras.jl) — the discrete
  side: quivers, seed mutation, exchange graphs, green sequences.
- [ExactWKB.jl](https://github.com/benedikt-nagler/ExactWKB.jl) — the bridge, which Borel
  sums the WKB series of a Schrödinger-type ODE with this package and matches the result
  against cluster combinatorics via the Iwaki–Nakanishi dictionary.

## License

MIT
