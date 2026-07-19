# Resurgence.jl

[![CI](https://github.com/benedikt-nagler/Resurgence.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/benedikt-nagler/Resurgence.jl/actions/workflows/CI.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

Resurgence in Julia: Borel transforms, Padé approximants, Borel–Padé–Laplace summation,
large-order analysis, transseries and alien calculus.

Perturbative expansions in quantum mechanics and quantum field theory usually have zero
radius of convergence. The coefficients grow factorially, so the series is not a function —
but it still determines one. [Resurgence
theory](https://en.wikipedia.org/wiki/Resurgent_function), due to Écalle, recovers it through
the **Borel transform**. For $\Phi(\hbar) = \sum_n a_n \hbar^{n+\beta}$,

$$\hat{B}(\zeta) = \sum_{n} \frac{a_n}{\Gamma(n+\beta)}\, \zeta^{\,n+\beta-1}$$

converges on a disk, and the Laplace integral $\int_0^\infty e^{-\zeta/\hbar}\hat{B}(\zeta)\,d\zeta$
resums the original series.

The point of the construction is what sits in the Borel plane. The singularities of $\hat B$
are at the instanton actions $A$ of the problem: they fix the growth
$a_n \sim S\,\Gamma(n+b)/A^n$ of the coefficients, and they correspond to the exponentially
small effects $e^{-A/\hbar}$ that the perturbative series cannot see. When a singularity lies
on the integration ray, the two lateral sums differ by such a term — the **Stokes
phenomenon**. So the divergence of the series is not a defect; it encodes the nonperturbative
content.

These structures turn up in anharmonic oscillators and double wells, in exact WKB and
spectral problems, and in matrix models and topological strings. The notebook [resurgence in
perturbation theory](examples/divergent_perturbation_theory.ipynb) walks through the toolbox
on the Euler, quartic-oscillator and Airy series.

Series with exact coefficients stay exact — `Rational` input runs through the whole pipeline
and precision follows the caller's argument type, so the same code runs in `Float64` or in
`BigFloat` at whatever `setprecision` is active.

## Installation

Not yet registered in General. From the Julia REPL:

```julia
pkg> add https://github.com/benedikt-nagler/Resurgence.jl
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

## Functionality

**Formal series.** `FormalSeries` over any coefficient type, with arithmetic, truncation,
differentiation and fractional power offsets. `LogSeries` adds a polynomial in $\log\hbar$,
which is what resonant problems produce.

**Borel machinery.** Exact `borel` and `inverse_borel` in rising-factorial normalization.

**Padé approximants.** Over any field, with pole and residue extraction. `borel_pade_poles`
locates the Borel singularities of a series, which is the practical route to the instanton
actions.

**Summation.** `borel_sum` and `laplace_sum` for Borel–Padé–Laplace sums, `lateral_sum` on
tilted rays, `stokes_discontinuity` for the jump across a Stokes line, and `median_sum` for
the real average. A guard detects poles sitting on the integration ray instead of silently
returning nonsense.

**Acceleration and conformal maps.** `accelerate` offers Shanks, Wynn-ε, Levin-u/t and
Richardson; `conformal_borel` maps the Borel plane for series whose singularities form a cut
rather than isolated poles, where Padé struggles.

**Large-order analysis.** `large_order_fit` recovers $(A, b, S)$ from
$a_n \sim S\,\Gamma(n+b)/A^n$, including multi-saddle fits with complex-conjugate action
pairs. `subtract_singularity` peels off a known singularity to expose the subleading one.

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

## Extensions

The core depends only on QuadGK; the rest is loaded on demand.

```julia
using CairoMakie                   # or GLMakie
plot_borel_plane(borel(Φ))         # pole scatter with Laplace rays
plot_large_order(Φ)
plot_optimal_truncation(Φ, 0.1)

using BaryRational
aaa_borel(borel(Φ))                # AAA rational approximation of the Borel function
```

## Related packages

[ExactWKB.jl](https://github.com/benedikt-nagler/ExactWKB.jl) uses this package to Borel sum
the Voros series of one-dimensional Schrödinger problems, and connects them to
[ClusterAlgebras.jl](https://github.com/benedikt-nagler/ClusterAlgebras.jl) through the
Iwaki–Nakanishi dictionary. Resurgence.jl itself has no dependency on either.

## References

- F. J. Dyson, *Divergence of perturbation theory in quantum electrodynamics*, Phys. Rev. **85** (1952), 631–632.
- C. M. Bender, T. T. Wu, *Anharmonic oscillator*, Phys. Rev. **184** (1969), 1231–1260.
- J. Écalle, *Les fonctions résurgentes* I–III, Publ. Math. d'Orsay (1981–1985).
- A. Voros, *The return of the quartic oscillator: the complex WKB method*, Ann. Inst. H. Poincaré A **39** (1983), 211–338.
- M. V. Berry, C. J. Howls, *Hyperasymptotics for integrals with saddles*, Proc. R. Soc. Lond. A **434** (1991), 657–675.
- I. Aniceto, G. Başar, R. Schiappa, *A primer on resurgent transseries and their asymptotics*, Phys. Rept. **809** (2019), 1–135.
- D. Dorigoni, *An introduction to resurgence, trans-series and alien calculus*, Ann. Phys. **409** (2019), 167914.

## License

MIT
