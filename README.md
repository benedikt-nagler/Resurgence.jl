# Resurgence.jl

A Julia package for resurgence: Borel summation, Stokes phenomena, and large-order analysis of divergent series.

## Introduction

Perturbative expansions in quantum mechanics and quantum field theory typically have zero radius of convergence. For example, for the quartic anharmonic oscillator the coefficients grow like $n!\,3^n$.

[Resurgence theory](https://en.wikipedia.org/wiki/Resurgent_function) (Écalle) handles such series through the **Borel transform**. For $\Phi(\hbar) = \sum_n a_n \hbar^{n+\beta}$,

$$\hat{B}(\zeta) = \sum_{n} \frac{a_n}{\Gamma(n+\beta)}\, \zeta^{\,n+\beta-1}$$

converges on a disk, and the Laplace integral $\int_0^\infty e^{-\zeta/\hbar}\,\hat{B}(\zeta)\,d\zeta$ resums the series. The singularities of $\hat{B}$ sit at the instanton actions $A$ of the problem: they determine the growth $a_n \sim S\,\Gamma(n+b)/A^n$ of the coefficients and correspond to the exponentially small effects $e^{-A/\hbar}$ that the series misses. A singularity on the integration ray makes the two lateral sums differ by such a term, and this is called the **Stokes phenomenon**.

These structures show up in many places: anharmonic oscillators and double wells, exact WKB and spectral problems, matrix models and topological strings, Painlevé equations. The notebook [Divergent series in physics](examples/divergent_perturbation_theory.ipynb) works through the basic toolbox of the package on the Euler, quartic-oscillator, and Airy series.

## Features

- **Formal series**: `FormalSeries` over any coefficient type, with arithmetic and fractional power offsets
- **Borel machinery**: exact Borel and inverse Borel transforms (`borel`, `inverse_borel`) in rising-factorial normalization
- **Padé approximants**: over any field, with pole and residue extraction; `borel_pade_poles` finds the Borel singularities of a series
- **Summation**: Borel–Padé–Laplace sums (`borel_sum`, `laplace_sum`), lateral sums on tilted rays, Stokes discontinuities, median summation; a guard detects poles on the integration ray
- **Acceleration and conformal maps**: Richardson/Wynn extrapolation, Borel-plane conformal maps (`conformal_borel`) for series with branch cuts
- **Large-order analysis**: `large_order_fit` recovers $(A, b, S)$ from $a_n \sim S\,\Gamma(n+b)/A^n$, including multi-saddle fits with complex-conjugate action pairs
- **Hyperasymptotics**: optimal truncation, Dingle terminants, hyperasymptotic summation
- **Transseries and alien calculus**: one-parameter transseries $\sum_n \sigma^n e^{-nA/\hbar}\,\Phi_n$, alien derivatives via the bridge equation, the Stokes automorphism, and Stokes constants from Borel-plane residues or from large-order growth
- **Named examples with exact recursions**: `FormalSeries(:euler, n)`, `:airy`, `:airy_bi`, `:quartic` (Bender–Wu ground-state energy), which also serve as test oracles
- **Optional extensions**: Borel-plane plotting via Makie; AAA rational approximation via BaryRational

## Installation

The package is not registered in the Julia General registry. Install it directly from the repository:

```julia
pkg> add https://github.com/benedikt-nagler/Resurgence.jl
```

## Quick start

```julia
using Resurgence

# The Euler series, exact rationals; divergent for every ħ ≠ 0
Φ = FormalSeries(:euler, 40)

# Borel transform, exactly −1/(1+ζ)
B = borel(Φ)
poles(pade(B; order = 1))         # ζ = −1, the Borel singularity

# Borel–Padé–Laplace summation; precision follows the argument type
borel_sum(Φ, 0.1; order = 1)      # ≈ −e^{10}E₁(10)
borel_sum(Φ, big"0.1"; order = 1) # same pipeline in BigFloat

# Large-order analysis: aₙ ~ S·Γ(n+b)/Aⁿ
fit = large_order_fit(Φ)           # (A, b, S) = (1, 0, 1) for Euler
fit = large_order_fit(FormalSeries(:quartic, 80); order = 8)
                                   # A = 1/3, the quartic instanton action

# Stokes phenomenon: lateral sums, discontinuity, median summation
Ψ = FormalSeries(:airy_bi, 40)     # Borel pole at +4/3, so ℝ₊ is a Stokes ray
lateral_sum(borel(Ψ), 0.3; side = :plus)
median_sum(Ψ, 0.3)

# Transseries and alien calculus
F = Transseries(:euler, 20)
alien_derivative(F; stokes = -2π * im)
```

## Optional extensions

Load additional packages to unlock extensions:

```julia
using CairoMakie                   # or GLMakie
plot_borel_plane(borel(Φ))         # Borel-plane pole scatter with Laplace rays

using BaryRational
aaa_borel(borel(Φ))                # AAA rational approximation of the Borel function
```

## References

- F. J. Dyson, *Divergence of perturbation theory in quantum electrodynamics*, Phys. Rev. **85** (1952), 631–632.
- C. M. Bender, T. T. Wu, *Anharmonic oscillator*, Phys. Rev. **184** (1969), 1231–1260.
- J. Écalle, *Les fonctions résurgentes* I–III, Publ. Math. d'Orsay (1981–1985).
- A. Voros, *The return of the quartic oscillator: the complex WKB method*, Ann. Inst. H. Poincaré A **39** (1983), 211–338.
- M. V. Berry, C. J. Howls, *Hyperasymptotics for integrals with saddles*, Proc. R. Soc. Lond. A **434** (1991), 657–675.
- I. Aniceto, G. Başar, R. Schiappa, *A primer on resurgent transseries and their asymptotics*, Phys. Rept. **809** (2019), 1–135.
- D. Dorigoni, *An introduction to resurgence, trans-series and alien calculus*, Ann. Phys. **409** (2019), 167914.
