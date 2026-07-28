# Changelog

All notable changes to this project are documented in this file. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the project adheres to
[semantic versioning](https://semver.org/spec/v2.0.0.html).

## [0.1.1] - 2026-07-28

Documentation-only release. No changes to the public API.

### Changed

- The README and the documentation home page now give `pkg> add Resurgence` as the
  installation instruction, since the package is in the General Registry.
- The README carries a registry version badge.

### Added

- This changelog.

## [0.1.0] - 2026-07-25

First release, registered on the Julia General Registry.

### Added

- **Formal series.** `FormalSeries{T}` over any coefficient type, with a rational
  `power_offset` for fractional sectors, exact `Rational` and `BigFloat` coefficients alike,
  and a small library of named oracles (`FormalSeries(:euler, n)`, `:airy`, `:airy_bi`,
  `:quartic`, `:exp`) used as worked examples and test oracles rather than as a library of
  expansions.
- **Borel transform and Padé.** Exact `borel` / `inverse_borel` in the rising-factorial
  normalization; `pade` over any field with `poles`, `residues` and `borel_pade_poles`, plus
  opt-in degeneracy reduction; AAA rational approximation (`aaa_approximant`, `aaa_borel`)
  through a `BaryRational` package extension.
- **Summation.** `laplace_sum`, `lateral_sum`, `stokes_discontinuity` and `borel_sum` over a
  common `AbstractBorelApproximant` seam, with precision following the caller's argument
  types. Sequence acceleration (`accelerate`: Shanks, Wynn-epsilon, Levin-u and Levin-t;
  `partial_sums`), conformal maps (`ConformalMap`, `ConformalPade`, `conformal_borel`), and
  superasymptotics plus level-1 hyperasymptotics (`dingle_terminant`, `optimal_truncation`,
  `hyper_sum`).
- **Large-order analysis.** `richardson` and `large_order_fit` recovering the action,
  exponent and Stokes constant from coefficient growth, including multi-saddle fits with
  oscillatory conjugate pairs, and Darboux peeling via `subtract_singularity`.
- **Transseries and alien calculus.** One-parameter `Transseries{T,A}` with the bridge
  equation (`alien_derivative`), the Stokes automorphism as a parameter shift
  (`stokes_automorphism`, `transseries_sum`), median summation (`median_sum`) and two routes
  to `stokes_constant`. Multi-parameter `MultiTransseries{T,A,K}` over an action lattice with
  `pointed_alien_derivative`, and resonant sectors carrying logarithms via `LogSeries{T}`
  (`resonance_lattice`, `resonance_depth`, `resonant_solve`).
- **Nonlinear ODEs.** `transseries_solve` for formal transseries solutions, with Painleve I
  (`painleve1`, `painleve1_action`) as the worked case.
- **Plotting** through a `Makie` package extension: `plot_borel_plane`, `plot_large_order`,
  `plot_optimal_truncation`.

[0.1.1]: https://github.com/benedikt-nagler/Resurgence.jl/releases/tag/v0.1.1
[0.1.0]: https://github.com/benedikt-nagler/Resurgence.jl/releases/tag/v0.1.0
