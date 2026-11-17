# Resurgence.jl

Analyze supplied asymptotic series with explicit Borel normalization, Padé approximation
and numerical Laplace summation. This development tree reduces the unreleased v0.2 scope.

## Installation

```julia
using Pkg
Pkg.add("Resurgence")
```

The registered release may contain APIs archived in this working tree.
Use `Pkg.develop(path="/path/to/Resurgence.jl")` for the reduced development version.

```julia
using Resurgence
phi = FormalSeries(:euler, 12)
B = borel(phi)
r = pade(B; order=1)
laplace_sum(B, r, 0.2)
```

- **Formal series:** truncated arithmetic with explicit powers and coefficient types.
- **Borel normalization:** rising factorials and a separately stored leading constant.
- **Approximants:** Padé, supplied conformal maps and optional BaryRational sampling.
- **Summation:** numerical Laplace rays, laterals and the arithmetic average of two laterals.
- **Large order:** a single factorial-growth fit with stated assumptions.
- **Truncation:** direct nonzero term magnitudes and explicit unknown remainder diagnostics.
- **References:** Euler, Airy, quartic ground-state and zero-dimensional phi4 coefficients.

Integration of an approximant does not certify the original Borel sum. Higher-degree pole
locations retain a Float64 seeding limitation. Term sizes and order changes are not error bounds.

See the [manual](https://benedikt-nagler.github.io/Resurgence.jl/).
Removed layers are [archived](archive/2026-10-26-reduction/STATUS.md) and unvalidated.
ExactWKB.jl consumes the series-analysis core. See [LICENSE](LICENSE).
