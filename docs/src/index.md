# Resurgence.jl

Analyze supplied asymptotic coefficients with formal arithmetic, explicit Borel
normalization, rational approximation and numerical Laplace integration.

The reduced Borel coefficients are `b_n = a_n/(beta)_n`. The Laplace integral restores
`zeta^(beta-1)/Gamma(beta)`. With leading power zero the constant is separated first.
Fractional powers use the principal branch. Unsupported winding around its cut is refused.

Padé, conformal and AAA approximants infer analytic structure from finite data.
Their poles need not be singularities of the original function. Quadrature error estimates
concern the approximant only. A lateral deformation must not cross its known singularities.
Use an explicit regular Laplace ray when choosing a different contour for an approximant.
`median_sum` is an arithmetic average of two single-series laterals.

`large_order_fit` assumes a single factorial growth law with nonzero coefficients.
It does not establish a complete action lattice. `optimal_truncation.error` is the
first omitted nonzero term's magnitude, not a remainder bound. It is `nothing` when
the stored range does not resolve a minimum. Structural zeros cannot set that minimum.

Linear and quadratic Padé poles use BigFloat formulas. Higher degrees still use Float64
companion seeds followed by optional Newton refinement. Poorly conditioned or clustered
higher-degree roots need independent checks. Small residuals do not certify root locations.
Angles supplied as Float64 retain that input precision even in a BigFloat calculation.

Named Euler, Airy, exponential, quartic and zero-dimensional phi4 series are coefficient
reference cases only. Removed models are preserved in `archive/2026-10-26-reduction/`,
explicitly experimental and unvalidated. They are excluded from loading, tests and this manual.
