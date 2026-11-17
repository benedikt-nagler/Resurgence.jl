# Archived, experimental, not validated

Snapshot before the 2026-10-26 reduction. `original/` preserves source, tests,
examples, public documentation and ignored planning/instruction files.
Head, status, local patch and SHA-256 digests record provenance.
Notebooks were copied without inspection. No archived file enters package loading,
routine tests or documentation builds. This snapshot is not a validated release.

Removed: transseries, logarithmic sectors, resonance, model-sector solvers, Painleve,
topological strings, alien operators, hyperasymptotic re-expansion and multi-saddle fitting.
`optimal_truncation` and single-series `median_sum` were separated before removal.
Generic conformal maps, acceleration, Darboux subtraction and AAA remain numerical tools.

## Reproduced faults

The operator `stokes_automorphism(MultiTransseries(Transseries(:euler,12));stokes=-2pi*im)`
evaluates at sigma=0, hbar=0.2 to -6.283185307i. Parameter translation instead gives
0.270766255491-0.021167884793i. The operator omits the exponential grading, and addition
loses the longer perturbative sector. An old test explicitly enshrined the incompatible
operator effect. The archived operator is not repaired.

The original `optimal_truncation(conifold_gap_series(12),0.1)` reports N=1 and error=0.
Structural zeros and subtraction of partial sums cause that false diagnostic.
The active generic routine is repaired without retaining the conifold application.

No complete verification of the removed model-dependent mathematics was done.
References and passing tests alone do not validate physical interpretations.
