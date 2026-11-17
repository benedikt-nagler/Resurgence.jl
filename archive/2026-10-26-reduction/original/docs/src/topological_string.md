```@meta
CurrentModule = Resurgence
```

# Topological strings

A topological-string free energy is the canonical divergent series of the subject, and
two Calabi-Yau geometries have one in closed form. Both are available here as ordinary
series, so the rest of the package applies to them.

`Resurgence.jl` is a **consumer** here. It never computes a Gromov-Witten or
Gopakumar-Vafa invariant, never runs topological recursion, and never runs the
holomorphic-anomaly recursion. Everything below is a quoted closed form. For a geometry
without one, generate the ``F_g`` elsewhere and feed them in as a
[`FormalSeries`](@ref) on the same grid.

The worked example, with figures, is `examples/topological_string_resurgence.ipynb`.

## Background

The topological string on a Calabi-Yau threefold ``X`` has a genus expansion in the
string coupling ``\lambda``,

```math
F(\lambda, t) = \sum_{g \ge 0} \lambda^{2g-2} F_g(t), \qquad
F_g(t) = \sum_{\beta} N_{g,\beta}\, Q^{\beta}, \qquad Q = e^{-t},
```

where ``N_{g,\beta}`` are the genus-``g`` Gromov-Witten invariants and ``t`` is the
complexified Kähler parameter. The sum over ``g`` diverges factorially, and the Borel
singularities of ``F`` in ``\lambda`` sit at ``2\pi i`` times periods of ``X``: the
central charges of D-branes. Reading those actions and their multiplicities off the
``F_g`` is what this package is for.

## Grid convention

Both constructors return a series in ``\lambda`` on the **unit** ``\lambda`` grid: power
offset 2, odd coefficients exactly zero, and the coefficient of ``\lambda^{2g-2}`` equal
to ``F_g``. This is not cosmetic. It makes [`borel`](@ref) divide genus ``g`` by
``\Gamma(2g-2) = (2g-3)!``, which is precisely the normalization in which the gap
series' Borel transform is elementary ([`conifold_gap_borel`](@ref)). Expressed in
``x = \lambda^2`` instead, the same series is Gevrey-2 and nothing closes.

The price is that [`large_order_fit`](@ref) does not apply: it requires nonzero
coefficients and Gevrey-1 growth. Extract the action from the genus subsequence with
[`richardson`](@ref) instead,

```math
\frac{F_g}{F_{g+1}}\,(2g-2)(2g-1) \;\longrightarrow\; A^2 = (2\pi i t)^2 .
```

## Free energies

```@docs
conifold_gap_series
resolved_conifold_series
```

## Borel-plane data

The gap series is the package's first example with **double** Borel poles, which has two
consequences. [`pade`](@ref) renders each as a coalescing pair of simple poles straddling the
true location, so the pair's midpoint is the singularity and the splitting is the diagnostic.
[`stokes_constant`](@ref)'s residue route does not apply, since it is exact for simple poles
only; the exact Laurent data is available analytically from
[`conifold_gap_stokes_data`](@ref) instead.

```@docs
conifold_gap_borel
conifold_gap_actions
conifold_gap_stokes_data
resolved_conifold_actions
```
