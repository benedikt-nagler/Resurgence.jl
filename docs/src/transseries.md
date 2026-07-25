```@meta
CurrentModule = Resurgence
```

# Transseries

One-parameter [`Transseries`](@ref) ``\sum_n \sigma^n e^{-nA/\hbar}\Phi_n``, and
[`MultiTransseries`](@ref) for several actions at once, stored sparsely over the action
lattice.

The one- and many-action types share most of their interface — the accessors below dispatch on
both.

```@autodocs
Modules = [Resurgence]
Pages = ["transseries.jl", "multi_transseries.jl"]
```

## Resonance

When the actions satisfy integer relations the transseries develops logarithms.
[`resonance_lattice`](@ref) computes the relations exactly over ``\mathbb{Z}``, and
[`resonant_solve`](@ref) is the primitive that generates the log tower.

```@autodocs
Modules = [Resurgence]
Pages = ["resonance.jl"]
```

## Solving ODEs

[`transseries_solve`](@ref) is a model-agnostic driver that builds a transseries solution
sector by sector from a user-supplied linear solve.

```@autodocs
Modules = [Resurgence]
Pages = ["ode.jl"]
```

## Painlevé I

[`painleve1`](@ref) applies the driver to Painlevé I, with both normalization conventions of
the instanton action available (`:gikm` by default, `:string` for the string-equation
convention).

```@autodocs
Modules = [Resurgence]
Pages = ["painleve.jl"]
```
