```@meta
CurrentModule = Resurgence
```

# Alien calculus

[`alien_derivative`](@ref) via the bridge equation (single- and multi-parameter, by lattice
charge or by singularity location), [`stokes_automorphism`](@ref) as both a ``\sigma``-shift
and an operator, and [`stokes_constant`](@ref) from either Borel-plane residues or large-order
growth.

The single- and multi-action routines share their interface — the functions below dispatch on
both [`Transseries`](@ref) and [`MultiTransseries`](@ref).

```@autodocs
Modules = [Resurgence]
Pages = ["alien.jl", "multi_alien.jl"]
```
