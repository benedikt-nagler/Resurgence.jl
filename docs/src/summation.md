```@meta
CurrentModule = Resurgence
```

# Summation

[`borel_sum`](@ref) and [`laplace_sum`](@ref) give Borel–Padé–Laplace sums,
[`lateral_sum`](@ref) works on tilted rays, [`stokes_discontinuity`](@ref) is the jump across
a Stokes line, and [`median_sum`](@ref) is the real average. A guard detects poles sitting on
the integration ray instead of silently returning nonsense — see [`PoleOnRay`](@ref).

```@autodocs
Modules = [Resurgence]
Pages = ["laplace.jl"]
```

## Sequence acceleration

[`accelerate`](@ref) offers Shanks, Wynn-ε, Levin-u/t and Richardson. These act on any
sequence, not only on a Borel pipeline — [`partial_sums`](@ref) turns a series into one — so
they are usable on their own as convergence accelerators.

```@autodocs
Modules = [Resurgence]
Pages = ["acceleration.jl"]
```

## Hyperasymptotics

Optimal truncation, Dingle terminants, and level-one Berry–Howls hyperasymptotic summation.

```@autodocs
Modules = [Resurgence]
Pages = ["hyperasymptotics.jl"]
```
