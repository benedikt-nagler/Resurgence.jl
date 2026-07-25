```@meta
CurrentModule = Resurgence
```

# Borel plane

[`borel`](@ref) and [`inverse_borel`](@ref) are exact, in rising-factorial normalization.
What you do next is approximate the Borel function from its Taylor coefficients: Padé is the
default, [`ConformalPade`](@ref) handles a cut rather than isolated poles, and `aaa_borel`
(see [Plotting](@ref) for the extension) is the AAA alternative. All of them are
[`AbstractBorelApproximant`](@ref)s and plug into the same summation routines.

```@autodocs
Modules = [Resurgence]
Pages = ["borel.jl"]
```

## Padé approximants

Padé over any field, with pole and residue extraction. [`borel_pade_poles`](@ref) locates the
Borel singularities of a series, which is the practical route to the instanton actions.

```@autodocs
Modules = [Resurgence]
Pages = ["pade.jl"]
```

## Conformal maps

For series whose Borel singularities form a cut rather than isolated poles, where Padé
struggles.

```@autodocs
Modules = [Resurgence]
Pages = ["conformal.jl"]
```
