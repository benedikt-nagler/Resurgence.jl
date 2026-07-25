```@meta
CurrentModule = Resurgence
```

# Plotting

The core depends only on QuadGK, PrecompileTools and the LinearAlgebra stdlib. Plotting and
AAA rational approximation live in package extensions, loaded on demand - calling one of these
before its trigger package is loaded raises an error telling you what to load.

```julia
using CairoMakie                   # or GLMakie
plot_borel_plane(borel(Φ))         # pole scatter with Laplace rays
plot_large_order(Φ)
plot_optimal_truncation(Φ, 0.1)

using BaryRational
aaa_borel(borel(Φ))                # AAA rational approximation of the Borel function
```

## Makie extension

```@docs
plot_borel_plane
plot_large_order
plot_optimal_truncation
```

## BaryRational extension

AAA runs in `Float64` - treat it as a pole and branch-cut hunter, not a high-precision
summation tool (that is Padé's job, in `BigFloat`).

```@docs
aaa_approximant
aaa_borel
```
