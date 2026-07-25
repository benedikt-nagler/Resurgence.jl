```@meta
CurrentModule = Resurgence
```

# Series

[`FormalSeries`](@ref) is a truncated power series over any coefficient type, with
arithmetic, truncation, differentiation and fractional power offsets. [`LogSeries`](@ref)
adds a polynomial in ``\log\hbar``, which is what resonant problems produce.

[`FormalSeries`](@ref)`(:name, n)` builds a named example from its exact recursion. These
double as test oracles - the suite checks mathematical identities rather than stored numbers.

```@autodocs
Modules = [Resurgence]
Pages = ["formal_series.jl"]
```

## Series with logarithms

```@autodocs
Modules = [Resurgence]
Pages = ["log_series.jl"]
```

## Named examples

```@autodocs
Modules = [Resurgence]
Pages = ["named_series.jl"]
```
