```@meta
CurrentModule = Resurgence
```

# Large-order analysis

[`large_order_fit`](@ref) recovers ``(A, b, S)`` from ``a_n \sim S\,\Gamma(n+b)/A^n``,
including multi-saddle fits with complex-conjugate action pairs.

```@autodocs
Modules = [Resurgence]
Pages = ["large_order.jl"]
```

## Peeling off a known singularity

[`subtract_singularity`](@ref) removes a known Borel singularity to expose the subleading one.

```@autodocs
Modules = [Resurgence]
Pages = ["singularities.jl"]
```
