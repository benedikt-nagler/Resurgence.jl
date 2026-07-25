```@meta
CurrentModule = Resurgence
```

# Errors

Every invalid operation throws a typed error from the [`ResurgenceError`](@ref) hierarchy,
each with an informative `showerror`. Structs are immutable and every operation returns a new
object, so an error never leaves a half-mutated value behind.

```@autodocs
Modules = [Resurgence]
Pages = ["errors.jl"]
```
