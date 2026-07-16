using Resurgence
using Test
using Aqua
import AbstractAlgebra
import SpecialFunctions
using QuadGK: quadgk

@testset "Resurgence.jl" begin
    @testset "Aqua" begin
        Aqua.test_all(Resurgence; ambiguities = false)
        Aqua.test_ambiguities(Resurgence)
    end
    include("test_errors.jl")
    include("test_formal_series.jl")
    include("test_log_series.jl")
    include("test_named_series.jl")
    include("test_borel.jl")
    include("test_pade.jl")
    include("test_laplace.jl")
    include("test_large_order.jl")
    include("test_acceleration.jl")
    include("test_conformal.jl")
    include("test_singularities.jl")
    include("test_hyperasymptotics.jl")
    include("test_transseries.jl")
    include("test_multi_transseries.jl")
    include("test_resonance.jl")
    include("test_alien.jl")
    include("test_multi_alien.jl")
    include("test_show.jl")
    include("test_makie_ext.jl")
    include("test_aaa_ext.jl")
end
