@testset "errors" begin
    for E in (InvalidArgument, IncompatibleSeries, DegeneratePade, PoleOnRay,
              NoSingularityFound)
        @test E <: ResurgenceError
    end
    @test ResurgenceError <: Exception
    # every showerror message starts with the error-type name
    @test startswith(sprint(showerror, InvalidArgument("x")), "InvalidArgument")
    @test startswith(sprint(showerror, IncompatibleSeries(:var, :ħ, :x)),
                     "IncompatibleSeries")
    @test startswith(sprint(showerror, DegeneratePade(2, 2)), "DegeneratePade")
    @test startswith(sprint(showerror, PoleOnRay(1.0 + 0im, 0.0)), "PoleOnRay")
    @test startswith(sprint(showerror, NoSingularityFound("x")), "NoSingularityFound")
end
