@testset "optimal truncation" begin
    setprecision(160) do
        phi = FormalSeries(:euler,40)
        x = big"0.1"
        exact = -exp(1/x)*first(quadgk(t -> exp(-t)/t,1/x,BigFloat(Inf);rtol=big"1e-35"))
        for rule in (:smallest_term,:action)
            result = optimal_truncation(phi,x;rule)
            @test 8 <= result.N <= 11
            # This bound follows from the alternating Euler integral in this example only.
            @test abs(result.value-exact) <= 2result.error
        end
    end
    sparse = FormalSeries([1.0,0.0,-10.0,0.0,200.0,0.0,-6000.0,0.0,240000.0])
    result = optimal_truncation(sparse,0.3)
    @test result.N == 2 && result.value == 1.0
    @test result.error ≈ 0.9
    @test optimal_truncation(FormalSeries([1e16,1.0,1e8]),1.0).error == 1.0
    for phi in (FormalSeries([0.0,0.0]),FormalSeries([1.0,0.0,0.001]))
        result = optimal_truncation(phi,0.1)
        @test result.N == n_terms(phi) && result.error === nothing
    end
    @test_throws InvalidArgument optimal_truncation(FormalSeries([1.0,1e-300]),1e-100)
    @test_throws InvalidArgument optimal_truncation(sparse,0.3;rule=:unknown)
end
