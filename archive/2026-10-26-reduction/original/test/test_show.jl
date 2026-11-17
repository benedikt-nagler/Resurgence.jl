@testset "show" begin
    Φ = FormalSeries([1//1, -2//1, 0//1, 3//1])
    s = sprint(show, Φ)
    @test startswith(s, "FormalSeries{Rational{Int64}}")
    @test occursin("O(ħ^4)", s)
    @test occursin("- 2", s)                       # sign folded into the term
    @test !occursin("+ -", s)
    # zero coefficients are skipped, at most a few terms shown
    @test !occursin("0//1*", s)
    long = sprint(show, FormalSeries(:euler, 100))
    @test length(long) < 200

    latex = sprint(show, MIME"text/latex"(), FormalSeries([1//2, 1//1]; power_offset = 1//2))
    @test occursin("\\frac{1}{2}", latex)
    @test occursin("ħ^{1/2}", latex)
    # LaTeX joins coefficient and power by juxtaposition, never a literal *
    latex2 = sprint(show, MIME"text/latex"(), FormalSeries(:euler, 5))
    @test !occursin("*", latex2)
    @test occursin("2\\,ħ^{3}", latex2)
    @test occursin("*", sprint(show, FormalSeries(:euler, 5)))   # text mode keeps it

    B = borel(FormalSeries(:euler, 5))
    sb = sprint(show, B)
    @test startswith(sb, "BorelSeries")
    @test occursin("β = 1", sb) && occursin("from ħ", sb) && occursin("ζ", sb)
    @test occursin("\\Gamma", sprint(show, MIME"text/latex"(), B))
    # split-off constant is displayed
    @test occursin("[const 1", sprint(show, borel(FormalSeries(:airy, 5))))

    r = pade(B.series.coeffs, 1, 1)
    @test occursin("[1/1]", sprint(show, r))

    m = conformal_map(-1 // 1)
    sm = sprint(show, m)
    @test startswith(sm, "ConformalMap") && occursin("ζ₀ = -1", sm)
    @test occursin("pair", sprint(show, conformal_map(2im; pair = true)))
    c = conformal_borel(B; zeta0 = -1 // 1)
    sc = sprint(show, c)
    @test startswith(sc, "ConformalPade") && occursin("ζ₀ = -1", sc)

    F = Transseries(:euler, 5)
    sf = sprint(show, F)
    @test startswith(sf, "Transseries")
    @test occursin("A = 1", sf) && occursin("2 sectors", sf) && occursin("Φ₀", sf)
    lf = sprint(show, MIME"text/latex"(), F)
    @test occursin("\\sigma^n", lf) && occursin("\\Phi_n", lf)
    # rational action rendered as a fraction
    @test occursin("\\frac{-4}{3}",
                   sprint(show, MIME"text/latex"(), Transseries(:airy, 3)))
end

@testset "show: log sectors" begin
    logh = LogSeries([FormalSeries([0//1]), FormalSeries([1//1])])
    s = sprint(show, logh)
    @test startswith(s, "LogSeries")
    @test occursin("log-degree 1", s) && occursin("log(ħ)", s)
    l = sprint(show, MIME"text/latex"(), logh)
    @test occursin("\\log ħ", l)
    # degree ≥ 2 renders the power
    L2 = LogSeries([FormalSeries([1//1]), FormalSeries([0//1]), FormalSeries([2//1])])
    @test occursin("log^2(ħ)", sprint(show, L2))
    @test occursin("\\log^{2} ħ", sprint(show, MIME"text/latex"(), L2))
    # a resonant MultiTransseries advertises its log degree; a non-resonant one does not
    R = MultiTransseries((1//1, -1//1), Dict((0, 0) => logh))
    sr = sprint(show, R)
    @test startswith(sr, "MultiTransseries") && occursin("log-degree ≤ 1", sr)
    P = MultiTransseries((1//1, -1//1), Dict((0, 0) => FormalSeries([1//1])))
    @test !occursin("log-degree", sprint(show, P))
    @test occursin("\\Phi_{(0, 0)}", sprint(show, MIME"text/latex"(), R))
end
