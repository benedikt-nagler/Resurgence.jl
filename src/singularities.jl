# Borel-plane singularity tools: Darboux peeling (subtract the dominant singularity
# so the subleading one becomes visible) and the multi-saddle large-order fit
# a_n ≈ Σ_j S_j Γ(p_n + b_j) A_j^{-p_n} for several competing actions.

# Γ(x) for real x (poles at nonpositive integers excluded), via _gamma and the
# recurrence for x ≤ 0.
function _gamma_real(x::Real)
    x > 0 && return _gamma(x)
    (isinteger(x) || iszero(x)) &&
        throw(InvalidArgument("Γ has a pole at $x"))
    m = Int(ceil(-float(x))) + 1
    g = _gamma(float(x) + m)
    z = float(x) + m
    for _ in 1:m
        z -= 1
        g /= z
    end
    g
end

"""
    subtract_singularity(B::BorelSeries, zeta0; beta = 1, c = nothing,
                         order = nothing, method = :wynn, accel_order = nothing)
        -> NamedTuple{(:peeled, :c, :zeta0, :beta)}
    subtract_singularity(B::BorelSeries; kwargs...)

Darboux peeling: model the dominant singularity of the *reduced* Borel function as
``c\\,(1 - ζ/ζ_0)^{-β}`` (Taylor coefficients ``s_n = c\\,(β)_n ζ_0^{-n}/n!``),
subtract it coefficient-by-coefficient, and return the peeled
[`BorelSeries`](@ref) (same `B.beta`, constant term, and source variable) — so
`poles(pade(out.peeled))` or `large_order_fit(inverse_borel(out.peeled))` can
resolve the **subleading** singularity/action that the dominant one masked.

Here `beta` is the singularity exponent ``β`` — **not** the `B.beta` of the Borel
normalization; `beta = 1` is a simple pole. The amplitude `c` is estimated when not
given: for `beta = 1` from the Padé residue nearest `zeta0`
(``c = -\\mathrm{res}/ζ_0``, exact for a rational Borel function); for general
`beta` by [`accelerate`](@ref)ing the sequence ``b_n ζ_0^n\\, n!/(β)_n → c``
(kwargs `method`, `accel_order`). Without an explicit `zeta0` the
smallest-modulus Borel–Padé pole is used (`order` is passed to [`pade`](@ref));
throws [`NoSingularityFound`](@ref) when there is none.
"""
function subtract_singularity(B::BorelSeries, zeta0::Number;
                              beta::Real = 1,
                              c::Union{Nothing,Number} = nothing,
                              order::Union{Nothing,Integer} = nothing,
                              method::Symbol = :wynn,
                              accel_order::Union{Nothing,Integer} = nothing)
    iszero(zeta0) && throw(InvalidArgument("the singularity position must be ≠ 0"))
    beta > 0 || throw(InvalidArgument("the singularity exponent beta must be > 0"))
    b = B.series.coeffs
    N = length(b)
    camp = c
    if camp === nothing
        if beta == 1
            rs = residues(pade(B; order, reduce = true))
            isempty(rs) &&
                throw(NoSingularityFound("subtract_singularity: the Borel–Padé " *
                                         "approximant has no poles; pass c"))
            _, res = argmin(p -> abs(first(p) - zeta0), rs)
            camp = -res / zeta0
        else
            # c_n = b_n ζ₀ⁿ n!/(β)_n → c, built by the recurrence
            # f_{n+1} = f_n ζ₀ (n+1)/(β+n) on the coefficient-free factor
            f = one(float(b[1])) * one(float(zeta0))
            cn = Vector{typeof(float(b[1]) * f)}(undef, N)
            for n in 0:(N - 1)
                cn[n + 1] = b[n + 1] * f
                f *= zeta0 * (n + 1) / (beta + n)
            end
            camp = accelerate(cn; method, order = accel_order)
        end
    end
    # peel: t_n = c (β)_n ζ₀^{-n}/n!, with t_{n+1} = t_n (β+n)/((n+1) ζ₀)
    t = camp * one(b[1])
    peeled = Vector{typeof(b[1] - t)}(undef, N)
    for n in 0:(N - 1)
        peeled[n + 1] = b[n + 1] - t
        t *= (beta + n) / ((n + 1) * zeta0)
    end
    out = BorelSeries(FormalSeries(peeled, variable(B);
                                   power_offset = B.beta - 1),
                      B.beta, oftype(peeled[1], B.constant_term), B.source_var)
    (peeled = out, c = camp, zeta0 = zeta0, beta = beta)
end

function subtract_singularity(B::BorelSeries;
                              order::Union{Nothing,Integer} = nothing, kwargs...)
    ζs = poles(pade(B; order, reduce = true))
    isempty(ζs) &&
        throw(NoSingularityFound("subtract_singularity: the Borel–Padé approximant " *
                                 "has no poles; pass zeta0 explicitly"))
    subtract_singularity(B, ζs[1]; order, kwargs...)
end

# ---- multi-saddle large-order fit --------------------------------------------------
#
# Model: a_n ≈ Σ_g S_g Γ(p_n + b_g) A_g^{-p_n} (+ the conjugate for pair groups),
# p_n = n + power_offset. Actions are initialized from the Borel–Padé poles (which
# *are* the actions), amplitudes S_g solved linearly at fixed (A_g, b_g) (variable
# projection), then (Re A, Im A, b) polished by damped Gauss–Newton with a numeric
# Jacobian. For real coefficient input, complex-conjugate pole pairs are fitted as
# a single group contributing 2 Re[S Γ(p_n+b) A^{-p_n}] — the oscillatory
# large-order case where single-saddle ratio methods fail.

# Model column g_n = Γ(p_n + b) A^{-p_n} over the window (0-based indices `win`).
function _saddle_column(A::Complex{BigFloat}, b::BigFloat, p0::BigFloat,
                        win::UnitRange{Int})
    g = Vector{Complex{BigFloat}}(undef, length(win))
    p1 = p0 + first(win)
    Γ = _gamma_real(p1 + b) * one(Complex{BigFloat})
    Ap = exp(-(p1) * log(A))                    # A^{-p_1}, principal branch
    Ainv = inv(A)
    for (i, n) in pairs(win)
        g[i] = Γ * Ap
        Γ *= (p0 + n + b)                       # Γ(p_{n+1} + b)
        Ap *= Ainv
    end
    g
end

# Real-arithmetic least squares for the amplitudes at fixed groups: unknowns
# [Re S_g, Im S_g]; pair groups contribute 2Re(S g) (real rows only get the pair
# doubling — their imaginary part is identically zero and stays in the residual).
function _solve_amplitudes(cols, pairflags, a::Vector{Complex{BigFloat}})
    Nw = length(a)
    G = length(cols)
    D = Matrix{BigFloat}(undef, 2Nw, 2G)
    for g in 1:G
        for i in 1:Nw
            gr, gi = reim(cols[g][i])
            if pairflags[g]
                D[i, 2g - 1] = 2gr;  D[i, 2g] = -2gi
                D[Nw + i, 2g - 1] = 0;  D[Nw + i, 2g] = 0
            else
                D[i, 2g - 1] = gr;   D[i, 2g] = -gi
                D[Nw + i, 2g - 1] = gi;  D[Nw + i, 2g] = gr
            end
        end
    end
    rhs = vcat(real.(a), imag.(a))
    # normal equations via the self-owned dense solver (windows are small and
    # BigFloat precision absorbs the squared conditioning)
    AtA = D' * D
    Atb = D' * rhs
    x = _solve_dense!(AtA, Atb)
    x === nothing && return nothing
    [Complex{BigFloat}(x[2g - 1], x[2g]) for g in 1:G]
end

# model values over the window for the given groups and amplitudes
function _saddle_model(cols, pairflags, S)
    m = zeros(Complex{BigFloat}, length(cols[1]))
    for g in eachindex(cols)
        if pairflags[g]
            m .+= 2 .* real.(S[g] .* cols[g])
        else
            m .+= S[g] .* cols[g]
        end
    end
    m
end

function _multi_saddle_fit(Φ::FormalSeries, J::Int;
                           order::Union{Nothing,Integer},
                           window::Union{Nothing,AbstractRange},
                           maxiter::Integer)
    N = n_terms(Φ)
    a = [Complex{BigFloat}(float(real(c)), float(imag(c))) for c in Φ.coeffs]
    any(iszero, a) &&
        throw(InvalidArgument("large_order_fit requires nonzero coefficients"))
    p0 = BigFloat(power_offset(Φ))
    isreal_input = all(iszero ∘ imag, a)

    # 1) initialize the actions from the Borel–Padé poles
    ζs = poles(pade(borel(Φ); order, reduce = true))
    length(ζs) ≥ J ||
        throw(NoSingularityFound("multi-saddle fit with saddles = $J needs ≥ $J " *
                                 "Borel–Padé poles, found $(length(ζs))"))
    A0 = ζs[1:J]

    # group conjugate pairs for real input
    tol = big"1e-6"
    used = falses(J)
    groups = NamedTuple{(:A, :b, :pair),Tuple{Complex{BigFloat},BigFloat,Bool}}[]
    for i in 1:J
        used[i] && continue
        paired = false
        if isreal_input && abs(imag(A0[i])) > tol * abs(A0[i])
            for j in (i + 1):J
                used[j] && continue
                if abs(A0[j] - conj(A0[i])) < tol * abs(A0[i])
                    used[j] = true
                    paired = true
                    break
                end
            end
        end
        used[i] = true
        Ai = paired && imag(A0[i]) < 0 ? conj(A0[i]) : A0[i]
        push!(groups, (A = Ai, b = zero(BigFloat), pair = paired))
    end
    G = length(groups)

    win = window === nothing ? ((N ÷ 2):(N - 1)) : (Int(first(window)):Int(last(window)))
    (first(win) ≥ 0 && last(win) ≤ N - 1 && length(win) ≥ 2G + 1) ||
        throw(InvalidArgument("fit window $win invalid for $N coefficients and " *
                              "$G saddle groups (needs ≥ $(2G + 1) points in 0:$(N - 1))"))
    aw = a[(first(win) + 1):(last(win) + 1)]
    scale = abs.(aw)

    # residual (weighted, stacked re/im) and amplitudes for a parameter vector
    # θ = [ReA_1, ImA_1, b_1, ReA_2, …]
    function eval_theta(θ::Vector{BigFloat})
        cols = [_saddle_column(Complex{BigFloat}(θ[3g - 2], θ[3g - 1]), θ[3g], p0, win)
                for g in 1:G]
        S = _solve_amplitudes(cols, [gr.pair for gr in groups], aw)
        S === nothing && return nothing, nothing
        m = _saddle_model(cols, [gr.pair for gr in groups], S)
        r = vcat(real.(m .- aw) ./ scale, imag.(m .- aw) ./ scale)
        r, S
    end

    θ = BigFloat[]
    for gr in groups
        push!(θ, real(gr.A), imag(gr.A), gr.b)
    end
    r, S = eval_theta(θ)
    r === nothing && throw(InvalidArgument("multi-saddle fit: singular amplitude " *
                                           "system at the initial actions"))
    cost = sum(abs2, r)
    λ = big"1e-3"
    converged = false
    for _ in 1:maxiter
        # numeric Jacobian of the residual (S re-solved per perturbation)
        Jm = Matrix{BigFloat}(undef, length(r), length(θ))
        ok = true
        for k in eachindex(θ)
            h = sqrt(eps(BigFloat)) * (1 + abs(θ[k]))
            θp = copy(θ)
            θp[k] += h
            rp, _ = eval_theta(θp)
            rp === nothing && (ok = false; break)
            Jm[:, k] = (rp .- r) ./ h
        end
        ok || break
        JtJ = Jm' * Jm
        Jtr = Jm' * r
        accepted = false
        for _ in 1:12                            # Levenberg damping sweeps
            Aλ = JtJ + λ * LinearAlgebra.Diagonal(max.(LinearAlgebra.diag(JtJ),
                                                       eps(BigFloat)))
            δ = _solve_dense!(Matrix(Aλ), copy(-Jtr))
            δ === nothing && (λ *= 10; continue)
            θn = θ .+ δ
            rn, Sn = eval_theta(θn)
            rn === nothing && (λ *= 10; continue)
            costn = sum(abs2, rn)
            if costn < cost
                θ, r, S, cost = θn, rn, Sn, costn
                λ = max(λ / 10, big"1e-30")
                accepted = true
                if sum(abs2, δ) ≤ eps(BigFloat) * (1 + sum(abs2, θ))
                    converged = true
                end
                break
            else
                λ *= 10
            end
        end
        if !accepted                              # local minimum reached
            converged = cost < big"1e-40" || converged
            break
        end
        converged && break
    end
    converged |= cost < big"1e-40"

    out = NamedTuple{(:A, :b, :S),Tuple{Complex{BigFloat},BigFloat,Complex{BigFloat}}}[]
    for g in 1:G
        A = Complex{BigFloat}(θ[3g - 2], θ[3g - 1])
        push!(out, (A = A, b = θ[3g], S = S[g]))
        groups[g].pair && push!(out, (A = conj(A), b = θ[3g], S = conj(S[g])))
    end
    sort!(out; by = x -> (abs(x.A), angle(x.A)))
    (saddles = out, converged = converged)
end
