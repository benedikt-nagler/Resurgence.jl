# Makie weakdep extension: Borel-plane pole scatter with Laplace-ray overlay.
module ResurgenceMakieExt

using Makie: Makie, Figure, Axis, scatter!, lines!, axislegend
using Resurgence: Resurgence, BorelSeries, PadeApproximant, pade, poles

"""
    plot_borel_plane(B::BorelSeries; order = nothing, rays = [0.0], refine = true)
    plot_borel_plane(r::PadeApproximant; rays = [0.0], refine = true)

Scatter the Borel–Padé poles in the complex ζ-plane, with the Laplace rays at the
angles in `rays` overlaid. Returns a `Makie.Figure`.
"""
function Resurgence.plot_borel_plane(B::BorelSeries;
                                     order::Union{Nothing,Integer} = nothing,
                                     kwargs...)
    Resurgence.plot_borel_plane(pade(B; order); kwargs...)
end

function Resurgence.plot_borel_plane(r::PadeApproximant;
                                     rays::AbstractVector{<:Real} = [0.0],
                                     refine::Bool = true)
    ζs = ComplexF64.(poles(r; refine))
    fig = Figure()
    ax = Axis(fig[1, 1]; xlabel = "Re ζ", ylabel = "Im ζ",
              title = "Borel plane: [$(Resurgence.numerator_degree(r))/" *
                      "$(Resurgence.denominator_degree(r))] Padé poles")
    R = isempty(ζs) ? 1.0 : 1.5 * maximum(abs, ζs)
    for θ in rays
        lines!(ax, [0.0, R * cos(θ)], [0.0, R * sin(θ)];
               color = :gray, linestyle = :dash)
    end
    isempty(ζs) ||
        scatter!(ax, real.(ζs), imag.(ζs); color = :crimson, label = "Padé poles")
    isempty(ζs) || axislegend(ax)
    fig
end

end # module ResurgenceMakieExt
