using Documenter
using Resurgence

DocMeta.setdocmeta!(Resurgence, :DocTestSetup, :(using Resurgence); recursive = true)

makedocs(;
    modules = [Resurgence],
    authors = "Benedikt Nagler <benedikt.nagler@protonmail.com>",
    sitename = "Resurgence.jl",
    format = Documenter.HTML(;
        canonical = "https://benedikt-nagler.github.io/Resurgence.jl",
        edit_link = "main",
        assets = String[],
        mathengine = Documenter.KaTeX(),
    ),
    pages = [
        "Home" => "index.md",
        "Manual" => [
            "Series" => "series.md",
            "Borel plane" => "borel.md",
            "Summation" => "summation.md",
            "Large-order analysis" => "large_order.md",
            "Transseries" => "transseries.md",
            "Alien calculus" => "alien.md",
            "Plotting" => "plotting.md",
        ],
        "Errors" => "errors.md",
        "API index" => "api.md",
    ],
    checkdocs = :exports,
)

deploydocs(; repo = "github.com/benedikt-nagler/Resurgence.jl", devbranch = "main")
