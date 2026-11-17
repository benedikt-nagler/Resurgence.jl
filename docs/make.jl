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
    pages = ["Home" => "index.md", "Series and summation" => "api.md", "Numerical helpers" => "numerics.md"],
    checkdocs = :exports,
)

deploydocs(; repo = "github.com/benedikt-nagler/Resurgence.jl", devbranch = "main")
