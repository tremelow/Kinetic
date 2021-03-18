push!(LOAD_PATH,joinpath(@__DIR__, ".."))
using Documenter, Kinetic

makedocs(
    modules = [StiffKinetic],
    format = Documenter.HTML(; prettyurls = get(ENV, "CI", nothing) == "true"),
    authors = "Léopold Trémant",
    sitename = "StiffKinetic.jl",
    pages = Any["index.md"]
    # strict = true,
    # clean = true,
    # checkdocs = :exports,
)

deploydocs(
    repo = "github.com/tremelow/StiffKinetic.jl.git",
)
