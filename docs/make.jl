using Documenter, MPAlign

makedocs(
    modules = [MPAlign],
    format = Documenter.HTML(; prettyurls = get(ENV, "CI", nothing) == "true"),
    authors = "Marc Lelarge",
    sitename = "MPAlign.jl",
    pages = Any["index.md"]
    # strict = true,
    # clean = true,
    # checkdocs = :exports,
)

deploydocs(
    repo = "github.com/mlelarge/MPAlign.jl.git",
    push_preview = true
)
