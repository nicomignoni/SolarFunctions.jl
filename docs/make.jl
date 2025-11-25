# push!(LOAD_PATH, "../src/")

using Documenter, Irradia

makedocs(
    sitename="Irradia.jl",
    pages = [
        "Home" => "index.md",
        "API" => "api.md"
    ],
    repo=Remotes.GitHub("nicomignoni", "Irradia.jl")
)
