push!(LOAD_PATH, "../src/")

using Documenter, Irradia

makedocs(
    sitename="Irradia.jl"
    pages = [
        "API" => "api.md"
    ]
)
