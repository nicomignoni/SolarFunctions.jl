using Documenter, DocumenterCitations, SolarFunctions

bib = CitationBibliography(joinpath(@__DIR__, "src", "references.bib"))

makedocs(;
    sitename="SolarFunctions.jl",
    pages = [
        "Home" => "index.md",
        "API" => "api.md",
        "References" => "references.md"
    ],
    format = Documenter.HTML(
        edit_link="master",
        assets=["assets/favicon.ico"]
    ),
    repo=Remotes.GitHub("nicomignoni", "SolarFunctions.jl"),
    plugins=[bib]
)

deploydocs(
    repo = "github.com/nicomignoni/SolarFunctions.jl.git",
)
