using BayesianDiscovery
using Documenter

DocMeta.setdocmeta!(BayesianDiscovery, :DocTestSetup, :(using BayesianDiscovery); recursive=true)

makedocs(;
    modules=[BayesianDiscovery],
    authors="Josh North <jsnowynorth@gmail.com>",
    repo="https://github.com/jsnowynorth/BayesianDiscovery.jl/blob/{commit}{path}#{line}",
    sitename="BayesianDiscovery.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)
