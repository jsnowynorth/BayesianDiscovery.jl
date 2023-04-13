using BayesianDiscovery
using Documenter

DocMeta.setdocmeta!(BayesianDiscovery, :DocTestSetup, :(using BayesianDiscovery); recursive=true)

makedocs(;
    modules=[BayesianDiscovery],
    authors="Josh North <jsnowynorth@gmail.com>",
    repo="https://github.com/jsnowynorth/BayesianDiscovery.jl/blob/{commit}{path}#{line}",
    sitename="BayesianDiscovery.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "nothing") == "true",
        canonical="https://jsnowynorth.github.io/BayesianDiscovery.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Examples" => [
            "Burgers Equation" => "burgers.md",
            "Heat Equation" => "heat.md",
            "Reaction Diffusion" => "reactiondiffusion.md",
            "Feature Library" => "featurelibrary.md"
        ],
        "API" => "api.md"
    ],
)


deploydocs(;
    repo="github.com/jsnowynorth/BayesianDiscovery.jl.git",
    devbranch="main",
)

# using DocumenterTools
# DocumenterTools.genkeys(user="jsnowynorth", repo="git@github.com:jsnowynorth/BayesianDiscovery.jl.git")