using HypersparseMatrices
using Documenter

DocMeta.setdocmeta!(HypersparseMatrices, :DocTestSetup, :(using HypersparseMatrices); recursive=true)

makedocs(;
    modules=[HypersparseMatrices],
    authors="Wimmerer <kimmerer@mit.edu> and contributors",
    repo="https://github.com/Wimmerer/HypersparseMatrices.jl/blob/{commit}{path}#{line}",
    sitename="HypersparseMatrices.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://Wimmerer.github.io/HypersparseMatrices.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/Wimmerer/HypersparseMatrices.jl",
    devbranch="main",
)
