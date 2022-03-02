using HyperSparseMatrices
using Documenter

DocMeta.setdocmeta!(HyperSparseMatrices, :DocTestSetup, :(using HyperSparseMatrices); recursive=true)

makedocs(;
    modules=[HyperSparseMatrices],
    authors="Wimmerer <kimmerer@mit.edu> and contributors",
    repo="https://github.com/Wimmerer/HyperSparseMatrices.jl/blob/{commit}{path}#{line}",
    sitename="HyperSparseMatrices.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://Wimmerer.github.io/HyperSparseMatrices.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/Wimmerer/HyperSparseMatrices.jl",
    devbranch="main",
)
