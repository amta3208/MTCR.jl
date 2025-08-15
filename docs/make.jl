using MTCR
using Documenter

DocMeta.setdocmeta!(MTCR, :DocTestSetup, :(using MTCR); recursive=true)

makedocs(;
    modules=[MTCR],
    authors="Amin Taziny",
    sitename="MTCR.jl",
    format=Documenter.HTML(;
        canonical="https://amta3208.github.io/MTCR.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/amta3208/MTCR.jl",
    devbranch="main",
)
