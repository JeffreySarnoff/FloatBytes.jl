using FloatBytes
using Documenter

DocMeta.setdocmeta!(FloatBytes, :DocTestSetup, :(using FloatBytes); recursive=true)

makedocs(;
    modules=[FloatBytes],
    authors="Jeffrey Sarnoff <jeffrey.sarnoff@gmail.com> and contributors",
    sitename="FloatBytes.jl",
    format=Documenter.HTML(;
        canonical="https://JeffreySarnoff.github.io/FloatBytes.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/JeffreySarnoff/FloatBytes.jl",
    devbranch="main",
)
