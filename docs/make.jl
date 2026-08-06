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
        prettyurls=get(ENV, "CI", "false") == "true",
        assets=String[],
        size_threshold_ignore=["reference.md"],   # 120 alias docstrings on one page
    ),
    pages=[
        "Home" => "index.md",
        "Manual" => [
            "Formats" => "formats.md",
            "Projection" => "projection.md",
            "Operations" => "operations.md",
            "Stochastic rounding" => "stochastic.md",
            "Julia interop" => "interop.md",
        ],
        "API reference" => "reference.md",
    ],
    doctest=true,          # doctests are part of the build, not optional
    checkdocs=:exports,    # every exported name must be documented
    warnonly=false,        # any docs error fails the build
)

deploydocs(;
    repo="github.com/JeffreySarnoff/FloatBytes.jl",
    devbranch="main",
)
