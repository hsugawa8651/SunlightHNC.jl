using SunlightHNC
using Documenter

DocMeta.setdocmeta!(SunlightHNC, :DocTestSetup, :(using SunlightHNC); recursive=true)

makedocs(;
    modules=[SunlightHNC],
    authors="Hiroharu Sugawara, Tokyo Metropolitan University <hsugawa@tmu.ac.jp>",
    sitename="SunlightHNC.jl",
    format=Documenter.HTML(;
        canonical="https://hsugawa8651.github.io/SunlightHNC.jl",
        edit_link="master",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Manual" => [
            "Wizard" => "wizard.md",
            "Models" => "models.md",
            "Plots" => "plots.md",
        ],
        "References" => "references.md",
    ],
)

deploydocs(;
    repo="github.com/hsugawa8651/SunlightHNC.jl",
    devbranch="master",
)
