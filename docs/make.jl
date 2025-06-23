using SunlightHNC
using Documenter

DocMeta.setdocmeta!(SunlightHNC, :DocTestSetup, :(using SunlightHNC); recursive=true)

makedocs(;
    modules=[SunlightHNC],
    authors="Hiroharu Sugawara <hsugawa@gmail.com> and contributors",
    sitename="SunlightHNC.jl",
    format=Documenter.HTML(;
        canonical="https://hsugawa8651.github.io/SunlightHNC.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/hsugawa8651/SunlightHNC.jl",
    devbranch="main",
)
