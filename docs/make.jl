using Aeroacoustics
using Documenter

DocMeta.setdocmeta!(Aeroacoustics, :DocTestSetup, :(using Aeroacoustics); recursive=true)

makedocs(;
    modules=[Aeroacoustics],
    authors="Matthew Priddin and contributors",
    repo="https://github.com/mjp98/Aeroacoustics.jl/blob/{commit}{path}#{line}",
    sitename="Aeroacoustics.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://mjp98.github.io/Aeroacoustics.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/mjp98/Aeroacoustics.jl",
    devbranch="main",
)
