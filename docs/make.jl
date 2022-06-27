using QuasiPeriodicHelmholtz
using Documenter

DocMeta.setdocmeta!(QuasiPeriodicHelmholtz, :DocTestSetup, :(using QuasiPeriodicHelmholtz); recursive=true)

makedocs(;
    modules=[QuasiPeriodicHelmholtz],
    authors="Matthew Priddin and contributors",
    repo="https://github.com/mjp98/QuasiPeriodicHelmholtz.jl/blob/{commit}{path}#{line}",
    sitename="QuasiPeriodicHelmholtz.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://mjp98.github.io/QuasiPeriodicHelmholtz.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/mjp98/QuasiPeriodicHelmholtz.jl",
    devbranch="main",
)
