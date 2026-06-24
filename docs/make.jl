using RSDP
using Documenter

DocMeta.setdocmeta!(RSDP, :DocTestSetup, :(using RSDP); recursive=true)

makedocs(;
    modules=[RSDP],
    authors="RSDP contributors",
    sitename="RSDP.jl",
    remotes=nothing,
    format=Documenter.HTML(;
        canonical="https://jeremypparker.github.io/RSDP.jl",
        edit_link="master",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Quick start" => "quickstart.md",
        "Manual" => [
            "Statuses and diagnostics" => "manual/statuses.md",
            "Exactification" => "manual/exactification.md",
            "Certificates" => "manual/certificates.md",
            "Cone checks" => "manual/cone_checks.md",
            "MOI extraction" => "manual/moi_extraction.md",
            "Facial reduction" => "manual/facial_reduction.md",
            "Limitations" => "manual/limitations.md",
        ],
        "Design" => [
            "Architecture" => "design/architecture.md",
            "Design spikes" => "design/design_spikes.md",
            "Weighted SOS" => "design/weighted_sos.md",
            "Solver selection" => "design/solver_selection.md",
            "Literature" => "design/literature.md",
            "Scaling" => "design/scaling.md",
        ],
        "API" => "api.md",
    ],
)

deploydocs(;
    repo="github.com/jeremypparker/RSDP.jl",
    devbranch="master",
)
