# Add ADRIA src directory to PATH so files are found
push!(LOAD_PATH, "../src/")

using Documenter, ADRIA

makedocs(;
    sitename="ADRIA.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", nothing) == "true",
        assets=["assets/favicon.ico"]
    ),
    pages=[
        "index.md",
        "Concepts" =>
            [
                "concepts/synopsis.md",
                "concepts/dMCDA.md",
                "concepts/disturbances.md"
            ],
        "Usage" => [
            "usage/getting_started.md",
            "usage/domain.md",
            "usage/results.md",
            "usage/scenario_generation.md",
            "usage/scenario_runs.md",
            "usage/analysis.md",
            "usage/cookbook.md"
            # "usage/scenario_discovery.md"
        ],
        "Architecture" => [
            "architecture/architecture.md",
            "architecture/domain_and_resultsets.md"
        ],
        "Development" => [
            "development/development_setup.md",
            "development/metrics.md",
            "development/release_guide.md",
            "development/building_docs.md"],
        "API.md"
    ]
)

deploydocs(;
    repo="github.com/open-AIMS/ADRIA.jl.git",
    devbranch="main",
    target="build",
    push_preview=false
)
