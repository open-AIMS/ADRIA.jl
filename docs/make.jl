push!(LOAD_PATH, "../src/")

using Documenter, ADRIA


makedocs(sitename="ADRIA Documentation",
    format=Documenter.HTML(
        prettyurls=get(ENV, "CI", nothing) == "true",
        sidebar_sitename=false,
        assets=["assets/favicon.ico"]
    ),
    pages=[
        "index.md",
        "synopsis.md",
        "Usage" => [
            "usage/getting_started.md",
            "usage/domain.md",
            "usage/scenario_generation.md",
            "usage/scenario_runs.md",
            "usage/analysis.md",
        ],
        "dMCDA.md",
        # "Examples" => [
        #     "Usage" => [
        #     ],
        #     "simple_showcase.md",
        #     "advanced_showcase.md"
        # ],
        "Development" => [
            "development/development_setup.md",
            "development/architecture.md",
            "development/metrics.md",
            "development/release_guide.md",
            "development/building_docs.md"],
        "API.md"
    ]
)

deploydocs(
    repo="github.com/open-AIMS/ADRIA.jl.git",
    devbranch="main",
    target="build",
    push_preview=false
)
