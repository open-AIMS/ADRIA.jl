push!(LOAD_PATH, "../src/")

using Documenter, ADRIA


makedocs(sitename="ADRIA Documentation",
    format=Documenter.HTML(
        prettyurls=get(ENV, "CI", nothing) == "true"
    ),
    pages=[
        "index.md",
        "synopsis.md",
        "dMCDA.md",
        # "Examples" => [
        #     "Usage" => [
        #     ],
        #     "simple_showcase.md",
        #     "advanced_showcase.md"
        # ],
        "Development" => [
            "development_setup.md",
            "architecture.md",
        ],
        "API.md"
    ]
)