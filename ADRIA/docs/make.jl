push!(LOAD_PATH, "../src/")

using Documenter, DocumenterVitepress, Literate, ADRIA, ADRIAanalysis

# Generate documentation pages from Literate sources
usage_dir = joinpath(@__DIR__, "src", "usage")
for name in [
    "getting_started",
    "loading_results",
    "scenario_runs",
    "scenario_discovery",
    "analysis",
    "cookbook"
]
    Literate.markdown(
        joinpath(usage_dir, "$(name).jl"),
        usage_dir;
        flavor=Literate.DocumenterFlavor()
    )
end

makedocs(;
    sitename="ADRIA.jl",
    format=DocumenterVitepress.MarkdownVitepress(;
        repo="github.com/open-AIMS/ADRIA.jl",
        devbranch="main"
    ),
    pages=[
        "index.md",
        "Concepts" => [
            joinpath("concepts", "dMCDA.md"),
            joinpath("concepts", "disturbances.md")
        ],
        "Usage" => [
            joinpath("usage", "getting_started.md"),
            joinpath("usage", "loading_a_domain.md"),
            joinpath("usage", "loading_results.md"),
            joinpath("usage", "generating_scenarios.md"),
            joinpath("usage", "scenario_runs.md"),
            joinpath("usage", "scenario_discovery.md"),
            joinpath("usage", "analysis.md"),
            joinpath("usage", "cookbook.md")
        ],
        "Architecture" => [
            joinpath("architecture", "architecture.md"),
            joinpath("architecture", "domain_and_resultsets.md")
        ],
        "Development" => [
            joinpath("development", "development_setup.md"),
            joinpath("development", "metrics.md"),
            joinpath("development", "release_guide.md"),
            joinpath("development", "building_docs.md")
        ],
        "API.md"
    ]
)

# On Windows, Documenter normalizes page paths with normpath(), producing backslash
# separators. DocumenterVitepress writes these directly into config.mts nav links and
# into cross-reference hrefs in the compiled markdown. JavaScript treats \g, \a, etc.
# as escape sequences, collapsing /usage\getting_started into /usagegetting_started.
# Fix by normalizing separators to forward slashes in all generated output files.
if Sys.iswindows()
    build_dir = joinpath(@__DIR__, "build", ".documenter")

    # Fix markdown cross-reference links: only replace \ inside ](/...) destinations
    # to avoid corrupting LaTeX math (e.g. \alpha) elsewhere in the content.
    for (root, _, files) in walkdir(build_dir)
        for file in files
            endswith(file, ".md") || continue
            path = joinpath(root, file)
            content = read(path, String)
            fixed = replace(content, r"\]\(/[^)]*\)" => m -> replace(m, '\\' => '/'))
            fixed != content && write(path, fixed)
        end
    end

    # Fix nav/sidebar link values in config.mts (safe global replace — no math here).
    config_path = joinpath(build_dir, ".vitepress", "config.mts")
    if isfile(config_path)
        content = read(config_path, String)
        write(config_path, replace(content, '\\' => '/'))
    end
end

deploydocs(;
    repo="github.com/open-AIMS/ADRIA.jl.git",
    devbranch="main",
    target="build",
    push_preview=false
)
