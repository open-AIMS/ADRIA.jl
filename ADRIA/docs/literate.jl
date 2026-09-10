using Literate

# Regenerate the Literate-derived usage pages. Run standalone for a fast preview
# (`julia --project=docs docs/literate.jl`); also included by `make.jl`.

usage_dir = joinpath(@__DIR__, "src", "usage")
usage_pages = [
    "getting_started",
    "loading_a_domain",
    "loading_results",
    "generating_scenarios",
    "scenario_runs",
    "scenario_discovery",
    "analysis",
    "cookbook",
    "exporting_to_rme"
]

for name in usage_pages
    Literate.markdown(
        joinpath(usage_dir, "$(name).jl"),
        usage_dir;
        flavor=Literate.DocumenterFlavor(),
        codefence=("```julia" => "```")
    )
end
