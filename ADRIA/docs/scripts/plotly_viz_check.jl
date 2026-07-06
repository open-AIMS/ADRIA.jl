# ADRIA/docs/scripts/plotly_viz_check.jl
#
# Plotly backend wrapper for viz_check_common.jl.
#
# Exercises every ADRIA visualization using the PlotlyBase backend.
# Each figure is written to ADRIA/docs/scripts/plotly_viz_output/<name>.html.
#
# Usage (from repo root):
#   ADRIA_TEST_DOMAIN=/path/to/domain julia --project=. ADRIA/docs/scripts/plotly_viz_check.jl
#
# Required environment variables:
#   ADRIA_TEST_DOMAIN  Path to the domain directory to load.
#
# Optional environment variables:
#   ADRIA_OUTPUT_DIR   Where ADRIA writes/reads result sets (default: ./Outputs).
#   OPEN_IN_BROWSER    Set to "true" to open each HTML figure in the browser.
#
# Heavy dependencies not in ADRIAviz test Project.toml:
#   MLJ, SIRUS, ADRIAanalysis
# These must be available in the active Julia project when running this script.

# ----------------------------------------------------------------------------
# Backend setup
# ----------------------------------------------------------------------------

const RESULTS = @NamedTuple{name::String, ok::Bool, msg::String, elapsed::Float64}[]

const OUT_DIR = joinpath(@__DIR__, "plotly_viz_output")
isdir(OUT_DIR) || mkpath(OUT_DIR)

const OPEN_IN_BROWSER = get(ENV, "OPEN_IN_BROWSER", "false") == "true"

const BACKEND_NAME = "plotly"

_activate_backend() = ADRIAviz.activate(:plotly)

function _save_fig(fig, name::String)
    fname = replace(name, r"[^A-Za-z0-9]+" => "_") * ".html"
    path = joinpath(OUT_DIR, fname)
    open(path, "w") do io
        show(io, MIME"text/html"(), fig)
    end
    OPEN_IN_BROWSER && ADRIA.viz.show_in_browser(fig)
    return path
end

# ----------------------------------------------------------------------------
# Backend-specific viz dispatch
# ----------------------------------------------------------------------------

_viz_tsa(rs, si) = ADRIA.viz.tsa(si)

_viz_rsa(scens, y, foi) = ADRIA.viz.rsa(scens, y, foi)
_viz_outcome_map(scens, y, foi) = ADRIA.viz.outcome_map(scens, y, foi)

_viz_rules_scatter(rs, X, tgt, rules) = ADRIA.viz.rules_scatter(X, tgt, rules)

_viz_tsc_map(rs, tac_sites, tac_clusters) = ADRIA.viz.map(rs, tac_sites, tac_clusters)

_viz_backend_extras() = nothing

# ----------------------------------------------------------------------------
# Run shared checks
# ----------------------------------------------------------------------------

include("viz_check_common.jl")

println("HTML output written to: $OUT_DIR")
