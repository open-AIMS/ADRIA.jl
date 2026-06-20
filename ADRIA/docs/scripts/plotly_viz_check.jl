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
#
# Plotly tsa takes Si directly (no ResultSet); rsa/outcome_map take stacked
# YAXArrays rather than the Dataset returned by the analysis functions.
# ----------------------------------------------------------------------------

import ADRIA.YAXArrays: At

function _stack_rsa(rsa_ds, factors::Vector{Symbol}, rs)
    si_q = collect(rsa_ds[factors[1]].axes[1])
    rows = [coalesce.(collect(rsa_ds[f][Si = At("Si")]), NaN) for f in factors]
    mat = permutedims(reduce(hcat, rows))
    Si = ADRIA.DataCube(Matrix{Float64}(mat); factors=factors, si_quantile=si_q)
    fvals = Matrix{Float64}(rs.inputs[:, factors])
    return Si, fvals
end

function _stack_outcome_map(om_ds, factors::Vector{Symbol}, rs)
    si_q = collect(om_ds[factors[1]].axes[1])
    CI = [:lower, :mean, :upper]
    arr = Array{Float64}(undef, length(factors), 3, length(si_q))
    for (i, f) in enumerate(factors)
        da = om_ds[f]
        arr[i, 1, :] = coalesce.(collect(da[CI = At("lower")]), NaN)
        arr[i, 2, :] = coalesce.(collect(da[CI = At("mean")]), NaN)
        arr[i, 3, :] = coalesce.(collect(da[CI = At("upper")]), NaN)
    end
    outcomes = ADRIA.DataCube(arr; factors=factors, CI=CI, si_quantile=si_q)
    fvals = Matrix{Float64}(rs.inputs[:, factors])
    return outcomes, fvals
end

_viz_tsa(rs, si) = ADRIA.viz.tsa(si)

function _viz_rsa(rs, rsa_ds, foi)
    Si, fvals = _stack_rsa(rsa_ds, foi, rs)
    return ADRIA.viz.rsa(Si, fvals)
end

function _viz_outcome_map(rs, om_ds, foi)
    outcomes, fvals = _stack_outcome_map(om_ds, foi, rs)
    return ADRIA.viz.outcome_map(outcomes, fvals)
end

_viz_rules_scatter(rs, scens, tgt, rules) = ADRIA.viz.rules_scatter(scens, tgt, rules)

_viz_tsc_map(rs, tac_sites, tac_clusters) = ADRIA.viz.map(rs, tac_sites, tac_clusters)

_viz_backend_extras() = nothing

# ----------------------------------------------------------------------------
# Run shared checks
# ----------------------------------------------------------------------------

include("viz_check_common.jl")

println("HTML output written to: $OUT_DIR")
