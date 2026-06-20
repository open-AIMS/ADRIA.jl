# ADRIA/docs/scripts/makie_viz_check.jl
#
# Makie (CairoMakie) backend wrapper for viz_check_common.jl.
#
# Exercises every ADRIA visualization using the CairoMakie backend, saving
# each figure as a PNG. Serves two purposes:
#
#   1. Visual smoke check — verify all Makie visualizations render without
#      errors against a real domain and result set.
#
#   2. Documentation figure generation — the output PNGs are committed to
#      docs/src/assets/imgs/analysis/ and used in the documentation pages.
#      Run this script offline whenever documentation figures need refreshing.
#
# This script should NOT be run on CI (compute-heavy; requires heavy deps).
#
# Usage (from repo root):
#   ADRIA_TEST_DOMAIN=/path/to/domain julia --project=. ADRIA/docs/scripts/makie_viz_check.jl
#
# Required environment variables:
#   ADRIA_TEST_DOMAIN  Path to the domain directory to load.
#
# Optional environment variables:
#   ADRIA_OUTPUT_DIR   Where ADRIA writes/reads result sets (default: ./Outputs).
#   ADRIA_FIGURE_DIR   Where to save output PNGs.
#                      Default: docs/src/assets/imgs/analysis/ (relative to
#                               the docs/ directory).
#   ADRIA_FIG_FORMAT   "png" (default) or "svg".
#
# Heavy dependencies not in ADRIAviz test Project.toml:
#   CairoMakie, MLJ, SIRUS, ADRIAanalysis
# These must be available in the active Julia project when running this script.
#
# Model outputs (Zarr result sets in Outputs/) are gitignored.
# The generated PNG figures should be committed when updated.

# ----------------------------------------------------------------------------
# Backend setup
# ----------------------------------------------------------------------------

using CairoMakie

const RESULTS = @NamedTuple{name::String, ok::Bool, msg::String, elapsed::Float64}[]

const FIG_FORMAT = get(ENV, "ADRIA_FIG_FORMAT", "png")

const DEFAULT_FIG_DIR = joinpath(@__DIR__, "makie_viz_output")
const FIG_DIR = get(ENV, "ADRIA_FIGURE_DIR", DEFAULT_FIG_DIR)
isdir(FIG_DIR) || mkpath(FIG_DIR)

const BACKEND_NAME = "makie"

_activate_backend() = ADRIAviz.activate("CairoMakie")

function _save_fig(fig, name::String)
    fname = replace(name, r"[^A-Za-z0-9]+" => "_") * "." * FIG_FORMAT
    path = joinpath(FIG_DIR, fname)
    save(path, fig)
    return path
end

# ----------------------------------------------------------------------------
# Backend-specific viz dispatch
#
# Makie tsa/rsa/outcome_map take a ResultSet as first arg; rules_scatter also
# takes rs before scens.
# ----------------------------------------------------------------------------

_viz_tsa(rs, si) = ADRIA.viz.tsa(rs, si)
_viz_rsa(rs, rsa_ds, foi) = ADRIA.viz.rsa(rs, rsa_ds, foi)
_viz_outcome_map(rs, om_ds, foi) = ADRIA.viz.outcome_map(rs, om_ds, foi)
_viz_rules_scatter(rs, scens, tgt, rules) = ADRIA.viz.rules_scatter(rs, scens, tgt, rules)
_viz_tsc_map(rs, tac_sites, tac_clusters) = ADRIA.viz.map(rs, tac_sites, tac_clusters)

function _viz_backend_extras()
    check("convergence_components_heatmap") do
        outcome = dropdims(mean(s_tac; dims=:timesteps); dims=:timesteps)
        conv_foi = Symbol[
            f for f in (:dhw_scenario, :guided) if string(f) in names(rs.inputs)
        ]
        isempty(conv_foi) && (conv_foi = foi[1:min(2, length(foi))])
        Si_conv = convergence(scens, outcome, conv_foi)
        ADRIA.viz.convergence(
            Si_conv, conv_foi; opts=Dict{Symbol,Any}(:viz_type => :heatmap)
        )
    end

    check("ranks_plot") do
        rank_freq = ADRIA.decision.ranks_to_frequencies(ADRIA.metrics.seed_ranks(rs))
        ADRIA.viz.ranks_to_frequencies(rs, rank_freq, 1)
    end
end

# ----------------------------------------------------------------------------
# Run shared checks
# ----------------------------------------------------------------------------

include("viz_check_common.jl")

println("Figures written to: $FIG_DIR")
if count(r -> !r[2], RESULTS) < length(RESULTS)
    println("\nTo commit updated docs figures:")
    println("  Copy the figures to ADRIA/docs/src/assets/imgs/analysis/")
    println("  git add ADRIA/docs/src/assets/imgs/analysis/*.$FIG_FORMAT")
    println("  git commit -m 'docs: refresh Makie example figures'")
end
