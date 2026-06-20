module ADRIAviz

using ADRIA
using ADRIA: ResultSet, RMEResultSet, AnnotatedOutcomes, attach_scenario_metadata
using OrderedCollections, Statistics, DataFrames, Distributions
using YAXArrays
using PrecompileTools

include("activate.jl")
include("analysis.jl")
include("_scenario_helpers.jl")
include("theme.jl")
include("viz/viz.jl")

# Import symbols from viz submodule to make them available at ADRIAviz level
using .viz:
    outcome_title, outcome_label, set_plot_opts!,
    OPT_TYPE, DEFAULT_OPT_TYPE,
    _time_labels, _calc_gridsize, set_typography_defaults!, timesteps,
    _loc_id_col, _get_site_ids, _site_ids, _haversine_km, _nice_length, validate_extent,
    GBR_COASTAL_PLACES, CoastalPlace, MapDecorationData, compute_map_decorations

# Re-export shared spatial utilities and viz symbols for extensions and downstream consumers
export _loc_id_col, _get_site_ids, _site_ids, _haversine_km, _nice_length, validate_extent
export GBR_COASTAL_PLACES, CoastalPlace, MapDecorationData, compute_map_decorations
export outcome_title, outcome_label, set_plot_opts!
export OPT_TYPE, DEFAULT_OPT_TYPE
export _time_labels, _calc_gridsize, set_typography_defaults!, timesteps

# Re-export theme and analysis helpers
export COLORS, labels
export _get_scenario_groups, _scenario_types, _scenario_rcps, _scenario_clusters
export relative_sensitivities, outcome_probability

@compile_workload begin
    # Scenario grouping helpers — pure DataFrame operations, no backend
    _scenario_rcps(DataFrame(:RCP => [45, 45, 85, 85]))
    _scenario_clusters(BitVector([true, false, true, false]))
    _scenario_clusters([1, 1, 2, 2])
end
end
