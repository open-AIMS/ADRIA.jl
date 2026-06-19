module ADRIAviz

using ADRIA
using ADRIA: ResultSet, RMEResultSet, AnnotatedOutcomes, attach_scenario_metadata
using OrderedCollections, Statistics, DataFrames, Distributions
using PrecompileTools

include("activate.jl")
include("analysis.jl")
include("_scenario_helpers.jl")
include("theme.jl")
include("viz/viz.jl")

# Re-export shared spatial utilities for extensions and downstream consumers
export _loc_id_col, _get_site_ids, _site_ids, _haversine_km, _nice_length
export GBR_COASTAL_PLACES, CoastalPlace, MapDecorationData, compute_map_decorations

@compile_workload begin
    # Scenario grouping helpers — pure DataFrame operations, no backend
    _scenario_rcps(DataFrame(:RCP => [45, 45, 85, 85]))
    _scenario_clusters(BitVector([true, false, true, false]))
    _scenario_clusters([1, 1, 2, 2])
end
end
