module decision

using InteractiveUtils: subtypes
using StatsBase
using YAXArrays
using ADRIA:
    component_params,
    DataCube,
    Domain,
    EcoModel,
    n_locations,
    loc_area,
    loc_k_area,
    timesteps

using ADRIA:
    DiscreteOrderedUniformDist,
    Factor

using
    Combinatorics,
    DataFrames,
    JMcDM

include("dMCDA.jl")
include("Criteria/DecisionPreferences.jl")
include("Criteria/DecisionWeights.jl")
include("location_selection.jl")
include("strategies/strategies.jl")

export
    # Intervention preferences
    SeedPreferences,
    FogPreferences,
    # SRMPreferences,
    # Helper Methods
    decision_matrix,
    filter_criteria,
    update_criteria_values!,
    select_locations,
    unguided_selection,
    decision_frequency,
    weighted_projection,
    identify_within_depth_bounds,
    summary_stat_env,
    mcda_methods,
    # Strategies
    DECISION_STRATEGY,
    is_reactive,
    is_periodic,
    is_decision_year,
    build_seed_strategy,
    build_fog_strategy,
    filter_candidate_locations,
    PeriodicStrategy,
    ReactiveStrategy
end
