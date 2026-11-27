module decision

using Random
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
    location_out_spread,
    timesteps

using ADRIA:
    DiscreteOrderedUniformDist,
    Factor,
    _factor_spec

using
    Combinatorics,
    DataFrames,
    JMcDM

include("dMCDA.jl")
include("Criteria/DecisionWeights.jl")
include("Criteria/DecisionPreferences.jl")
include("location_selection.jl")
include("strategies/strategies.jl")

export
    # Intervention preferences
    DecisionPreferences,
    # ShdPreferences,
    # Criteria definitions
    Criterion,
    criterion_names,
    criterion_fields,
    # Helper Methods
    decision_matrix,
    filter_criteria,
    update_criteria_values!,
    select_locations,
    unguided_selection,
    decision_frequency,
    weighted_projection,
    build_decay,
    identify_within_depth_bounds,
    summary_stat_env,
    mcda_methods,
    # Strategies
    DECISION_STRATEGY,
    is_reactive,
    is_periodic,
    is_decision_year,
    build_CAq_strategy,
    build_Fog_strategy,
    build_LvM_strategy,
    build_state,
    filter_candidate_locations,
    PeriodicStrategy,
    ReactiveStrategy
end
