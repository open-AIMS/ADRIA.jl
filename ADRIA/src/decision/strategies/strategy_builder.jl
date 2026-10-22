"""
    build_CAq_strategy(
        params::YAXArray, domain::Domain, locations::Vector{String}
    )::DecisionStrategy

Construct appropriate coral aquaculture strategy from scenario parameters.

# Arguments
- `params`: Scenario parameters (either DataFrameRow or YAXArray)
- `domain`: Domain object containing location information and timeframe
- `locations` : Location IDs to consider for deployment

# Returns
DecisionStrategy object (PeriodicStrategy or ReactiveStrategy)
"""
function build_CAq_strategy(
    params::YAXArray, domain::Domain, locations::Vector{String}
)::DecisionStrategy
    strategy_idx = Int64(params[At("iv_CAq_strategy")])
    strategy_params = build_strategy_params("iv_CAq", params, domain, locations)
    return build_strategy(strategy_type(strategy_idx), strategy_params)
end

"""
    build_Fog_strategy(
        params::YAXArray, domain::Domain, locations::Vector{String}
    )::DecisionStrategy

Construct appropriate fogging strategy from scenario parameters.

# Arguments
- `params`: Scenario parameters (either DataFrameRow or YAXArray)
- `domain`: Domain object containing location information and timeframe
- `locations` : Location IDs to consider for deployment

# Returns
DecisionStrategy object (PeriodicStrategy or ReactiveStrategy)
"""
function build_Fog_strategy(
    params::YAXArray, domain::Domain, locations::Vector{String}
)::DecisionStrategy
    strategy_idx = Int64(params[At("iv_Fog_strategy")])
    strategy_params = build_strategy_params("iv_Fog", params, domain, locations)
    return build_strategy(strategy_type(strategy_idx), strategy_params)
end

"""
    build_LvM_strategy(
        params::YAXArray, domain::Domain, locations::Vector{String}
    )::DecisionStrategy

Construct appropriate larval methods strategy from scenario parameters.

# Arguments
- `params`: Scenario parameters (either DataFrameRow or YAXArray)
- `domain`: Domain object containing location information and timeframe
- `locations` : Location IDs to consider for deployment

# Returns
DecisionStrategy object (PeriodicStrategy or ReactiveStrategy)
"""
function build_LvM_strategy(
    params::YAXArray, domain::Domain, locations::Vector{String}
)::DecisionStrategy
    strategy_idx = Int64(params[At("iv_LvM_strategy")])
    strategy_params = build_strategy_params("iv_LvM", params, domain, locations)
    return build_strategy(strategy_type(strategy_idx), strategy_params)
end

function build_strategy_params(prefix::String, params::YAXArray, domain::Domain, locations)
    return (
        locations=locations,
        timesteps=length(timesteps(domain)),
        iv_year_start=Int64(params[At("$(prefix)_year_start")]),
        iv_years=Int64(params[At("$(prefix)_years")]),
        iv_deployment_freq=Int64(params[At("$(prefix)_deployment_freq")]),
        reactive_absolute_threshold=Float64(params[At("reactive_absolute_threshold")]),
        reactive_loss_threshold=Float64(params[At("reactive_loss_threshold")]),
        reactive_min_cover_remaining=Float64(params[At("reactive_min_cover_remaining")]),
        reactive_response_delay=Int64(params[At("reactive_response_delay")]),
        revisit_cadence=Int64(params[At("$(prefix)_revisit_cadence")]),
        min_locations=Int64(
            params[At(prefix == "iv_LvM" ? "iv_LvM_min_iv_locations" : "min_iv_locations")]
        )
    )
end

function build_strategy(
    ::Type{PeriodicStrategy}, params::NamedTuple
)::PeriodicStrategy
    return PeriodicStrategy(
        params.locations,
        params.iv_year_start,
        params.iv_years,
        params.iv_deployment_freq,
        params.timesteps,
        params.revisit_cadence,
        params.min_locations
    )
end
function build_strategy(
    ::Type{ReactiveStrategy}, params::NamedTuple
)::ReactiveStrategy
    return ReactiveStrategy(
        params.locations,
        params.iv_year_start,
        params.iv_years,
        params.reactive_absolute_threshold,
        params.reactive_loss_threshold,
        params.reactive_min_cover_remaining,
        params.reactive_response_delay,
        params.revisit_cadence
    )
end
