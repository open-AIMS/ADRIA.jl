"""
    build_seed_strategy(
        params::YAXArray, domain::Domain, locations::Vector{String}
    )::DecisionStrategy

Construct appropriate seeding strategy from scenario parameters.

# Arguments
- `params`: Scenario parameters (either DataFrameRow or YAXArray)
- `domain`: Domain object containing location information and timeframe
- `locations` : Location IDs to consider for deployment

# Returns
DecisionStrategy object (PeriodicStrategy or ReactiveStrategy)
"""
function build_seed_strategy(
    params::YAXArray, domain::Domain, locations::Vector{String}
)::DecisionStrategy
    strategy_type = Int64(params[At("seed_strategy")])
    if is_periodic(strategy_type)
        return PeriodicStrategy(
            locations,
            Int64(params[At("seed_year_start")]),
            Int64(params[At("seed_years")]),
            Int64(params[At("seed_deployment_freq")]),
            length(timesteps(domain))
        )
    elseif is_reactive(strategy_type)
        return ReactiveStrategy(
            locations,
            Int64(params[At("seed_year_start")]),
            Int64(params[At("seed_years")]),
            Float64(params[At("reactive_absolute_threshold")]),
            Float64(params[At("reactive_loss_threshold")]),
            Float64(params[At("reactive_min_cover_remaining")]),
            Int64(params[At("reactive_response_delay")]),
            Int64(params[At("reactive_cooldown_period")])
        )
    end

    return throw(ArgumentError("Unknown seed strategy type: $strategy_type"))
end

"""
    build_fog_strategy(
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
function build_fog_strategy(
    params::YAXArray, domain::Domain, locations::Vector{String}
)::DecisionStrategy
    strategy_type = Int64(params[At("fog_strategy")])

    if is_periodic(strategy_type)
        # Periodic Strategy
        return PeriodicStrategy(
            locations,
            Int64(params[At("fog_year_start")]),
            Int64(params[At("fog_years")]),
            Int64(params[At("fog_deployment_freq")]),
            length(timesteps(domain))
        )
    elseif is_reactive(strategy_type)
        # Reactive Strategy
        return ReactiveStrategy(
            locations,
            Int64(params[At("fog_year_start")]),
            Int64(params[At("fog_years")]),
            Float64(params[At("reactive_absolute_threshold")]),
            Float64(params[At("reactive_loss_threshold")]),
            Float64(params[At("reactive_min_cover_remaining")]),
            Int64(params[At("reactive_response_delay")]),
            Int64(params[At("reactive_cooldown_period")])
        )
    end

    return throw(ArgumentError("Unknown fog strategy type: $strategy_type"))
end
