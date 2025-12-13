"""
    build_seed_strategy(params::DataFrameRow, domain::Domain, locations::Vector{String})::DecisionStrategy
    build_seed_strategy(params::YAXArray, domain::Domain, locations::Vector{String})::DecisionStrategy

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
    strategy_type = try
        params[At("seed_strategy")]
    catch
        throw(ArgumentError("No strategy type defined"))
    end

    # NamedTuple(zip(Symbol.(params.factors), params))
    return return_seed_strategy(Int64(strategy_type), params, domain, locations)
end
function build_seed_strategy(
    params::DataFrameRow, domain::Domain, locations::Vector{String}
)::DecisionStrategy
    # Not sure I need this implementation for DataFrameRow...
    strategy_type = try
        params.seed_strategy
    catch
        throw(ArgumentError("No strategy type defined"))
    end

    return return_seed_strategy(
        Int64(strategy_type), params, domain, locations
    )
end

function return_seed_strategy(
    strategy_type::Int64, params::YAXArray, domain::Domain, locations::Vector{String}
)::DecisionStrategy
    if strategy_type == 0
        return PeriodicStrategy(
            locations,
            Int64(params[At("seed_year_start")]),
            Int64(params[At("seed_years")]),
            Int64(params[At("seed_deployment_freq")]),
            length(timesteps(domain))
        )
    elseif strategy_type == 1
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

    throw(ArgumentError("Unknown seed strategy type: $strategy_type"))
end

"""
    build_fog_strategy(params::DataFrameRow, domain::Domain)::DecisionStrategy
    build_fog_strategy(params::YAXArray, domain::Domain)::DecisionStrategy

Construct appropriate fogging strategy from scenario parameters.

# Arguments
- `params`: Scenario parameters (either DataFrameRow or YAXArray)
- `domain`: Domain object containing location information and timeframe

# Returns
DecisionStrategy object (PeriodicStrategy or ReactiveStrategy)
"""
function build_fog_strategy(
    params::Union{DataFrameRow,YAXArray},
    domain::Domain,
    locations::Vector{String}
)::DecisionStrategy
    strategy_type = Int64(params[At("fog_strategy")])

    if strategy_type == 0
        # Periodic Strategy
        return PeriodicStrategy(
            locations,
            Int64(params[At("fog_year_start")]),
            length(timesteps(domain)),
            Int64(params[At("fog_years")]),
            Int64(params[At("fog_deployment_freq")])
        )
    elseif strategy_type == 1
        # Reactive Strategy
        return ReactiveStrategy(
            locations,
            Float64(params[At("reactive_absolute_threshold")]),
            Float64(params[At("reactive_loss_threshold")]),
            Float64(params[At("reactive_min_cover_remaining")]),
            Int64(params[At("reactive_response_delay")]),
            Int64(params[At("reactive_cooldown_period")])
        )
    end

    throw(ArgumentError("Unknown fog strategy type: $strategy_type"))
end
