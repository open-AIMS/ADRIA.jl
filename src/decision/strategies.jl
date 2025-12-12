abstract type DecisionStrategy end

"""
    is_decision_year(strategy::DecisionStrategy, timestep::Int64)::Bool

Fallback method that always returns true. Should be overridden by specific strategies.
"""
function is_decision_year(
    strategy::DecisionStrategy,
    timestep::Int64
)::Bool
    return true
end

"""
    PeriodicStrategy <: DecisionStrategy

Deploys at regular intervals starting from a specified year.

# Fields
- `target_locations::Vector{String}`: Optional list of location IDs to target
- `start_year::Int64`: Year to begin deployments (relative to simulation start)
- `duration::Int64`: Number of years to deploy for
- `frequency::Int64`: Frequency (in years) between deployments (0 = deploy once)
- `decision_years::Vector{Bool}`: Precomputed vector indicating deployment years
"""
struct PeriodicStrategy <: DecisionStrategy
    target_locations::Vector{String}
    start_year::Int64
    duration::Int64
    frequency::Union{Float64,Int64}  # in years
    decision_years::Vector{Bool}  # precomputed decision years

    function PeriodicStrategy(
        target_locations::Vector{String},
        start_year::Int64,
        duration::Int64,
        frequency::Union{Float64,Int64},
        timeframe::Int64
    )
        decision_years = decision_frequency(start_year, timeframe, duration, frequency)
        return new(target_locations, start_year, duration, frequency, decision_years)
    end
end

"""
    is_decision_year(strategy::PeriodicStrategy, timestep::Int64)::Bool

Check if deployment should occur at this timestep for periodic strategy.

# Arguments
- `strategy`: PeriodicStrategy with precomputed decision years
- `timestep`: Current simulation timestep

# Returns
Boolean indicating whether this is a deployment year
"""
function is_decision_year(
    strategy::PeriodicStrategy,
    timestep::Int64
)::Bool
    return strategy.decision_years[timestep]
end

"""
    filter_candidate_locations(
        strategy::PeriodicStrategy,
        timestep::Int64,
        state::NamedTuple
    )::Vector{String}

Return all target locations for periodic deployment strategy.

Periodic strategies deploy to all specified target locations whenever it is a decision
year, without filtering based on environmental conditions.

# Arguments
- `strategy`: PeriodicStrategy containing target locations
- `timestep`: Current simulation timestep
- `_`: Unused `state` argument to maintain compatibility

# Returns
Vector of location IDs eligible for deployment (all target locations)
"""
function filter_candidate_locations(
    strategy::PeriodicStrategy,
    timestep::Int64,
    _::Union{Nothing,NamedTuple}
)::Vector{String}
    if !is_decision_year(strategy, timestep)
        return String[]
    end

    return strategy.target_locations
end

"""
    AdaptiveStrategy <: DecisionStrategy

Deploys based on coral cover conditions - triggers on either:
1. Absolute threshold: when cover falls below a minimum level
2. Loss threshold: when cover drops by more than a percentage

Evaluate conditions every timestep.

# Fields
- `target_locations::Vector{String}`: Locations to consider (empty = all locations)
- `start_year::Int64`: Year to begin deployments (relative to simulation start)
- `duration::Int64`: Number of years to deploy for
- `absolute_threshold::Float64`: Deploy when cover < this value (e.g., 0.10 = 10%)
- `loss_threshold::Float64`: Deploy when proportional loss > this (e.g., 0.30 = 30% loss)
- `min_cover_remaining::Float64`: Don't deploy if cover below this (location unviable)
- `response_delay::Int64`: Timesteps to wait after trigger before deployment
"""
struct AdaptiveStrategy <: DecisionStrategy
    target_locations::Vector{String}
    start_year::Int64
    duration::Int64
    absolute_threshold::Float64
    loss_threshold::Float64
    min_cover_remaining::Float64
    response_delay::Int64
end

"""
    is_decision_year(strategy::AdaptiveStrategy, timestep::Int64)::Bool

Check if deployment should occur at this timestep for periodic strategy.

# Arguments
- `strategy`: AdaptiveStrategy with precomputed decision years
- `timestep`: Current simulation timestep

# Returns
Boolean indicating whether this is a deployment year
"""
function is_decision_year(
    strategy::AdaptiveStrategy,
    timestep::Int64
)::Bool
    if strategy.start_year .<= timestep .<= (strategy.start_year + strategy.duration - 1)
        return true
    end

    return false
end

"""
    filter_candidate_locations(
        strategy::AdaptiveStrategy,
        timestep::Int64,
        state::NamedTuple
    )::Vector{String}

Filter target locations to those meeting adaptive deployment criteria.

Returns locations where either:
1. Current cover is below absolute_threshold AND above min_cover_remaining
2. Recent cover loss exceeds loss_threshold AND current cover above min_cover_remaining

# Arguments
- `strategy`: AdaptiveStrategy with threshold parameters and target locations
- `timestep`: Current simulation timestep
- `state`: NamedTuple containing data for target locations only:
    - `current_cover::Vector{Float64}`: Current coral cover at each target location
    - `recent_cover_losses::Matrix{Float64}`: Cover losses [timestep ⋅ target_location]

# Returns
Vector of location IDs meeting deployment criteria (subset of target locations)

# Notes
Loss threshold check only occurs after (response_delay + lookback_window) timesteps
have elapsed to ensure sufficient historical data is available.
"""
function filter_candidate_locations(
    strategy::AdaptiveStrategy,
    timestep::Int64,
    state::NamedTuple
)::Vector{String}
    if !is_decision_year(strategy, timestep)
        return String[]
    end

    # Absolute threshold condition
    absolute_mask = (
        (state.current_cover .< strategy.absolute_threshold) .&
        (state.current_cover .> strategy.min_cover_remaining)
    )

    # Loss threshold condition
    recent_losses = if isempty(state.recent_cover_losses)
        trues(length(state.current_cover))
    else
        (first(state.recent_cover_losses) .>= strategy.loss_threshold)
    end

    above_min_cover = (state.current_cover .> strategy.min_cover_remaining)
    loss_mask = recent_losses .& above_min_cover

    # Combine with OR logic
    candidate_mask = absolute_mask .| loss_mask
    return strategy.target_locations[candidate_mask]
end
