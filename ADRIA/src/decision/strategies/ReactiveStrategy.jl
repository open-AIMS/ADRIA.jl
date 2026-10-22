"""
    ReactiveStrategy <: DecisionStrategy

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
- `revisit_cadence::Int64`: Timesteps before location becomes eligible again (0 = no restriction)
"""
struct ReactiveStrategy <: DecisionStrategy
    target_locations::Vector{String}
    start_year::Int64
    duration::Int64
    absolute_threshold::Float64
    loss_threshold::Float64
    min_cover_remaining::Float64
    response_delay::Int64
    revisit_cadence::Int64
end

"""
    is_decision_year(strategy::ReactiveStrategy, timestep::Int64)::Bool

Check if deployment should occur at this timestep for periodic strategy.

# Arguments
- `strategy`: ReactiveStrategy with precomputed decision years
- `timestep`: Current simulation timestep

# Returns
Boolean indicating whether this is a deployment year
"""
function is_decision_year(
    strategy::ReactiveStrategy,
    timestep::Int64
)::Bool
    if strategy.start_year .<= timestep .<= (strategy.start_year + strategy.duration - 1)
        return true
    end

    return false
end

"""
    filter_candidate_locations(
        strategy::ReactiveStrategy,
        timestep::Int64,
        state::NamedTuple
    )::Vector{String}

Filter target locations to those meeting deployment criteria.

Returns locations where EITHER condition is met:
1. Current cover is below absolute_threshold AND above min_cover_remaining
2. Recent cover loss exceeds loss_threshold AND current cover above min_cover_remaining

If revisit_cadence > 0, excludes locations deployed within the cadence period.

# Arguments
- `strategy`: ReactiveStrategy with threshold parameters and target locations
- `timestep`: Current simulation timestep
- `state`: NamedTuple containing data for target locations only:
    - `current_cover::Vector{Float64}`: Current coral cover at each target location
    - `recent_cover_losses::Matrix{Float64}`: Cover losses [timestep × target_location]
    - `last_deployment::Vector{Int64}`: Timestep of most recent deployment (0 if never), only required if revisit_cadence > 0

# Returns
Vector of location IDs meeting deployment criteria (subset of target locations)

# Notes
Loss threshold check only occurs after response_delay timesteps have elapsed.
When revisit_cadence = 0, behaves as pure reactive strategy (no spatial spreading).
When revisit_cadence > 0, prevents repeated deployments to same locations. Unlike
`PeriodicStrategy`, there is no backfill here — a reactive strategy should only ever
act on locations meeting its trigger conditions; cadence purely restricts that set
further.
"""
function filter_candidate_locations(
    strategy::ReactiveStrategy,
    timestep::Int64,
    state::NamedTuple
)::Vector{String}
    # Check if we are within the deployment timeframe
    if !is_decision_year(strategy, timestep)
        return String[]
    end

    targets = unique(strategy.target_locations)

    # Absolute threshold condition
    absolute_mask = (
        (state.current_cover .< strategy.absolute_threshold) .&
        (state.current_cover .> strategy.min_cover_remaining)
    )

    # Loss threshold condition
    recent_losses = if isempty(state.recent_cover_losses)
        trues(length(state.current_cover))
    else
        (state.recent_cover_losses .>= strategy.loss_threshold)
    end

    above_min_cover = (state.current_cover .> strategy.min_cover_remaining)
    loss_mask = recent_losses .& above_min_cover

    # Combine reactive criteria with OR logic
    reactive_mask = absolute_mask .| loss_mask

    # Apply cadence filter if revisit_cadence > 0
    if strategy.revisit_cadence > 0
        not_on_cooldown = revisit_cadence_mask(
            state.last_deployment, timestep, strategy.revisit_cadence
        )
        candidate_mask = reactive_mask .& not_on_cooldown
    else
        # No cadence restriction - pure reactive
        candidate_mask = reactive_mask
    end

    return targets[candidate_mask]
end
