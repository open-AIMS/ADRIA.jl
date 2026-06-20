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
    make_decision = try
        strategy.decision_years[timestep]
    catch err
        if !(err isa BoundsError)
            rethrow(err)
        end

        false
    end

    return make_decision
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
