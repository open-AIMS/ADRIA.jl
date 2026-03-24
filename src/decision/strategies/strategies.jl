abstract type DecisionStrategy end

include("PeriodicStrategy.jl")
include("ReactiveStrategy.jl")
# include("PreventativeStrategy.jl")
include("strategy_builder.jl")

const DECISION_STRATEGY = Dict(
    :periodic => 1,
    :reactive => 2,
    1 => :periodic,
    2 => :reactive
)

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

function strategy_type(strategy_idx::Int64)
    if is_reactive(strategy_idx)
        return ReactiveStrategy
    elseif is_periodic(strategy_idx)
        return PeriodicStrategy
    end
    throw(ArgumentError("Unknown mc strategy type: $strategy_type"))
end

function is_reactive(strategy_idx::T)::Bool where {T<:Real}
    return strategy_idx == DECISION_STRATEGY[:reactive]
end
function is_reactive(strategy_idx::AbstractVector{T})::BitVector where {T<:Real}
    return is_reactive.(strategy_idx)
end

function is_periodic(strategy_idx::T)::Bool where {T<:Real}
    return strategy_idx == DECISION_STRATEGY[:periodic]
end
function is_periodic(strategy_idx::AbstractVector{T})::BitVector where {T<:Real}
    return is_periodic.(strategy_idx)
end

function build_state(domain, strategy, states)
    target_loc_indices = findall(
        in.(domain.loc_ids, Ref(strategy.target_locations))
    )
    return strategy_status(strategy, states, target_loc_indices)
end

function strategy_status(::DecisionStrategy, states, idx)
    return (
        current_cover=states.current_cover[idx],
        recent_cover_losses=first(states.recent_cover_losses)[idx]
    )
end

function strategy_status(::ReactiveStrategy, states, idx)
    return (
        current_cover=states.current_cover[idx],
        recent_cover_losses=first(states.recent_cover_losses)[idx],
        last_deployment=states.last_deployment[idx]
    )
end
