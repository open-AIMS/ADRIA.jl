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

function build_state(
    domain, strategy, current_loc_cover, recent_cover_losses, last_deployment
)
    # Get target location indices for this strategy
    target_loc_indices = findall(
        in.(domain.loc_ids, Ref(strategy.target_locations))
    )

    return if isa(strategy, ReactiveStrategy)
        # Need to track deployment history
        (
            current_cover=current_loc_cover[target_loc_indices],
            recent_cover_losses=first(recent_cover_losses)[target_loc_indices],
            last_deployment=last_deployment[target_loc_indices]
        )
    else
        (
            current_cover=current_loc_cover[target_loc_indices],
            recent_cover_losses=first(recent_cover_losses)[target_loc_indices]
        )
    end
end
