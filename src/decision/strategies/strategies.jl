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
