abstract type DecisionStrategy end

include("PeriodicStrategy.jl")
include("ReactiveStrategy.jl")
# include("PreventativeStrategy.jl")
include("strategy_builder.jl")

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
