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

"""
    strategy_type(strategy_idx::Int64)
    strategy_type(param_set::YAXArray{Float64,1}, iv_type::String)

#  Arguments
- `strategy_idx` : Can be either `1` (periodic) or `2` (reactive).
- `param_set` : Single scenario param_set.
- `iv_type` : Can be either `"seed"`, `"mc"` or `"fog"`.
"""
function strategy_type(strategy_idx::Int64)
    if is_reactive(strategy_idx)
        return ReactiveStrategy
    elseif is_periodic(strategy_idx)
        return PeriodicStrategy
    end
    throw(ArgumentError("Unknown mc strategy type: $strategy_type"))
end
function strategy_type(param_set::YAXArray{Float64,1}, iv_type::String)
    return strategy_type(Int64(param_set[factors = At(["$(iv_type)_strategy"])][1]))
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

"""
    revisit_cadence_mask(last_deployment, timestep, cadence)::BitVector

Locations are eligible (`true`) when never deployed to, or when at least `cadence`
timesteps have elapsed since their last deployment. `cadence <= 0` disables the
restriction entirely (all locations eligible).
"""
function revisit_cadence_mask(last_deployment, timestep, cadence)::BitVector
    cadence <= 0 && return trues(length(last_deployment))
    return (timestep .- last_deployment .>= cadence) .| (last_deployment .== 0)
end

"""
    build_state(domain, strategy, states)

Build per-target-location state for `strategy`, aligned to `unique(strategy.target_locations)`'s
own order (not `domain.loc_ids`'s order). This matters because `filter_candidate_locations`
indexes the returned state's arrays against the deduplicated target-location list directly
(see each strategy's `filter_candidate_locations`); indexing state computed in `domain.loc_ids`
order against a target list in a different order would silently return the wrong locations
whenever `strategy.target_locations` is not itself ordered consistently with `domain.loc_ids`
(e.g. a multi-share config listing locations out of domain order).
"""
function build_state(domain, strategy, states)
    targets = unique(strategy.target_locations)
    idx = indexin(targets, domain.loc_ids)
    return strategy_status(strategy, states, idx)
end

function strategy_status(::DecisionStrategy, states, idx)
    return (
        current_cover=@view(states.current_cover[idx]),
        recent_cover_losses=@view(first(states.recent_cover_losses)[idx])
    )
end

function strategy_status(::ReactiveStrategy, states, idx)
    return (
        current_cover=@view(states.current_cover[idx]),
        recent_cover_losses=@view(first(states.recent_cover_losses)[idx]),
        last_deployment=@view(states.last_deployment[idx])
    )
end

function strategy_status(::PeriodicStrategy, states, idx)
    return (
        current_cover=@view(states.current_cover[idx]),
        recent_cover_losses=@view(first(states.recent_cover_losses)[idx]),
        last_deployment=@view(states.last_deployment[idx])
    )
end
