"""
    PeriodicStrategy <: DecisionStrategy

Deploys at regular intervals starting from a specified year.

Locations are filtered by `revisit_cadence`: a location deployed to within the last
`revisit_cadence` timesteps is excluded from candidates, forcing spatial spread of
deployments. If cadence filtering drops the eligible candidate count below
`min_locations`, the pool is backfilled with the excluded locations that were deployed
to longest ago (soonest to become eligible again), breaking ties by stable-sort order
(i.e. `target_locations`'s original order), until `min_locations` is met or all target
locations are included.

Note: backfill operates on the aggregate `target_locations` pool, not per deployment
share. A backfilled location — one still within its own cadence window — is included in
the returned candidate list and so is legitimately eligible for selection by *any* share
that contains it, not just a share short of its own `min_iv_locations`. This is the
intended effect of aggregate-pool backfill, not a bug. Because `min_iv_locations`/
`mc_min_iv_locations` is enforced per deployment share downstream (in `scenario.jl`, for
seed/mc only — fog has no per-share loop), aggregate backfill does not guarantee any
individual share clears its own minimum; that is a pre-existing limitation of
`min_iv_locations` in the multi-share case, unrelated to and not worsened by cadence.

# Fields
- `target_locations::Vector{String}`: Optional list of location IDs to target
- `start_year::Int64`: Year to begin deployments (relative to simulation start)
- `duration::Int64`: Number of years to deploy for
- `frequency::Int64`: Frequency (in years) between deployments (0 = deploy once)
- `decision_years::Vector{Bool}`: Precomputed vector indicating deployment years
- `revisit_cadence::Int64`: Minimum timesteps before a location can be re-selected (0 = no restriction)
- `min_locations::Int64`: Minimum candidate locations to return; triggers backfill if cadence filtering drops below this
"""
struct PeriodicStrategy <: DecisionStrategy
    target_locations::Vector{String}
    start_year::Int64
    duration::Int64
    frequency::Union{Float64,Int64}  # in years
    decision_years::Vector{Bool}  # precomputed decision years
    revisit_cadence::Int64
    min_locations::Int64

    function PeriodicStrategy(
        target_locations::Vector{String},
        start_year::Int64,
        duration::Int64,
        frequency::Union{Float64,Int64},
        timeframe::Int64,
        revisit_cadence::Int64,
        min_locations::Int64
    )
        decision_years = decision_frequency(start_year, timeframe, duration, frequency)
        return new(
            target_locations, start_year, duration, frequency, decision_years,
            revisit_cadence, min_locations
        )
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

Return target locations for periodic deployment strategy, filtered by revisit cadence.

Periodic strategies deploy to all specified target locations whenever it is a decision
year, subject to `revisit_cadence`: locations deployed to within the last
`revisit_cadence` timesteps are excluded, unless doing so would drop the candidate count
below `min_locations`, in which case the pool is backfilled with the longest-idle
excluded locations (see struct docstring).

# Arguments
- `strategy`: PeriodicStrategy containing target locations, revisit_cadence and min_locations
- `timestep`: Current simulation timestep
- `state`: NamedTuple with `last_deployment::Vector{Int64}`, aligned to
  `unique(strategy.target_locations)`'s own order (see `build_state`)

# Returns
Vector of location IDs eligible for deployment
"""
function filter_candidate_locations(
    strategy::PeriodicStrategy,
    timestep::Int64,
    state::Union{Nothing,NamedTuple}
)::Vector{String}
    if !is_decision_year(strategy, timestep)
        return String[]
    end

    targets = unique(strategy.target_locations)

    if strategy.revisit_cadence <= 0 || isnothing(state)
        return targets
    end

    eligible = revisit_cadence_mask(
        state.last_deployment, timestep, strategy.revisit_cadence
    )

    if count(eligible) >= strategy.min_locations
        return targets[eligible]
    end

    # Backfill with excluded locations, longest-idle first (stable sort preserves
    # `targets`'s original order as the tie-break).
    excluded_idx = findall(.!eligible)
    idle_time = timestep .- state.last_deployment[excluded_idx]
    order = sortperm(idle_time; rev=true)
    backfill_idx = excluded_idx[order]

    n_needed = strategy.min_locations - count(eligible)
    n_backfill = min(n_needed, length(backfill_idx))

    return vcat(targets[eligible], targets[backfill_idx[1:n_backfill]])
end
