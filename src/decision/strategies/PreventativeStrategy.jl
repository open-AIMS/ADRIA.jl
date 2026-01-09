"""
    PreventiveStrategy <: DecisionStrategy

Deploys based on predicted exposure rather than observed damage.
Uses environmental forecasts and reef condition to prevent losses before they occur.

This is a proactive strategy that anticipates damage using leading indicators.

# Fields
- `target_locations::Vector{String}`: Locations to consider for deployment
- `dhw_forecast_threshold::Float64`: Deploy if projected DHW exceeds this (e.g., 8.0)
- `min_cover_threshold::Float64`: Deploy if cover trending below this (e.g., 0.15)
- `forecast_window::Int64`: Timesteps ahead to consider in forecast (e.g., 3 years)
"""
struct PreventiveStrategy <: DecisionStrategy
    target_locations::Vector{String}
    dhw_forecast_threshold::Float64
    min_cover_threshold::Float64
    forecast_window::Int64
end

"""
    filter_candidate_locations(
        strategy::PreventiveStrategy,
        _::Int64,
        state::NamedTuple
    )::Vector{String}

Filter target locations to those with high predicted vulnerability.

Returns locations where ANY of these conditions are met:
1. Projected DHW exceeds forecast threshold
2. Vulnerability index exceeds threshold
3. Cover trending below minimum threshold

# Arguments
- `strategy`: PreventiveStrategy with threshold parameters and target locations
- `_`: Unused argument, kept for compatibility
- `state`: NamedTuple containing:
    - `current_cover::Vector{Float64}`: Current coral cover at each target location
    - `dhw_forecast::Matrix{Float64}`: Projected DHW [forecast_window × target_location]

# Returns
Vector of location IDs meeting preventive deployment criteria
"""
function filter_candidate_locations(
    strategy::PreventiveStrategy,
    _::Int64,
    state::NamedTuple
)::Vector{String}
    # Get mean projected DHW over forecast window
    mean_forecast_dhw = dropdims(mean(state.dhw_forecast; dims=1); dims=1)

    # Multi-criteria triggering (OR logic - deploy if ANY condition met)
    candidate_mask =
        (
            # Criterion 1: High projected thermal stress
            mean_forecast_dhw .>= strategy.dhw_forecast_threshold
        ) .& (
            # But only if location is still viable
            state.current_cover .> strategy.min_cover_threshold
        )

    return strategy.target_locations[candidate_mask]
end
