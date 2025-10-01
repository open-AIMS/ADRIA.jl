using Bootstrap
using Random

"""Scenario-level summaries.

Note: Aggregates across the `site` dimension so trajectories over time for each scenario are
      returned. The difference between these and `temporal` metrics is that these methods
      keep the scenario dimension.

TODO: Produce summary stats. Currently returns just the mean.
"""

"""
    scenario_trajectory(data::AbstractArray; metric=mean)::YAXArray{<:Real}

Produce scenario trajectories using the provided metric/aggregation function.

# Arguments
- `data` : Results to aggregate
- `metric` : Function or Callable used to summarize data

# Returns
Matrix[timesteps ⋅ scenarios]
"""
function scenario_trajectory(data::AbstractArray; metric=mean)::YAXArray{<:Real}
    tf_labels = axis_labels(data, :timesteps)

    s::Matrix{eltype(data)} =
        metric.(
            JuliennedArrays.Slices(
                data[timesteps=At(tf_labels)], axis_index(data, :locations)
            )
        )

    return DataCube(s; timesteps=tf_labels, scenarios=1:size(s, 2))
end

"""
    scenario_total_cover(rs::ResultSet; kwargs...)::AbstractArray{<:Real}

Calculate the mean absolute coral for each scenario for the entire domain.

# Arguments
- `tac` : Total absolute cover
- `rs` : ResultSet
"""
function _scenario_total_cover(tac::AbstractArray; kwargs...)::AbstractArray{<:Real}
    return dropdims(sum(slice_results(tac; kwargs...); dims=:locations); dims=:locations)
end
function _scenario_total_cover(rs::ResultSet; kwargs...)::AbstractArray{<:Real}
    tac = total_absolute_cover(rs)
    return _scenario_total_cover(tac::AbstractArray; kwargs...)
end
scenario_total_cover = Metric(
    _scenario_total_cover,
    (:timesteps, :locations, :scenarios),
    (:timesteps, :scenarios),
    "Cover",
    IS_NOT_RELATIVE,
    UNIT_AREA
)

"""
    scenario_relative_cover(rs::ResultSet; kwargs...)::AbstractArray{<:Real}

Calculate the mean relative coral cover for each scenario for the entire domain.
"""
function _scenario_relative_cover(rs::ResultSet; kwargs...)::AbstractArray{<:Real}
    target_sites = haskey(kwargs, :locations) ? kwargs[:locations] : (:)
    target_area = sum(loc_k_area(rs)[target_sites])

    return _scenario_total_cover(rs; kwargs...) ./ target_area
end
scenario_relative_cover = Metric(
    _scenario_relative_cover,
    (:timesteps, :scenarios, :locations),
    (:timesteps, :scenarios),
    "Relative Cover",
    IS_RELATIVE
)

"""
    scenario_relative_juveniles(X::YAXArray{<:Real,3}, coral_spec::DataFrame, k_area::AbstractVector{<:Real}; kwargs...)::AbstractArray{<:Real}
    scenario_relative_juveniles(rs::ResultSet; kwargs...)::YAXArray

Calculate the mean relative juvenile population for each scenario for the entire domain.

# Arguments
- `X` : Raw data for a single scenario.
- `rs` : Resultset.
- `coral_spec` : Coral spec DataFrame.
- `k_area` : K_area.

# Examples
```
num_scens = 2^5
scens = ADRIA.sample(dom, num_scens)

_coral_spec = ADRIA.to_coral_spec(scens[1,:])
_k_area = loc_k_area(dom)

# X contains raw coral cover results for a single scenario
ADRIA.metrics.scenario_relative_juveniles(X, _coral_spec, _k_area)
```
"""
function _scenario_relative_juveniles(
    aj::YAXArray{<:Real,3},
    k_area::AbstractVector{<:Real}
)::AbstractArray{<:Real}
    return dropdims(sum(aj; dims=:locations); dims=:locations) ./ sum(k_area)
end
function _scenario_relative_juveniles(rs::ResultSet; kwargs...)::YAXArray
    # Calculate relative domain-wide cover based on absolute values
    aj = absolute_juveniles(rs)
    return _scenario_relative_juveniles(aj, loc_k_area(rs))
end
scenario_relative_juveniles = Metric(
    _scenario_relative_juveniles,
    (:timesteps, :locations, :scenarios),
    (:timesteps, :scenarios),
    "Relative Juveniles",
    IS_RELATIVE
)

"""
    scenario_absolute_juveniles(data::YAXArray, coral_spec::DataFrame, k_area::AbstractVector{<:Real}; kwargs...)::AbstractArray{<:Real}
    scenario_absolute_juveniles(rs::ResultSet; kwargs...)::AbstractArray{<:Real}

Calculate the mean absolute juvenile population for each scenario for the entire domain.

# Arguments
- `aj` : Raw data for a single scenario.
- `k_area` : K_area.
- `rs` : Resultset.
"""
function _scenario_absolute_juveniles(
    aj::YAXArray
)::AbstractArray{<:Real}
    return dropdims(sum(aj; dims=:locations); dims=:locations)
end
function _scenario_absolute_juveniles(rs::ResultSet; kwargs...)::AbstractArray{<:Real}
    aj = absolute_juveniles(rs)
    return _scenario_absolute_juveniles(aj)
end
scenario_absolute_juveniles = Metric(
    _scenario_absolute_juveniles,
    (:timesteps, :locations, :scenarios),
    (:timesteps, :scenarios),
    "Number of Juveniles",
    IS_NOT_RELATIVE,
    UNIT_AREA
)

"""
    scenario_juvenile_indicator(data::YAXArray, coral_spec::DataFrame, k_area::AbstractVector{<:Real}; kwargs...)::AbstractArray{<:Real}
    scenario_juvenile_indicator(rs::ResultSet; kwargs...)::AbstractArray{<:Real}

Determine juvenile indicator ∈ [0, 1], where 1 indicates maximum mean juvenile density (51.8) has been achieved.

# Arguments
- `ji` : Juvenile Indicator for each location.
- `rs` : Resultset.
"""
function _scenario_juvenile_indicator(
    ji::YAXArray{<:Real,3}
)::AbstractArray{<:Real}
    return dropdims(mean(ji; dims=:locations); dims=:locations)
end
function _scenario_juvenile_indicator(rs::ResultSet; kwargs...)::AbstractArray{<:Real}
    ji = juvenile_indicator(rs)
    return _scenario_juvenile_indicator(ji)
end
scenario_juvenile_indicator = Metric(
    _scenario_juvenile_indicator,
    (:timesteps, :locations, :scenarios),
    (:timesteps, :scenarios),
    "Juvenile Indicator",
    IS_NOT_RELATIVE
)

"""
    scenario_asv(sv::YAXArray; kwargs...)::AbstractArray{<:Real}
    scenario_asv(rs::ResultSet; kwargs...)::AbstractArray{<:Real}

Calculate the mean absolute shelter volumes for each scenario for the entire domain.

# Arguments
- `asv` : Absolute shelter volume.
- `rs` : Resultset.
"""
function _scenario_asv(asv::YAXArray; kwargs...)::AbstractArray{<:Real}
    sv_sliced = slice_results(asv; kwargs...)
    return dropdims(sum(sv_sliced; dims=:locations); dims=:locations)
end
function _scenario_asv(rs::ResultSet; kwargs...)::AbstractArray{<:Real}
    return _scenario_asv(rs.outcomes[:absolute_shelter_volume]; kwargs...)
end
scenario_asv = Metric(
    _scenario_asv,
    (:timesteps, :locations, :scenarios),
    (:timesteps, :scenarios),
    "Volume",
    IS_NOT_RELATIVE,
    "$UNIT_VOLUME/$UNIT_AREA"
)

"""
    scenario_rsv(sv::YAXArray; kwargs...)::AbstractArray{<:Real}
    scenario_rsv(rs::ResultSet; kwargs...)::AbstractArray{<:Real}

Calculate the mean relative shelter volumes for each scenario for the entire domain.
"""
function _scenario_rsv(sv::YAXArray; kwargs...)::AbstractArray{<:Real}
    sv_sliced = slice_results(sv; kwargs...)
    return dropdims(mean(sv_sliced; dims=:locations); dims=:locations)
end
function _scenario_rsv(rs::ResultSet; kwargs...)::AbstractArray{<:Real}
    return _scenario_rsv(rs.outcomes[:relative_shelter_volume]; kwargs...)
end
scenario_rsv = Metric(
    _scenario_rsv,
    (:timesteps, :locations, :scenarios),
    (:timesteps, :scenarios),
    "Relative Volume",
    IS_RELATIVE
)

"""
    scenario_evenness(ev::YAXArray; kwargs...)::AbstractArray{<:Real}
    scenario_evenness(rs::ResultSet; kwargs...)::AbstractArray{<:Real}

Calculate the mean coral evenness for each scenario for the entire domain.
"""
function _scenario_evenness(ev::YAXArray; kwargs...)::AbstractArray{<:Real}
    ev_sliced = slice_results(ev; kwargs...)
    return scenario_trajectory(ev_sliced)
    # return dropdims(mean(ev_sliced, dims=:sites), dims=:sites)
end
function _scenario_evenness(rs::ResultSet; kwargs...)::AbstractArray{<:Real}
    return _scenario_evenness(rs.outcomes[:coral_evenness]; kwargs...)
end
scenario_evenness = Metric(
    _scenario_evenness,
    (:timesteps, :locations, :scenarios),
    (:timesteps, :scenarios),
    "Evenness Indicator",
    IS_NOT_RELATIVE
)

"""
    scenario_outcomes(rs::ResultSet, metrics::Vector{Metric})::YAXArray

Get outcomes for a given list of metrics and a result set.

# Arguments
- `rs` : ResultSet
- `metrics` : Vector of scenario Metrics (the ones that start with `scenario_`)

# Returns
YAXArray with (:timesteps, :scenarios, :outcomes)

# Examples
```
metrics::Vector{ADRIA.metrics.Metric} = [
    ADRIA.metrics.scenario_total_cover,
    ADRIA.metrics.scenario_asv,
    ADRIA.metrics.scenario_absolute_juveniles,
]

# 3-dimensional Array of outcomes
outcomes = ADRIA.metrics.scenario_outcomes(rs, metrics)
```
"""
function scenario_outcomes(rs::ResultSet, metrics::Vector{<:Metric})::YAXArray
    n_scenarios = size(rs.inputs, 1)

    scen_outcomes = ZeroDataCube(;
        timesteps=timesteps(rs),
        scenarios=1:n_scenarios,
        outcomes=to_symbol.(metrics)
    )

    for (i, metric) in enumerate(metrics)
        scen_outcomes[:, :, i] = metric(rs)
    end

    return scen_outcomes
end
