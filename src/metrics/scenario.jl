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

    s::Matrix{eltype(data)} = metric.(
        JuliennedArrays.Slices(data[timesteps=At(tf_labels)], axis_index(data, :sites))
    )

    return DataCube(s; timesteps=tf_labels, scenarios=1:size(s, 2))
end

"""
    scenario_total_cover(rs::ResultSet; kwargs...)::AbstractArray{<:Real}

Calculate the mean absolute coral for each scenario for the entire domain.
"""
function _scenario_total_cover(X::AbstractArray; kwargs...)::AbstractArray{<:Real}
    return dropdims(sum(slice_results(X; kwargs...), dims=:sites), dims=:sites)
end
function _scenario_total_cover(rs::ResultSet; kwargs...)::AbstractArray{<:Real}
    tac = total_absolute_cover(rs)
    return dropdims(sum(slice_results(tac; kwargs...), dims=:sites), dims=:sites)
end
scenario_total_cover = Metric(_scenario_total_cover, (:timesteps, :scenarios), "m²")

"""
    scenario_relative_cover(rs::ResultSet; kwargs...)::AbstractArray{<:Real}

Calculate the mean relative coral cover for each scenario for the entire domain.
"""
function _scenario_relative_cover(rs::ResultSet; kwargs...)::AbstractArray{<:Real}
    target_sites = haskey(kwargs, :sites) ? kwargs[:sites] : (:)
    target_area = sum(site_k_area(rs)[target_sites])

    return _scenario_total_cover(rs; kwargs...) ./ target_area
end
scenario_relative_cover = Metric(_scenario_relative_cover, (:timesteps, :scenarios))

"""
    scenario_relative_juveniles(data::YAXArray, coral_spec::DataFrame, k_area::AbstractVector{<:Real}; kwargs...)::AbstractArray{<:Real}
    scenario_relative_juveniles(rs::ResultSet; kwargs...)::YAXArray

Calculate the mean relative juvenile population for each scenario for the entire domain.
"""
function _scenario_relative_juveniles(
    data::YAXArray,
    coral_spec::DataFrame,
    k_area::AbstractVector{<:Real};
    kwargs...
)::AbstractArray{<:Real}
    ajuv = call_metric(absolute_juveniles, data, coral_spec; kwargs...)
    return dropdims(sum(ajuv, dims=:sites), dims=:sites) / sum(k_area)
end
function _scenario_relative_juveniles(rs::ResultSet; kwargs...)::YAXArray
    # Calculate relative domain-wide cover based on absolute values
    aj = absolute_juveniles(rs)
    return dropdims(sum(aj, dims=:sites), dims=:sites) ./ sum(site_k_area(rs))
end
scenario_relative_juveniles = Metric(_scenario_relative_juveniles, (:timesteps, :scenarios))

"""
    scenario_absolute_juveniles(data::YAXArray, coral_spec::DataFrame, k_area::AbstractVector{<:Real}; kwargs...)::AbstractArray{<:Real}
    scenario_absolute_juveniles(rs::ResultSet; kwargs...)::AbstractArray{<:Real}

Calculate the mean absolute juvenile population for each scenario for the entire domain.
"""
function _scenario_absolute_juveniles(
    data::YAXArray,
    coral_spec::DataFrame,
    k_area::AbstractVector{<:Real};
    kwargs...
)::AbstractArray{<:Real}
    juv = call_metric(absolute_juveniles, data, coral_spec; kwargs...)
    return dropdims(sum(juv, dims=:sites), dims=:sites) / sum(k_area)
end
function _scenario_absolute_juveniles(rs::ResultSet; kwargs...)::AbstractArray{<:Real}
    # Calculate relative domain-wide cover based on absolute values
    return dropdims(sum(absolute_juveniles(rs), dims=:sites), dims=:sites)
end
scenario_absolute_juveniles = Metric(_scenario_absolute_juveniles, (:timesteps, :scenarios))

"""
    scenario_juvenile_indicator(data::YAXArray, coral_spec::DataFrame, k_area::AbstractVector{<:Real}; kwargs...)::AbstractArray{<:Real}
    scenario_juvenile_indicator(rs::ResultSet; kwargs...)::AbstractArray{<:Real}

Determine juvenile indicator ∈ [0, 1], where 1 indicates maximum mean juvenile density (51.8) has been achieved.
"""
function _scenario_juvenile_indicator(
    data::YAXArray,
    coral_spec::DataFrame,
    k_area::AbstractVector{<:Real};
    kwargs...
)::AbstractArray{<:Real}
    juv = call_metric(juvenile_indicator, data, coral_spec, k_area; kwargs...)
    return dropdims(mean(juv, dims=:sites), dims=:sites) / sum(k_area)
end
function _scenario_juvenile_indicator(rs::ResultSet; kwargs...)::AbstractArray{<:Real}
    return dropdims(mean(juvenile_indicator(rs), dims=:sites), dims=:sites)
end
scenario_juvenile_indicator = Metric(_scenario_juvenile_indicator, (:timesteps, :scenarios))

"""
    scenario_asv(sv::YAXArray; kwargs...)::AbstractArray{<:Real}
    scenario_asv(rs::ResultSet; kwargs...)::AbstractArray{<:Real}

Calculate the mean absolute shelter volumes for each scenario for the entire domain.
"""
function _scenario_asv(sv::YAXArray; kwargs...)::AbstractArray{<:Real}
    sv_sliced = slice_results(sv; kwargs...)
    return dropdims(sum(sv_sliced, dims=:sites), dims=:sites)
end
function _scenario_asv(rs::ResultSet; kwargs...)::AbstractArray{<:Real}
    return _scenario_asv(rs.outcomes[:absolute_shelter_volume]; kwargs...)
end
scenario_asv = Metric(_scenario_asv, (:timesteps, :scenarios), "m³/m²")

"""
    scenario_rsv(sv::YAXArray; kwargs...)::AbstractArray{<:Real}
    scenario_rsv(rs::ResultSet; kwargs...)::AbstractArray{<:Real}

Calculate the mean relative shelter volumes for each scenario for the entire domain.
"""
function _scenario_rsv(sv::YAXArray; kwargs...)::AbstractArray{<:Real}
    sv_sliced = slice_results(sv; kwargs...)
    return dropdims(mean(sv_sliced, dims=:sites), dims=:sites)
end
function _scenario_rsv(rs::ResultSet; kwargs...)::AbstractArray{<:Real}
    return _scenario_rsv(rs.outcomes[:relative_shelter_volume]; kwargs...)
end
scenario_rsv = Metric(_scenario_rsv, (:timesteps, :scenarios))

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
scenario_evenness = Metric(_scenario_evenness, (:timesteps, :scenarios))

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
        outcomes=to_symbol.(metrics),
    )

    for (i, metric) in enumerate(metrics)
        scen_outcomes[:, :, i] = metric(rs)
    end

    return scen_outcomes
end
