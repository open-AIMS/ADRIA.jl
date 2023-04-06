"""Scenario-level summaries.

Note: Aggregates across the `location` dimension so trajectories over time for each scenario are returned.
      The difference between these and `temporal` metrics is that these methods keep the scenario dimension.

TODO: Produce summary stats. Currently returns just the mean.
"""


"""
    scenario_total_cover(rs::ResultSet; kwargs...)

Calculate the cluster-wide total absolute coral cover for each scenario.
"""
function _scenario_total_cover(rs::ResultSet; kwargs...)
    return dropdims(sum(slice_results(total_absolute_cover(rs); kwargs...), dims=:locations), dims=:locations)
end
scenario_total_cover = Metric(_scenario_total_cover, (:timesteps, :scenarios), "m²")


"""
    scenario_relative_cover(rs::ResultSet; kwargs...)

Calculate the cluster-wide relative coral cover for each scenario.
"""
function _scenario_relative_cover(rs::ResultSet; kwargs...)
    target_locations = haskey(kwargs, :locations) ? kwargs[:locations] : (:)
    target_area = sum(((rs.location_max_coral_cover./100.0).*rs.location_area)[target_locations])

    return _scenario_total_cover(rs; kwargs...) ./ target_area
end
scenario_relative_cover = Metric(_scenario_relative_cover, (:timesteps, :scenarios))


"""
    scenario_juveniles(data::NamedDimsArray; kwargs...)

Calculate the cluster-wide relative juvenile population for individual scenarios.

!!! warning DEPRECATED.
    This function is now deprecated. Use `scenario_relative_juveniles()` instead.
"""
function _scenario_juveniles(data::NamedDimsArray, coral_spec::DataFrame, area::AbstractVector{<:Real}; kwargs...)
    @warn "`scenario_juveniles()` is deprecated and will be removed in future versions. Use `scenario_relative_juveniles()` instead."
    return _scenario_rjuves(data, coral_spec, area; kwargs...)
end
function _scenario_juveniles(rs::ResultSet; kwargs...)
    @warn "`scenario_juveniles()` is deprecated and will be removed in future versions. Use `scenario_relative_juveniles()` instead."
    return scenario_rjuves(rs)
end
scenario_juveniles = Metric(_scenario_juveniles, (:timesteps, :scenario))


"""
    scenario_rjuves(data::NamedDimsArray, coral_spec::DataFrame, area::AbstractVector{<:Real}; kwargs...)

Calculate the relative cluster-wide juvenile population for individual scenarios.
"""
function _scenario_relative_juveniles(data::NamedDimsArray, coral_spec::DataFrame, area::AbstractVector{<:Real}; kwargs...)::NamedDimsArray
    ajuv = call_metric(absolute_juveniles, data, coral_spec; kwargs...)
    return dropdims(sum(ajuv, dims=:locations), dims=:locations) / sum(area)
end
function _scenario_relative_juveniles(rs::ResultSet; kwargs...)::NamedDimsArray
    # Calculate relative domain-wide cover based on absolute values
    # Note: element-wise operator (./ sum(...)) required to maintain NamedDimArray
    #       otherwise it will get returned as an AxisKey array for some reason
    return dropdims(sum(absolute_juveniles(rs), dims=:locations), dims=:locations) ./ sum(rs.location_area)
end
scenario_relative_juveniles = Metric(_scenario_relative_juveniles, (:timesteps, :scenario))


"""
    scenario_absolute_juveniles(data::NamedDimsArray, coral_spec::DataFrame, area::AbstractVector{<:Real}; kwargs...)
    scenario_absolute_juveniles(rs::ResultSet; kwargs...)::AbstractArray

Calculate the absolute cluster-wide juvenile population for individual scenarios.
"""
function _scenario_absolute_juveniles(data::NamedDimsArray, coral_spec::DataFrame, area::AbstractVector{<:Real}; kwargs...)
    juv = call_metric(absolute_juveniles, data, coral_spec; kwargs...)
    return dropdims(sum(juv, dims=:locations), dims=:locations) / sum(area)
end
function _scenario_absolute_juveniles(rs::ResultSet; kwargs...)::AbstractArray
    # Calculate relative domain-wide cover based on absolute values
    return dropdims(sum(absolute_juveniles(rs), dims=:locations), dims=:locations)
end
scenario_absolute_juveniles = Metric(_scenario_absolute_juveniles, (:timesteps, :scenario))


"""
    _scenario_juvenile_indicator(data::NamedDimsArray, coral_spec::DataFrame, area::V, k_area::V; kwargs...) where {V<:AbstractVector{<:Real}}
    _scenario_juvenile_indicator(rs::ResultSet; kwargs...)::AbstractArray

Determine juvenile indicator ∈ [0, 1], where 1 indicates maximum mean juvenile density (51.8) has been achieved.
"""
function _scenario_juvenile_indicator(data::NamedDimsArray, coral_spec::DataFrame, area::V, k_area::V; kwargs...) where {V<:AbstractVector{<:Real}}
    juv = call_metric(juvenile_indicator, data, coral_spec, area, k_area; kwargs...)
    return dropdims(mean(juv, dims=:locations), dims=:locations) / sum(area)
end
function _scenario_juvenile_indicator(rs::ResultSet; kwargs...)::AbstractArray
    return dropdims(mean(juvenile_indicator(rs), dims=:locations), dims=:locations)
end
scenario_juvenile_indicator = Metric(_scenario_juvenile_indicator, (:timesteps, :scenario))


"""
    scenario_asv(sv::NamedDimsArray; kwargs...)
    scenario_asv(rs::ResultSet; kwargs...)

Calculate the cluster-wide absolute shelter volume for each scenario.
"""
function _scenario_asv(sv::NamedDimsArray; kwargs...)
    sv_sliced = slice_results(sv; kwargs...)
    return dropdims(sum(sv_sliced, dims=:locations), dims=:locations)
end
function _scenario_asv(rs::ResultSet; kwargs...)
    return _scenario_asv(rs.outcomes[:absolute_shelter_volume]; kwargs...)
end
scenario_asv = Metric(_scenario_asv, (:scenario, :timesteps), "m³/m²")


"""
    scenario_rsv(sv::NamedDimsArray; kwargs...)
    scenario_rsv(rs::ResultSet; kwargs...)

Calculate the cluster-wide mean relative shelter volumes for each scenario.
"""
function _scenario_rsv(sv::NamedDimsArray; kwargs...)
    sv_sliced = slice_results(sv; kwargs...)
    return dropdims(mean(sv_sliced, dims=:locations), dims=:locations)
end
function _scenario_rsv(rs::ResultSet; kwargs...)
    return _scenario_rsv(rs.outcomes[:relative_shelter_volume]; kwargs...)
end
scenario_rsv = Metric(_scenario_rsv, (:scenario, :timesteps))
