"""Scenario-level summaries.

Note: Aggregates across the `site` dimension so trajectories over time for each scenario are returned.
      The difference between these and `temporal` metrics is that these methods keep the scenario dimension.

TODO: Produce summary stats. Currently returns just the mean.
"""


"""
    scenario_total_cover(data::NamedDimsArray, areas; kwargs...)

Calculate the cluster-wide total absolute coral cover for each individual scenario.
"""
function scenario_total_cover(data::NamedDimsArray, areas; kwargs...)
    sites = haskey(kwargs, :sites) ? kwargs[:sites] : (:)
    tac = call_metric(_total_absolute_cover, data, areas[sites]; kwargs...)

    return dropdims(sum(tac, dims=:sites), dims=:sites)
end
function scenario_total_cover(rs::ResultSet; kwargs...)
    return scenario_total_cover(rs.outcomes[:relative_cover], rs.site_area; kwargs...)
end


"""
    scenario_juveniles(data::NamedDimsArray; kwargs...)

Calculate the cluster-wide juvenile population for each individual scenario.
"""
function scenario_juveniles(data::NamedDimsArray; kwargs...)
    juv = call_metric(juveniles, data; kwargs...)
    return dropdims(sum(juv, dims=:sites), dims=:sites)
end
function scenario_juveniles(rs::ResultSet; kwargs...)
    return scenario_juveniles(rs.raw; kwargs...)
end


"""
    scenario_shelter_volume(data::NamedDimsArray, inputs; kwargs...)

Calculate the cluster-wide absolute shelter volume for each individual scenario.
"""
function scenario_shelter_volume(data::NamedDimsArray, area, max_cover, inputs; kwargs...)
    sv = call_metric(absolute_shelter_volume, data, area, max_cover, inputs; kwargs...)
    return dropdims(sum(sv, dims=:sites), dims=:sites)
end

function scenario_shelter_volume(sv::NamedDimsArray; kwargs...)
    sv_sliced = slice_results(sv; kwargs...)
    return dropdims(sum(sv_sliced, dims=:sites), dims=:sites)
end
function scenario_shelter_volume(rs::ResultSet; kwargs...)
    return scenario_shelter_volume(rs.outcomes[:absolute_shelter_volume]; kwargs...)
end
