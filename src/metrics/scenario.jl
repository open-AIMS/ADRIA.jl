"""Scenario-level summaries.

Note: Aggregates across the `site` dimension so trajectories over time for each scenario are returned.
      The difference between these and `temporal` metrics is that these methods keep the scenario dimension.

TODO: Produce summary stats. Currently returns just the mean.
"""


"""
    scenario_total_cover(rs::ResultSet; kwargs...)

Calculate the cluster-wide total absolute coral cover for each individual scenario.
"""
function scenario_total_cover(rs::ResultSet; kwargs...)
    return dropdims(sum(slice_results(total_absolute_cover(rs); kwargs...), dims=:sites), dims=:sites)
end


"""
    scenario_relative_cover(rs::ResultSet; kwargs...)

Calculate the cluster-wide mean relative coral cover for each individual scenario.
"""
function scenario_relative_cover(rs::ResultSet; kwargs...)
    return dropdims(mean(slice_results(relative_cover(rs); kwargs...), dims=:sites), dims=:sites)
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
    scenario_asv(sv::NamedDimsArray; kwargs...)
    scenario_asv(rs::ResultSet; kwargs...)

Calculate the cluster-wide absolute shelter volume for each individual scenario.
"""
function scenario_asv(sv::NamedDimsArray; kwargs...)
    sv_sliced = slice_results(sv; kwargs...)
    return dropdims(sum(sv_sliced, dims=:sites), dims=:sites)
end
function scenario_asv(rs::ResultSet; kwargs...)
    return scenario_asv(rs.outcomes[:absolute_shelter_volume]; kwargs...)
end


"""
    scenario_rsv(sv::NamedDimsArray; kwargs...)
    scenario_rsv(rs::ResultSet; kwargs...)

Calculate the cluster-wide mean relative shelter volumes for each individual scenario.
"""
function scenario_rsv(sv::NamedDimsArray; kwargs...)
    sv_sliced = slice_results(sv; kwargs...)
    return dropdims(mean(sv_sliced, dims=:sites), dims=:sites)
end
function scenario_rsv(rs::ResultSet; kwargs...)
    return scenario_rsv(rs.outcomes[:relative_shelter_volume]; kwargs...)
end
