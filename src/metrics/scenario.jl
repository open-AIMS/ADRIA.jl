"""Scenario-level summaries.

Collapses the `species` and `site` dimension so it is just trajectories over time for each scenario.

Note: The difference between this and `temporal` metrics is that these methods
keep the scenario dimension.

TODO: Produce summary stats. Currently returns just the mean.
"""


function scenario_total_cover(data::NamedDimsArray, areas; kwargs...)
    squash = (:sites, :reps)

    sites = haskey(kwargs, :sites) ? kwargs[:sites] : (:)
    tac = call_metric(total_cover, data, areas[sites]; kwargs...)

    return dropdims(mean(tac, dims=squash), dims=squash)
end
function scenario_total_cover(rs::ResultSet; kwargs...)
    return scenario_total_cover(rs.raw, rs.site_area; kwargs...)
end


function scenario_juveniles(data::NamedDimsArray; kwargs...)
    squash = (:sites, :reps)

    juv = call_metric(juveniles, data; kwargs...)
    return dropdims(mean(juv, dims=squash), dims=squash)
end
function scenario_juveniles(rs::ResultSet; kwargs...)
    return scenario_juveniles(rs.raw; kwargs...)
end


function scenario_shelter_volume(data::NamedDimsArray, inputs; kwargs...)
    squash = (:sites, :reps)
    sv = call_metric(shelter_volume, data, inputs; kwargs...)

    return dropdims(mean(sv, dims=squash), dims=squash)
end
function scenario_shelter_volume(rs::ResultSet)
    return scenario_shelter_volume(rs.raw, rs.inputs)
end
