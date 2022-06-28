"""Scenario-level summaries.

Collapses the `species` and `site` dimension so it is just trajectories over time for each scenario.

Note: The difference between this and `temporal` metrics is that these methods
keep the scenario dimension.

TODO: Produce summary stats. Currently returns just the mean.
"""


function scenario_total_cover(data, areas, timesteps=(:))
    target_dims = (:sites, :reps)
    return dropdims(mean(total_cover(data, areas), dims=target_dims), dims=target_dims)[timesteps=timesteps]
end
function scenario_total_cover(rs::ResultSet)
    return scenario_total_cover(rs.raw, rs.site_area)
end


function scenario_juveniles(data, timesteps=(:))
    target_dims = (:sites, :reps)
    return dropdims(mean(juveniles(data), dims=target_dims), dims=target_dims)[timesteps=timesteps]
end
function scenario_juveniles(rs::ResultSet)
    return scenario_juveniles(rs.raw)
end


function scenario_shelter_volume(data, inputs, timesteps=(:))
    target_dims = (:sites, :reps)
    return dropdims(mean(shelter_volume(data, inputs), dims=target_dims), dims=target_dims)[timesteps=timesteps]
end
function scenario_shelter_volume(rs::ResultSet)
    return scenario_shelter_volume(rs.raw, rs.inputs)
end
