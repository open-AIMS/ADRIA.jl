using ADRIA: YAXArray, DataCube, ZeroDataCube

"""Raw model coral cover results, proportional to k-area."""
function coral_cover(;
    n_timesteps::Int64=10,
    n_groups::Int64=5,
    n_sizes::Int64=7,
    n_locations::Int64=5,
    n_scenarios::Int64=32
)::YAXArray{Float64,5}
    coral_cover = MockData.full_coral_cover(;
        n_timesteps, n_groups, n_sizes, n_locations, n_scenarios
    )

    location_weights::Array{Float64,5} = reshape(
        rand(n_timesteps, n_locations, n_scenarios),
        (n_timesteps, 1, 1, n_locations, n_scenarios)
    )

    return coral_cover .* location_weights
end

"""Cover for each location at each timestep and scenario sums up to 1."""
function full_coral_cover(
    ;
    n_timesteps::Int64=10,
    n_groups::Int64=5,
    n_sizes::Int64=7,
    n_locations::Int64=5,
    n_scenarios::Int64=32
)::YAXArray{Float64,5}
    cover_weights::Array{Float64,5} = rand(
        n_timesteps,
        n_groups,
        n_sizes,
        n_locations,
        n_scenarios
    )

    # For each timestep and scenario the total cover of each location sums up to 1
    cover_data::Array{Float64,5} =
        (
            cover_weights ./ sum(cover_weights; dims=(2, 3))
        ) .* 0.999999

    return DataCube(cover_data, (:timesteps, :groups, :sizes, :locations, :scenarios))
end

"""Cover for all timesteps, locations and scenarios is zero"""
function zero_coral_cover(
    ;
    n_timesteps::Int64=10,
    n_groups::Int64=5,
    n_sizes::Int64=7,
    n_locations::Int64=5,
    n_scenarios::Int64=32
)::YAXArray{Float64,5}
    return ZeroDataCube(
        (:timesteps, :groups, :sizes, :locations, :scenarios),
        (n_timesteps, n_groups, n_sizes, n_locations, n_scenarios)
    )
end
