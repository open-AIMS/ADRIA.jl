using ADRIA: YAXArray, DataCube, ZeroDataCube

"""Raw model coral cover results, proportional to k-area."""
function coral_cover(;
    n_timesteps::Int64=10,
    n_group_sizes::Int64=35,
    n_locations::Int64=5,
    n_scenarios::Int64=32
)::YAXArray{Float64,4}
    coral_cover = Factories.full_coral_cover(;
        n_timesteps, n_group_sizes, n_locations, n_scenarios
    )

    location_weights::Array{Float64,4} = repeat(
        rand(n_timesteps, n_locations, n_scenarios), 1, 1, 1, n_group_sizes
    )

    return coral_cover .* permutedims(location_weights, [1, 4, 2, 3])
end

"""Cover for each location at each timestep and scenario sums up to 1."""
function full_coral_cover(
    ;
    n_timesteps::Int64=10,
    n_group_sizes::Int64=35,
    n_locations::Int64=5,
    n_scenarios::Int64=32
)::YAXArray{Float64,4}
    cover_weights::Array{Float64,4} = rand(
        n_timesteps,
        n_group_sizes,
        n_locations,
        n_scenarios
    )

    # For each timestep and scenario the total cover of each location sums up to 1
    cover_data::Array{Float64,4} = (cover_weights ./ sum(cover_weights; dims=2))

    # This is to prevent the cover to grow above 1 due to floating point precision error
    cover_data[:, end, :, :] .=
        1 .- dropdims(
            sum(cover_data[:, 1:(end - 1), :, :]; dims=2); dims=2
        )

    return DataCube(cover_data, (:timesteps, :species, :locations, :scenarios))
end

"""Cover for all timesteps, locations and scenarios is zero"""
function zero_coral_cover(
    ;
    n_timesteps::Int64=10,
    n_group_sizes::Int64=35,
    n_locations::Int64=5,
    n_scenarios::Int64=32
)::YAXArray{Float64,4}
    return ZeroDataCube(
        (:timesteps, :species, :locations, :scenarios),
        (n_timesteps, n_group_sizes, n_locations, n_scenarios)
    )
end
