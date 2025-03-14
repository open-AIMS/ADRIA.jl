"""
    distance_selected_locations(selected_locations1::Vector{String}, selected_locations2::Vector{String}, loc_data::DataFrame)::Float64

Given two vectors of location ids, it computes the average pairwise distance between locations for each possible pair
of locations.
It approximates each location's position as the location's polygon centroid.
The Haversine formula is used to compute the distance.

# Arguments
- `selected_locations1` : Vector with location id (reef_siteid or UNIQUE_ID) of the first set of locations.
- `selected_locations2` : Vector with location id (reef_siteid or UNIQUE_ID) of the second set of locations.
- `loc_data` : Location data for all locations.

# Returns
Average pairwise distance between locations at selected_locations1 and selected_locations2.
"""
function distance_selected_locations(
        selected_locations1::Vector{String},
        selected_locations2::Vector{String},
        loc_data::DataFrame
)::Float64
    reefs1 = filter_locations(loc_data, selected_locations1)
    reefs2 = filter_locations(loc_data, selected_locations2)

    reefs1_coord = ADRIA.centroids(reefs1)
    reefs2_coord = ADRIA.centroids(reefs2)

    pairwise_combination = Iterators.product(reefs1_coord, reefs2_coord)
    distance = mapreduce(coords -> Distances.haversine(coords[1], coords[2]), +, pairwise_combination)
    return distance / length(pairwise_combination)
end

"""
    filter_locations(loc_data::DataFrame, selected_locations::Vector{String})

Filter locations based on the reef_siteid or UNIQUE_ID.
# Arguments
- `loc_data` : Location data for all locations.
- `selected_locations` : Vector with reef_siteid of the first set of locations.
"""
function filter_locations(loc_data::DataFrame, selected_locations::Vector{String})
    if startswith(selected_locations[1], "reef_")
        return loc_data[loc_data.reef_siteid .∈ [selected_locations], :]
    else
        return loc_data[loc_data.UNIQUE_ID .∈ [selected_locations], :]
    end
end

