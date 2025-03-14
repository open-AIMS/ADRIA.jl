"""
    option_similarity(selected_locations1::Vector{String}, selected_locations2::Vector{String}, loc_data::DataFrame, max_distance::Float64=100000.0)::Float64

Given two DataFrames of locations, it computes a measure of how similar the locations are regarding geographic
distribution. We use the complement of the average pairwise distance between locations for each possible pair of
locations as a similarity measure.
It approximates each location's position as the location's polygon centroid.
The Haversine formula is used to compute the distance.

# Arguments
- `selected_locations1` : DataFrame with the first set of locations.
- `selected_locations2` : DataFrame with the second set of locations.
- `max_distance` : Maximum distance used to normalize the result. Defaults to 60000.0.

# Returns
Normalized mean pairwise distance between locations at selected_locations1 and selected_locations2.
"""
function option_similarity(
    selected_locations1::DataFrame,
    selected_locations2::DataFrame,
    max_distance::Float64=1000000.0
)::Float64
    locations1_coord = ADRIA.centroids(selected_locations1)
    locations2_coord = ADRIA.centroids(selected_locations2)

    pairwise_combination = Iterators.product(locations1_coord, locations2_coord)
    distance = mapreduce(
        coords -> Distances.haversine(coords[1], coords[2]), +, pairwise_combination
    )
    distance /= length(pairwise_combination)
    if distance > max_distance
        throw(ArgumentError("Distance bigger than max_distance: Distance: $(distance)"))
    end

    return 1 - distance / max_distance
end

"""
    _filter_locations(loc_data::DataFrame, selected_locations::Vector{String})
    _filter_locations(loc_data::DataFrame, selected_locations::Vector{Int64})

Filter locations based on the index, reef_siteid or UNIQUE_ID.
"""
function _filter_locations(loc_data::DataFrame, selected_locations::Vector{String})
    if all(isdigit, selected_locations[1])
        return loc_data[loc_data.UNIQUE_ID .∈ [selected_locations], :]
    end
    return loc_data[loc_data.reef_siteid .∈ [selected_locations], :]
end
function _filter_locations(loc_data::DataFrame, selected_locations::Vector{Int64})
    return loc_data[selected_locations, :]
end
