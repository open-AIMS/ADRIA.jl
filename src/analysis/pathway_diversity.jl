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
    cost_index(locations::DataFrame, ports::DataFrame; weight::Float64=0.6, max_distance_port::Float64=250000.0, max_dispersion::Float64=800000.0)
Compute a simplified cost index of intervening in selected locations. The cost index is a value between zero and one,
considering the distance to the closest port for each given location and the average dispersion of the locations.

# Arguments
- `locations` : Dataframe with locations that will receive intervention.
- `ports` : Dataframe with name and coordinate of ports.
- `weight` : Weight for the distance to port compared with dispersion. If weight=1 only distance to port is considered,
and if weight=0 only dispersion is considered. Defaults to 0.6.
- `max_distance_port` : Value used to normalize the distance to port measure. Defaults to 200000.0.
- `max_dispersion` : Value used to normalize dispersion of locations. Defaults to 10000.0.

# Returns
Cost index of intervening in selected locations.
"""
function cost_index(
    locations::DataFrame,
    ports::DataFrame;
    weight::Float64=0.6,
    max_distance_port::Float64=250000.0,
    max_dispersion::Float64=800000.0
)
    distance_port = _distance_port(locations, ports)
    if distance_port > max_distance_port
        throw(
            ArgumentError(
                "distance_port bigger than max_distance_port: distance_port: $(distance_port)"
            )
        )
    end
    distance_port /= max_distance_port

    dispersion = _dispersion(locations)
    if dispersion > max_dispersion
        throw(
            ArgumentError(
                "Dispersion bigger than max_dispersion: dispersion: $(dispersion)"
            )
        )
    end
    dispersion /= max_dispersion

    return distance_port * weight + dispersion * (1 - weight)
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

"""
    _distance_port(locations::DataFrame, ports::Dict{Symbol, Tuple{Float64, Float64}})

Compute the mean distance to the closest port across all locations in a given GeoDataFrame.
"""
function _distance_port(locations::DataFrame, ports::DataFrame)
    locations_coord = ADRIA.centroids(locations)
    distances = [
        minimum(Distances.haversine.([location_coord], ports.coordinates))
        for location_coord in locations_coord
    ]
    return mean(distances)
end

"""
    _ports()

Return dataframe with simplified list ports in the GBR and their coordinates.
Once we get a definitive list this will be moved to the datapackage spatial data.

# References
https://elibrary.gbrmpa.gov.au/jspui/retrieve/5391c720-b846-4fae-93c7-0ff53f829ca2/Ports%20and%20Shipping%20Information%20sheet-29May2013.pdf

"""
function _ports()
    return DataFrame([
        (name=:yeppoon, coordinates=(150.7864023, -23.1602957)),
        (name=:airlie_beach, coordinates=(148.7112726, -20.2661893)),
        (name=:cairns, coordinates=(145.76625, -16.92304)),
        (name=:townsville, coordinates=(146.80569, -19.26639)),
        (name=:gladstone, coordinates=(151.25775, -23.84852))
    ])
end

"""
    _dispersion(locations::DataFrame)

Compute the mean pairwise distance across locations in a given GeoDataFrame.
"""
function _dispersion(locations::DataFrame)
    return mean(ADRIA.mean_distance(locations))
end
