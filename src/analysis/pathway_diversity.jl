# Simplified list of ports in the GBR
# ref: https://elibrary.gbrmpa.gov.au/jspui/retrieve/5391c720-b846-4fae-93c7-0ff53f829ca2/Ports%20and%20Shipping%20Information%20sheet-29May2013.pdf
const PORTS = Dict(
    :quintell_beach => (143.5444588, -12.8437284),
    :cooktown => (145.2475149, -15.4610395),
    :cairns => (145.7808046, -16.9213696),
    :townsville => (146.8332585, -19.2529659),
    :abbot_point => (148.0948982, -19.9220318),
    :hay_point => (149.2737503, -21.2934747),
    :gladstone => (151.2993586, -23.8740713)
)

"""
    option_similarity(selected_locations1::Vector{String}, selected_locations2::Vector{String}, loc_data::DataFrame, max_distance::Float64=100000.0)::Float64

Given two vectors of location ids, it computes a measure of how similar, in terms of geografic distribution, the
locations are. We use as a similarity measure the complement of the average pairwise distance between locations for each
possible pair of locations.
It approximates each location's position as the location's polygon centroid.
The Haversine formula is used to compute the distance.

# Arguments
- `selected_locations1` : Vector with location id (reef_siteid or UNIQUE_ID) of the first set of locations.
- `selected_locations2` : Vector with location id (reef_siteid or UNIQUE_ID) of the second set of locations.
- `loc_data` : Location data for all locations.
- `max_distance` : Maximum distance used to normalize the result. Defaults to 60000.0.

# Returns
Average pairwise distance between locations at selected_locations1 and selected_locations2.
"""
function option_similarity(
    selected_locations1::Vector{String},
    selected_locations2::Vector{String},
    loc_data::DataFrame,
    max_distance::Float64=60000.0
)::Float64
    reefs1_coord = ADRIA.centroids(_filter_locations(loc_data, selected_locations1))
    reefs2_coord = ADRIA.centroids(_filter_locations(loc_data, selected_locations2))

    pairwise_combination = Iterators.product(reefs1_coord, reefs2_coord)
    distance = mapreduce(
        coords -> Distances.haversine(coords[1], coords[2]), +, pairwise_combination
    )
    distance /= length(pairwise_combination)

    return 1 - distance / max_distance
end

"""
    cost_index(selected_locations::Vector{String}, loc_data::DataFrame, ports::Dict{Symbol, Tuple{Float64, Float64}}; weigh::Float64 = 0.6, max_distance_port::Float64 = 200000.0, max_dispersion::Float64 = 10000.0)

Compute a simplified cost index of intervening in selected locations. The index is a value between zero and one and
consider the distance to port for the given locations and the average dispersion of the locations.

# Arguments
- `reefs` : Dataframe with locations that will receive intervention.
- `loc_data` : Location data for all locations.
- `weigh` : Weigh given to distance to port compared with dispersion. If weigh=1 only distance to port is considered and
if weigh=0 only dispersion is considered. Defaults to 0.6.
- `max_distance_port` : Value used to normalize the distance to port measure. Defaults to 200000.0.
- `max_dispersion` : Value used to normalize dispersion of locations. Defaults to 10000.0.

# Returns
Cost index of intervening in selected locations.
"""
function cost_index(
    reefs::DataFrame;
    weigh::Float64=0.6,
    max_distance_port::Float64=200000.0,
    max_dispersion::Float64=10000.0
)
    distance_port = _distance_port(reefs) / max_distance_port
    dispersion = _dispersion(reefs) / max_dispersion
    return distance_port * weigh + dispersion * (1 - weigh)
end
function cost_index(
    selected_locations::Vector{String},
    loc_data::DataFrame;
    weigh::Float64=0.6,
    max_distance_port::Float64=200000.0,
    max_dispersion::Float64=10000.0
)
    reefs = _filter_locations(loc_data, selected_locations)
    return cost_index(reefs; weigh, max_distance_port, max_dispersion)
end

"""
    _filter_locations(loc_data::DataFrame, selected_locations::Vector{String})

Filter locations based on the reef_siteid or UNIQUE_ID.
"""
function _filter_locations(loc_data::DataFrame, selected_locations::Vector{String})
    if all(isdigit, selected_locations[1])
        return loc_data[loc_data.UNIQUE_ID .∈ [selected_locations], :]
    end
    return loc_data[loc_data.reef_siteid .∈ [selected_locations], :]
end

"""
    _distance_port(reefs::DataFrame, ports::Dict{Symbol, Tuple{Float64, Float64}})

Average distance to port for a given geodataframe of reefs.
"""
function _distance_port(
    reefs::DataFrame,
    ports::Dict{Symbol,Tuple{Float64,Float64}}=PORTS
)
    reefs_coord = ADRIA.centroids(reefs)
    distances = [
        minimum(Distances.haversine.([reef_cood], values(ports))) for
        reef_cood in reefs_coord
    ]
    return mean(distances)
end

"""
    _dispersion(reefs::DataFrame)

Mean pairwise distance of reefs in the geodataframe.
"""
function _dispersion(reefs::DataFrame)
    return mean(ADRIA.mean_distance(reefs))
end
