"""
    pathway_diversity(rs::ResultSet, scens::DataFrame)::DataFrame
    pathway_diversity(rs::ResultSet, scens::DataFrame, option::Symbol)::Vector{Float64}

Compute pathway diversity for all options or only one option. If one option is passed it returns
the vector of probabilities and if no option is passed it returns a dataframe with the probability
and pathway diversity value.
"""
function pathway_diversity(
    rs::ResultSet, scens::DataFrame, idx_scens::Vector{Int64}
)::DataFrame
    options = ADRIA.analysis.option_seed_preference()
    options.probabilities = fill(Float64[], size(options, 1))
    options.pathway_diversity = zeros(size(options, 1))

    for (idx_option, start_option) in enumerate(options.option_name)
        options[idx_option, :probabilities] = pathway_diversity(
            rs, scens, idx_scens, start_option
        )
        options[idx_option, :pathway_diversity] = sum(
            _entropy.(options[idx_option, :probabilities])
        )
    end
    return options[:, [:option_name, :pathway_diversity]]
end
function pathway_diversity(
    rs::ResultSet, scens::DataFrame, idx_scens::Vector{Int64}, option::Symbol
)::Vector{Float64}
    start_time::Int64 = rs.inputs.seed_year_start[1] + rs.inputs.pd_frequency[1]
    end_time::Int64 =
        rs.inputs.seed_year_start[1] + rs.inputs.seed_years[1] - rs.inputs.pd_frequency[1]
    min_locs::Int64 = rs.inputs.min_iv_locations[1]
    mcda_method = ADRIA.mcda_methods()[Int64(rs.inputs.guided[1])]

    idx_option_scens = findall(
        option_ts -> option_ts[Int64(rs.inputs.seed_year_start[1])] == option,
        scens.option_ts
    )
    idx_scens = intersect(idx_scens, idx_option_scens)
    probs = ones(length(idx_scens))

    for (idx_prob, idx_scen) in enumerate(idx_scens)
        option_ts = scens.option_ts[idx_scen]
        for tstep in start_time:Int64(rs.inputs.pd_frequency[1]):end_time
            decision_matrix = rs.decision_matrix_log[timesteps=tstep, scenarios=idx_scen]
            probs[idx_prob] *= ADRIA.analysis.switching_probability(
                option_ts[tstep - 1], decision_matrix, rs.loc_data, mcda_method, min_locs,
                option_ts[tstep]
            )
        end
    end
    return probs
end

"""
    switching_probability(past_option::Symbol, decision_matrix::YAXArray)::Dict{Symbol, Float64}

Compute switching probability for all defined pathway diversity options.

# Arguments
- `past_option` : Symbol with past option.
- `decision_matrix` : Decision matrix with criterias used for MCDA analysis.
- `loc_data` : Location data for all locations.
- `mcda_method` : MCDA method used in selec_locations.
- `min_locs` : Minimum number of locations to be selected by the MCDA algorithm
- `ports` : Dataframe with name and coordinate of ports.
- `weight` : Weigh used to cost_index. The complementary is used to weight the option_similarity.
Defaults to 0.5

# Returns
DataFrame with probability of switching to each option or the probability value for one option if the option is passed.
"""
function switching_probability(
    past_option::Symbol,
    decision_matrix::YAXArray,
    loc_data::DataFrame,
    mcda_method::Union{Function,DataType},
    min_locs::Int64;
    ports::Union{DataFrame,Nothing}=nothing,
    weight::Float64=0.5
)::DataFrame
    options = option_seed_preference()
    options.probability = zeros(size(options, 1))
    valid_locs = collect(1:size(loc_data, 1))

    if ports == nothing
        ports = _ports()
    end

    past_locations = ADRIA.select_locations(
        options[options.option_name .== past_option, :preference][1], decision_matrix,
        mcda_method, valid_locs, min_locs
    )
    past_locations = _filter_locations(loc_data, past_locations)

    for row in eachrow(options)
        if row.option_name == past_option
            option_locations = past_locations
        else
            option_locations = ADRIA.select_locations(
                row.preference, decision_matrix, mcda_method, valid_locs, min_locs
            )
            option_locations = _filter_locations(loc_data, option_locations)
        end
        row.probability += (1 - cost_index(option_locations, ports)) * weight
        row.probability +=
            option_similarity(option_locations, past_locations) * (1 - weight)
    end

    # Normalize probabilities
    sum_probs = sum(options.probability)
    options.probability = options.probability / sum_probs

    return options[:, [:option_name, :probability]]
end
function switching_probability(
    past_option::Symbol,
    decision_matrix::YAXArray,
    loc_data::DataFrame,
    mcda_method::Union{Function,DataType},
    min_locs::Int64,
    option::Symbol;
    ports::Union{DataFrame,Nothing}=nothing,
    weight::Float64=0.5
)::Float64
    result = switching_probability(
        past_option, decision_matrix, loc_data, mcda_method, min_locs; ports, weight
    )
    return first(result[result.option_name .== option, :probability])
end

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
    max_distance::Float64=1500000.0
)::Float64
    if selected_locations1 == selected_locations2
        return 1
    end

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
    max_distance_port::Float64=700000.0,
    max_dispersion::Float64=990000.0
)
    distance_port = _distance_port(locations, ports)
    if distance_port > max_distance_port
        @warn "distance_port bigger than max_distance_port: distance_port: $(distance_port)."
        max_distance_port = distance_port
    end
    distance_port /= max_distance_port

    dispersion = _dispersion(locations)
    if dispersion > max_dispersion
        @warn "Dispersion bigger than max_dispersion: dispersion: $(dispersion)"
        max_dispersion = dispersion
    end
    dispersion /= max_dispersion

    return distance_port * weight + dispersion * (1 - weight)
end

"""
    option_seed_preference(; include_weights::Bool = false)::DataFrame

Return a dataframe with all options, including their names and SeedPreferences object. If include_weights = true it
also returns the weights given to all MCDA criteria.
"""
function option_seed_preference(; include_weights::Bool=false)::DataFrame
    criteria_names = collect(fieldnames(ADRIA.SeedCriteriaWeights))
    column_names = vcat([:option_name], criteria_names)
    option_weights = [
        :heat_stress 1.0 0.0 0.1 0.8 0.2 0.2 0.2 0.1;
        :geographic_spread 0.1 0.0 0.1 0.1 0.2 0.5 0.8 0.1;
        :connectivity 0.1 0.4 1.0 0.1 0.9 0.2 0.2 0.1;
        :functional_diversity 0.1 0.0 0.1 0.1 0.2 0.2 0.2 1.0;
        :balanced 1.0 0.85 0.9 0.95 0.7 0.5 0.5 0.5
    ]
    options = DataFrame(option_weights, column_names)

    scw = ADRIA.SeedCriteriaWeights()
    directions = map(field -> getfield(scw, field).direction, criteria_names)

    options.preference = map(
        row -> ADRIA.SeedPreferences(criteria_names, collect(row[2:end]), directions),
        eachrow(options)
    )

    if include_weights
        return options
    end
    return options[:, [1, end]]
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

"""
    _entropy(x::Float64)::Float64

Return causal entropy for a given value.
"""
function _entropy(x::Float64)::Float64
    if x == 0
        return 0
    end
    return -x * log(x)
end
