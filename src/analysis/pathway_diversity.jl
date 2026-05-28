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
    idx_scen = idx_scens[1] # index of scenario used to extract model parameters
    start_time::Int64 =
        rs.inputs.seed_year_start[idx_scen] + rs.inputs.pd_frequency[idx_scen]
    end_time::Int64 =
        rs.inputs.seed_year_start[idx_scen] + rs.inputs.seed_years[idx_scen] -
        rs.inputs.pd_frequency[idx_scen]
    min_locs::Int64 = rs.inputs.min_iv_locations[idx_scen]
    mcda_method = ADRIA.mcda_methods()[Int64(rs.inputs.guided[idx_scen])]

    # Find scenarios that start with $option on first seeding step
    idx_option_scens = findall(eachindex(scens.option_ts)) do i
        option_ts = decode_option_ts(
            rs.inputs.option_ts[i], rs.inputs.seed_year_start[i],
            rs.inputs.seed_years[i], rs.inputs.pd_frequency[i],
            size(rs.seed_log, :timesteps)
        )
        option_ts[Int(rs.inputs.seed_year_start[i])] == option
    end
    idx_scens = intersect(idx_scens, idx_option_scens)
    probs = ones(length(idx_scens))

    for (idx_prob, idx_scen) in enumerate(idx_scens)
        option_ts = decode_option_ts(
            rs.inputs.option_ts[idx_scen], rs.inputs.seed_year_start[idx_scen],
            rs.inputs.seed_years[idx_scen],
            rs.inputs.pd_frequency[idx_scen], size(rs.seed_log, :timesteps)
        )
        for tstep in start_time:Int64(rs.inputs.pd_frequency[idx_scen]):end_time
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
- `weights` : Tuple of weights for `(distance_to_port, dispersion, option_similarity)`. Defaults to `(0.3, 0.2, 0.5)`.

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
    weights::NTuple{3,Float64}=(0.3, 0.3, 0.4)
)::DataFrame
    options = option_seed_preference()
    options.probability = zeros(size(options, 1))
    options.selected_locations = [
        Vector{eltype(decision_matrix.location)}() for _ in 1:size(options, 1)
    ]
    valid_locs = collect(1:size(loc_data, 1))

    if ports == nothing
        ports = _ports()
    end

    for row in eachrow(options)
        row.selected_locations = ADRIA.select_locations(
            row.preference, decision_matrix, mcda_method, valid_locs, min_locs
        )
    end

    past_locations = _filter_locations(
        loc_data, options[options.option_name .== past_option, :selected_locations][1]
    )

    for row in eachrow(options)
        option_locations = _filter_locations(loc_data, row.selected_locations)

        common_ids = intersect(option_locations.UNIQUE_ID, past_locations.UNIQUE_ID)
        unique_option_locs = option_locations[option_locations.UNIQUE_ID .∉ [common_ids], :]
        unique_past_locs = past_locations[past_locations.UNIQUE_ID .∉ [common_ids], :]

        row.probability +=
            distance_port_score(unique_option_locs, unique_past_locs, ports) * weights[1]
        row.probability += dispersion_score(option_locations, past_locations) * weights[2]
        row.probability +=
            option_similarity(unique_option_locs, unique_past_locs) * weights[3]
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
    weights::NTuple{3,Float64}=(0.3, 0.2, 0.5)
)::Float64
    result = switching_probability(
        past_option, decision_matrix, loc_data, mcda_method, min_locs; ports, weights
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
    if isempty(selected_locations1) || isempty(selected_locations2)
        return 1.0
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
    option_seed_preference(; include_weights::Bool = false)::DataFrame

Return a dataframe with all options, including their names and SeedPreferences object. If include_weights = true it
also returns the weights given to all MCDA criteria.
"""
function option_seed_preference(; include_weights::Bool=false)::DataFrame
    criteria_names = collect(fieldnames(ADRIA.SeedCriteriaWeights))
    column_names = vcat([:option_name], criteria_names)
    option_weights = [
        :heat_stress 0.99 0.0 0.1 0.8 0.2 0.2 0.2 0.1;
        :geographic_spread 0.1 0.0 0.1 0.1 0.2 0.5 0.8 0.1;
        :connectivity 0.1 0.4 0.99 0.1 0.9 0.2 0.2 0.1;
        :functional_diversity 0.1 0.0 0.1 0.1 0.2 0.2 0.2 0.99;
        :balanced 0.99 0.85 0.9 0.95 0.7 0.5 0.5 0.5
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
    _distance_port(locations::DataFrame, ports::DataFrame)

Compute the mean distance (meters) to the closest port across all locations.
"""
function _distance_port(locations::DataFrame, ports::DataFrame)::Float64
    if isempty(locations)
        return 0.0
    end

    locations_coord = ADRIA.centroids(locations)
    distances = [
        minimum(Distances.haversine.([loc], ports.coordinates))
        for loc in locations_coord
    ]
    return mean(distances)
end

"""
    distance_port_score(option_locs, past_locs, ports; amplify_ranges=true)

Scale-invariant score in [0, 1] comparing how close option locations are to port relative to past
locations. Values above 0.5 mean option is closer to port than past; below 0.5 means farther.
If `amplify_ranges=true` (default), uses `(d²_past − d²_option)/(d²_past + d²_option)` which
amplifies moderate ratios; otherwise uses the basic `(d_past − d_option)/(d_past + d_option)`.
Returns 0.5 when both distances are zero.
"""
function distance_port_score(
    option_locs::DataFrame, past_locs::DataFrame, ports::DataFrame;
    amplify_ranges::Bool=true
)::Float64
    d_option = _distance_port(option_locs, ports)
    d_past = _distance_port(past_locs, ports)
    if d_option + d_past == 0.0
        return 0.5
    end
    if amplify_ranges
        score = (d_past^2 - d_option^2) / (d_past^2 + d_option^2)
    else
        score = (d_past - d_option) / (d_past + d_option)
    end
    return (score + 1) / 2
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

Compute the mean pairwise distance (meters) across locations.
"""
function _dispersion(locations::DataFrame)::Float64
    if isempty(locations)
        return 0.0
    else
        return mean(ADRIA.mean_distance(locations))
    end
end

"""
    dispersion_score(option_locs, past_locs; amplify_ranges=true)

Scale-invariant score in [0, 1] comparing the dispersion of option locations relative to past
locations. Values above 0.5 mean option is less dispersed than past; below 0.5 means more dispersed.
If `amplify_ranges=true` (default), uses `(d²_past − d²_option)/(d²_past + d²_option)` which
amplifies moderate ratios; otherwise uses the basic `(d_past − d_option)/(d_past + d_option)`.
Returns 0.5 when both dispersions are zero.
"""
function dispersion_score(
    option_locs::DataFrame, past_locs::DataFrame; amplify_ranges::Bool=true
)::Float64
    d_option = _dispersion(option_locs)
    d_past = _dispersion(past_locs)
    if amplify_ranges
        score = (d_past^2 - d_option^2) / (d_past^2 + d_option^2)
    else
        score = (d_past - d_option) / (d_past + d_option)
    end
    return (score + 1) / 2
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

"""
    encode_option_ts(combination::Tuple)::Int

Encode a tuple of option symbols as a base-5 integer.
Each position maps to an index 0–4 matching the order of `option_seed_preference()`.
"""
function encode_option_ts(combination::Tuple)::Int
    option_names = ADRIA.analysis.option_seed_preference().option_name
    option_to_idx = Dict(name => i - 1 for (i, name) in enumerate(option_names))
    return sum(option_to_idx[opt] * 5^(i - 1) for (i, opt) in enumerate(combination))
end

"""
    decode_option_ts(encoded, seed_year_start, seed_years, pd_frequency, max_time)

Decode a base-5 integer back into a full option time-series vector of length `max_time`.
Positions outside the seeding window `[seed_year_start, seed_year_start+seed_years)` are `:nothing`.
"""
function decode_option_ts(
    encoded::Real, seed_year_start::Real, seed_years::Real, pd_frequency::Real,
    max_time::Real
)::Vector{Symbol}
    encoded, seed_year_start, seed_years, pd_frequency, max_time =
        Int.((encoded, seed_year_start, seed_years, pd_frequency, max_time))
    option_names = ADRIA.analysis.option_seed_preference().option_name
    number_changes = seed_years ÷ pd_frequency
    decoded_combo = [option_names[(encoded ÷ 5^(i - 1)) % 5 + 1] for i in 1:number_changes]
    ts = fill(:nothing, max_time)
    for t in seed_year_start:(seed_year_start + seed_years - 1)
        ts[t] = decoded_combo[(t - seed_year_start) ÷ pd_frequency + 1]
    end
    return ts
end
