"""
    pathway_diversity(rs::ResultSet, idx_scens::Vector{Int64})::DataFrame
    pathway_diversity(rs::ResultSet, idx_scens::Vector{Int64}, option::Symbol)::Vector{Float64}

Compute pathway diversity for all options or only one option. If one option is passed it returns
the vector of probabilities and if no option is passed it returns a dataframe with the probability
and pathway diversity value.
"""
function pathway_diversity(
    rs::ResultSet, idx_scens::Vector{Int64}, removed_pathways::Int64
)::DataFrame
    max_time = size(rs.seed_log, :timesteps)
    decoded_ts = Dict(
        i => decode_option_ts(
            rs.inputs.option_ts[i], rs.inputs.seed_year_start[i],
            rs.inputs.seed_years[i], rs.inputs.pd_frequency[i], max_time
        )
        for i in idx_scens
    )

    options = copy(PD_OPTIONS())
    options.probabilities = fill(Float64[], size(options, 1))
    options.pathway_diversity = zeros(size(options, 1))

    for (idx_option, start_option) in enumerate(options.option_name)
        probs = pathway_diversity(rs, idx_scens, decoded_ts, start_option)
        options[idx_option, :probabilities] = probs
        if removed_pathways > 0
            n_total = length(probs) + removed_pathways
            resampled = ADRIA.StatsBase.sample(probs, n_total; replace=true)
            resampled = resampled ./ sum(resampled)
            options[idx_option, :pathway_diversity] = sum(_entropy.(resampled))
        else
            probs = probs ./ sum(probs)
            options[idx_option, :pathway_diversity] = sum(_entropy.(probs))
        end
    end
    return options[:, [:option_name, :pathway_diversity]]
end
function pathway_diversity(
    rs::ResultSet, idx_scens::Vector{Int64}, option::Symbol
)::Vector{Float64}
    max_time = size(rs.seed_log, :timesteps)
    decoded_ts = Dict(
        i => decode_option_ts(
            rs.inputs.option_ts[i], rs.inputs.seed_year_start[i],
            rs.inputs.seed_years[i], rs.inputs.pd_frequency[i], max_time
        )
        for i in idx_scens
    )
    return pathway_diversity(rs, idx_scens, decoded_ts, option)
end
function pathway_diversity(
    rs::ResultSet, idx_scens::Vector{Int64},
    decoded_ts::Dict{Int64,Vector{Symbol}}, option::Symbol
)::Vector{Float64}
    idx_scen = idx_scens[1] # index of scenario used to extract model parameters
    start_time::Int64 =
        rs.inputs.seed_year_start[idx_scen] + rs.inputs.pd_frequency[idx_scen]
    end_time::Int64 =
        rs.inputs.seed_year_start[idx_scen] + rs.inputs.seed_years[idx_scen] -
        rs.inputs.pd_frequency[idx_scen]
    min_locs::Int64 = rs.inputs.min_iv_locations[idx_scen]
    mcda_method = ADRIA.mcda_methods()[Int64(rs.inputs.guided[idx_scen])]

    idx_option_scens = filter(
        i -> decoded_ts[i][Int(rs.inputs.seed_year_start[i])] == option, idx_scens
    )
    probs = ones(length(idx_option_scens))

    # Pre-fetch all decision matrices for selected scenarios and timesteps to avoid repeated I/O
    pd_freq = Int64(rs.inputs.pd_frequency[idx_scen])
    tsteps = collect(start_time:pd_freq:end_time)
    cached_decision_matrices = readcubedata(
        rs.decision_matrix_log[timesteps=tsteps, scenarios=idx_option_scens]
    )

    # Cache outcome metrics and pre-compute per-tstep option performance (read-only in threads)
    rel_cover_data = Array(rs.outcomes[:relative_cover][scenarios=idx_scens])
    fd_data = Array(rs.outcomes[:coral_evenness][scenarios=idx_scens])
    max_time = size(rel_cover_data, 1)
    scen_to_idx = Dict(s => i for (i, s) in enumerate(idx_scens))

    option_perf_by_tstep = Dict(
        t => _compute_option_perf(
            rel_cover_data, fd_data, idx_scens, scen_to_idx,
            decoded_ts, t, min(t + pd_freq - 1, max_time)
        )
        for t in tsteps
    )

    Threads.@threads for idx_prob in eachindex(idx_option_scens)
        idx_scen = idx_option_scens[idx_prob]
        option_ts = decoded_ts[idx_scen]
        for tstep in tsteps
            decision_matrix = cached_decision_matrices[
                timesteps=At(tstep), scenarios=At(idx_scen)
            ]
            probs[idx_prob] *= ADRIA.analysis.switching_probability(
                option_ts[tstep - 1], decision_matrix, rs.loc_data, mcda_method, min_locs,
                option_ts[tstep], option_perf_by_tstep[tstep]
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
- `weights` : 5-tuple of weights for `(cum_rel_tac_diff, cum_fd_diff, distance_to_port, dispersion, option_similarity)`. Defaults to `(0.35, 0.25, 0.1, 0.1, 0.2)`.
- `option_perf` : Pre-computed per-option performance metrics for the relevant timestep (from `_compute_option_perf`), keyed by option symbol.

# Returns
DataFrame with probability of switching to each option or the probability value for one option if the option is passed.
"""
function switching_probability(
    past_option::Symbol,
    decision_matrix::YAXArray,
    loc_data::DataFrame,
    mcda_method::Union{Function,DataType},
    min_locs::Int64,
    option_perf::Dict{Symbol,YAXArray};
    ports::Union{DataFrame,Nothing}=nothing,
    weights::NTuple{5,Float64}=(0.35, 0.25, 0.1, 0.1, 0.2),
)::DataFrame
    options = copy(PD_OPTIONS())
    options.probability = zeros(size(options, 1))
    options.selected_locations = [
        Vector{eltype(decision_matrix.location)}() for _ in 1:size(options, 1)
    ]
    valid_locs = collect(1:size(loc_data, 1))

    if ports === nothing
        ports = PD_PORTS
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

        # Outcome-based metrics: cum_rel_tac_diff, cum_fd_diff (weights[1..2])
        # Compares the per-location performance of a scenario holding past_option against
        # a (possibly different) scenario holding the candidate option, paired by location.
        if haskey(option_perf, past_option) && haskey(option_perf, row.option_name)
            past_perf = option_perf[past_option]
            cand_perf = option_perf[row.option_name]

            δ_rel_tac = cand_perf[metric=At(:rel_tac)] ./ past_perf[metric=At(:rel_tac)] .- 1
            δ_fd = cand_perf[metric=At(:fd)] ./ past_perf[metric=At(:fd)] .- 1

            row.probability += weights[1] * two_sided_cvar(δ_rel_tac; σ=_σ_rel_tac)
            row.probability += weights[2] * two_sided_cvar(δ_fd; σ=_σ_fd)
        end

        # Spatial metrics: distance_to_port, dispersion, option_similarity (weights[3..5])
        row.probability +=
            distance_port_score(unique_option_locs, unique_past_locs, ports) * weights[3]
        row.probability += dispersion_score(option_locations, past_locations) * weights[4]
        row.probability +=
            option_similarity(unique_option_locs, unique_past_locs) * weights[5]
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
    option::Symbol,
    option_perf::Dict{Symbol,YAXArray};
    ports::Union{DataFrame,Nothing}=nothing,
    weights::NTuple{5,Float64}=(0.35, 0.25, 0.1, 0.1, 0.2),
)::Float64
    result = switching_probability(
        past_option, decision_matrix, loc_data, mcda_method, min_locs, option_perf;
        ports, weights
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
const _PD_OPTIONS = Ref{Union{Nothing,DataFrame}}(nothing)
function PD_OPTIONS()
    if isnothing(_PD_OPTIONS[])
        _PD_OPTIONS[] = option_seed_preference()
    end
    return _PD_OPTIONS[]
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
const PD_PORTS = _ports()

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
    option_names = PD_OPTIONS().option_name
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
    option_names = PD_OPTIONS().option_name
    number_changes = seed_years ÷ pd_frequency
    decoded_combo = [option_names[(encoded ÷ 5^(i - 1)) % 5 + 1] for i in 1:number_changes]
    ts = fill(:nothing, max_time)
    for t in seed_year_start:(seed_year_start + seed_years - 1)
        ts[t] = decoded_combo[(t - seed_year_start) ÷ pd_frequency + 1]
    end
    return ts
end

# Sigma sensitivity parameters for two_sided_cvar, tuned per metric's relative-change scale
const _σ_rel_tac = 0.004   # cum_rel_tac_diff
const _σ_fd = 0.003        # cum_fd_diff

"""
    two_sided_cvar(delta; tail_fraction, σ)

Map the two-sided CVaR of `delta` to [0, 1] via tanh. Returns > 0.5 when the upper
tail dominates (net beneficial), < 0.5 when the lower tail dominates (net detrimental),
and 0.5 when `delta` is empty.

`delta` is a per-location relative change vector (`candidate ./ past .- 1`).
`tail_fraction` is the fraction of locations in each tail (default 3%).
"""
function two_sided_cvar(
    delta::AbstractVector;
    tail_fraction::Float64=0.03,
    σ::Float64=0.001,
)::Float64
    isempty(delta) && return 0.5
    N = length(delta)
    k = max(1, min(floor(Int, tail_fraction * N), N ÷ 2))
    s = sort(delta)
    cvar_lower = mean(s[1:k])
    cvar_upper = mean(s[(end - k + 1):end])
    raw = cvar_upper - abs(cvar_lower)
    return (tanh(raw / σ) + 1) / 2
end

"""
    _compute_option_perf(rel_cover_data, fd_data, all_scens, scen_to_idx,
                         decoded_ts, tstep, t_end)

For each option, find one representative scenario that holds that option at `tstep`
(the first match found in `all_scens`) and compute its cumulative per-location
performance over the window [tstep, t_end]. Both the past (counterfactual) and
candidate option representatives are looked up at the same timestep `tstep`, so that
performance is compared over an identical future window.

Returns a `Dict{Symbol, YAXArray}` mapping option symbol to its representative
scenario's YAXArray with dims `(metric=[:rel_tac, :fd], location=1:n_locs)`.
"""
function _compute_option_perf(
    rel_cover_data::AbstractArray{<:Real,3},
    fd_data::AbstractArray{<:Real,3},
    all_scens::AbstractVector{Int64},
    scen_to_idx::Dict{Int64,Int64},
    decoded_ts::Dict{Int64,Vector{Symbol}},
    tstep::Int64,
    t_end::Int64
)::Dict{Symbol,YAXArray}
    n_locs = size(rel_cover_data, 2)
    representative_scen = Dict{Symbol,Int64}()
    for scen in all_scens
        opt = decoded_ts[scen][tstep]
        opt === :nothing && continue
        get!(representative_scen, opt, scen)
    end

    result = Dict{Symbol,YAXArray}()

    for (opt, scen) in representative_scen
        scen_idx = scen_to_idx[scen]
        rel_cover_window = rel_cover_data[tstep:t_end, :, scen_idx]   # (n_ts, n_locs)
        evenness_window = fd_data[tstep:t_end, :, scen_idx]           # (n_ts, n_locs)

        # cum_rel_tac: sum over time, per location
        rel_tac_vals = vec(sum(rel_cover_window; dims=1))
        # cum_fd: sum over time, per location
        fd_vals = vec(sum(evenness_window; dims=1))

        result[opt] = DataCube(
            permutedims(hcat(rel_tac_vals, fd_vals));
            metric=[:rel_tac, :fd], location=1:n_locs
        )
    end

    return result
end
