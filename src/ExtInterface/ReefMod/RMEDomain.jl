using Distributions
using Statistics

using CSV, DataFrames, ModelParameters

using Glob

using ADRIA: SimConstants, Domain, GDF
using ADRIA.decision:
    DecisionThresholds,
    DecisionWeights,
    DepthThresholds

using ADRIA.decision:
    SeedCriteriaWeights,
    FogCriteriaWeights

abstract type AbstractReefModDomain <: Domain end

mutable struct RMEDomain <: AbstractReefModDomain
    const name::String
    RCP::String
    env_layer_md::ADRIA.EnvLayer
    scenario_invoke_time::String  # time latest set of scenarios were run
    const conn::YAXArray{Float64}
    const loc_data::DataFrames.DataFrame
    const loc_id_col::String
    const cluster_id_col::String
    init_coral_cover::YAXArray{Float64}
    const coral_growth::CoralGrowth
    const loc_ids::Vector{String}
    dhw_scens::YAXArray{Float64}

    # `wave_scens` holds empty wave data to maintain compatibility with
    # ADRIA's dMCDA methods
    wave_scens::YAXArray{Float64}

    # `cyclone_mortality_scens` holds dummy cyclone mortality data to maintain compatibility
    # with ADRIA's dMCDA methods
    cyclone_mortality_scens::YAXArray{Float64}

    # Strategy target locations
    seed_target_locations::Vector{String}       # locations eligible for seeding
    fog_target_locations::Vector{String}        # locations eligible for fogging
    shade_target_locations::Vector{String}    # locations eligible for shading

    model::ModelParameters.Model
    sim_constants::SimConstants
end

"""
    _standardize_cluster_ids!(spatial_data::DataFrame)::Nothing

Standardize cluster id column name
["LOC_NAME_S" or "reef_name] to ["cluster_id"]
"""
function _standardize_cluster_ids!(spatial_data::DataFrame)::Nothing
    try
        rename!(spatial_data, Dict("LOC_NAME_S" => "cluster_id"))
    catch err
        if !(err isa ArgumentError)
            rethrow(err)
        end

        rename!(spatial_data, Dict("reef_name" => "cluster_id"))
    end

    return nothing
end

"""
    _manual_id_corrections!(spatial_data::DataFrame)::Nothing

Manual corrections to ensure correct alignment of reefs by their order.
"""
function _manual_id_corrections!(spatial_data::DataFrame, id_list::DataFrame)::Nothing
    try
        # Re-order spatial data to match RME dataset
        # MANUAL CORRECTION
        spatial_data[spatial_data.LABEL_ID .== "20198", :LABEL_ID] .= "20-198"
        id_order = [
            first(findall(x .== spatial_data.LABEL_ID)) for x in string.(id_list[:, 1])
        ]
        spatial_data[!, :] = spatial_data[id_order, :]

        # Check that the two lists of location ids are identical
        @assert isempty(findall(spatial_data.LABEL_ID .!= id_list[:, 1]))
    catch err
        if !(err isa ArgumentError)
            # Updated dataset already has corrections/modifications included.
            # If it is not an ArgumentError, something else has happened.
            rethrow(err)
        end

        if contains(string(err), ":LABEL_ID")
            # Hacky manual catch for canonical dataset.
            # Geopackage does not contain LABEL_ID, so check RME_GBRMPA_ID column.
            # The original file had "20198", but later revisions of the RME reef id list
            # seems to correct it.
            spatial_data[spatial_data.RME_GBRMPA_ID .== "20198", :RME_GBRMPA_ID] .= "20-198"
        end
    end

    return nothing
end

"""
    load_domain(::Type{RMEDomain}, fn_path::String, RCP::String; timeframe::Tuple{Int64, Int64}=(2022, 2100))::RMEDomain

Load a ReefMod Engine dataset.

# Arguments
- `RMEDomain` : DataType
- `fn_path` : path to ReefMod Engine dataset
- `RCP` : Representative Concentration Pathway scenario ID
- `timeframe` : Timeframe for simulations to be run. Defaults to (2022, 2100)
- `force_single_reef` : If true use a single reef dataset with user defined total and
habitable areas.
- `single_reef_total_area` : Reef total area in m² to be used when `force_single_reef` is `true`.
- `single_reef_k` : Reef k ∈ [0.0, 1.0] ro be used when `force_single_reef` is `true`.

# Returns
RMEDomain
"""
function load_domain(
    ::Type{RMEDomain},
    fn_path::String,
    RCP::String;
    timeframe::Union{Nothing,Tuple{Int64,Int64}}=nothing,
    force_single_reef::Bool=false,
    force_single_reef_id::String="",
    single_reef_total_area::Float64=1_000_000.0,
    single_reef_k::Float64=0.5
)::RMEDomain
    isdir(fn_path) ? true : error("Path does not exist or is not a directory.")

    data_files = joinpath(fn_path, "data_files")

    # Load the geopackage file
    gpkg_path = joinpath(data_files, "region", "reefmod_gbr.gpkg")
    spatial_data = GDF.read(gpkg_path)

    # Adjust spatial data if `force_single_reef == true`
    single_reef_idx = collect(1:nrow(spatial_data))  # select all locations by default
    if force_single_reef
        single_reef_idx = if isempty(force_single_reef_id)
            [1]
        else
            findall(spatial_data.UNIQUE_ID .== force_single_reef_id)
        end

        spatial_data = spatial_data[single_reef_idx, :]
    end

    # Find initial coral cover start year
    initial_csv_files = readdir(joinpath(data_files, "initial_csv"))
    icc_filename = initial_csv_files[occursin.("coral_sp", initial_csv_files)][1]
    icc_start_year::Int64 = parse(Int64, split(icc_filename, r"[_.]")[end - 1])

    # Load dhw_scens
    dhw_scens = load_DHW(RMEDomain, data_files, RCP; timeframe=timeframe)

    if timeframe == nothing
        timeframe = (icc_start_year, dhw_scens.timesteps[end])
        dhw_scens = dhw_scens[timesteps=At(timeframe[1]:timeframe[2])]
    end

    force_single_reef && (dhw_scens = dhw_scens[:, [1], :])
    loc_ids::Vector{String} = collect(dhw_scens.locs)

    # Standardize IDs to use for site/reef and cluster
    _standardize_cluster_ids!(spatial_data)

    reef_id_col = if "RME_UNIQUE_ID" ∈ names(spatial_data)
        "RME_UNIQUE_ID"
    elseif "UNIQUE_ID" ∈ names(spatial_data)
        "UNIQUE_ID"
    else
        error("No reef id column identified")
    end
    cluster_id_col = reef_id_col

    reef_ids::Vector{String} = spatial_data[:, reef_id_col]

    # Load accompanying ID list
    # TODO: Create canonical geopackage file that aligns all IDs.
    #       Doing this removes the need for the manual correction below and removes the
    #       dependency on this file.

    id_file_path = glob(glob"*.csv", joinpath(data_files, "id"))[1]
    id_list = CSV.read(
        id_file_path,
        DataFrame;
        header=false,
        comment="#"
    )

    if force_single_reef
        id_list = id_list[spatial_data.RME_GBRMPA_ID .== id_list[:, 1], :]
    end

    _manual_id_corrections!(spatial_data, id_list)

    # Check that the two lists of location ids are identical
    try
        @assert isempty(findall(spatial_data.LABEL_ID .!= id_list[:, 1]))
    catch
        @assert isempty(findall(spatial_data.RME_GBRMPA_ID .!= id_list[:, 1]))
    end

    # Overwrite spatial area when force_single_reef is true
    spatial_data[:, :area] .= if force_single_reef
        [single_reef_total_area]
    else
        # Convert area in km² to m²
        id_list[:, 2] .* 1e6
    end

    # Calculate `k` area (1.0 - "ungrazable" area)
    spatial_data[:, :k] .= if force_single_reef
        [single_reef_k]
    else
        1.0 .- id_list[:, 3]
    end

    dist_matrix = distance_matrix(spatial_data)
    spatial_data.mean_to_neighbor .= nearest_neighbor_distances(dist_matrix, 10)

    # Need to load initial coral cover after we know `k` area.
    init_coral_cover::YAXArray{Float64} = load_initial_cover(
        RMEDomain, data_files, loc_ids, spatial_data
    )

    conn_data::YAXArray{Float64} = load_connectivity_csv(
        RMEDomain, data_files, loc_ids; force_single_reef=single_reef_idx
    )

    # Set all site depths to 7m below sea level
    # (ReefMod does not account for depth)
    # Ensure column is of float type
    spatial_data[:, :depth_med] .= 7.0
    spatial_data[!, :depth_med] = convert.(Float64, spatial_data[!, :depth_med])

    # Add GBRMPA zone type info as well
    gbr_zt_path = joinpath(data_files, "region", "gbrmpa_zone_type.csv")
    gbr_zone_types = CSV.read(gbr_zt_path, DataFrame; types=String)
    missing_rows = ismissing.(gbr_zone_types[:, "GBRMPA Zone Types"])
    gbr_zone_types[missing_rows, "GBRMPA Zone Types"] .= ""
    zones = gbr_zone_types[:, "GBRMPA Zone Types"]
    zones = replace.(zones, "Zone" => "", " " => "")
    spatial_data[:, :zone_type] .= force_single_reef ? zones[single_reef_idx] : zones

    # This loads cyclone categories, not mortalities, so ignoring for now.
    # cyc_scens = load_cyclones(RMEDomain, data_files, loc_ids, timeframe)

    timeframe_range = timeframe[1]:timeframe[2]

    wave_scens::YAXArray{Float64} = ZeroDataCube(;
        T=Float64, timesteps=timeframe_range, locs=loc_ids, scenarios=[1]
    )

    functional_groups = functional_group_names()
    cyc_scens::YAXArray{Float64} = ZeroDataCube(;
        T=Float64,
        timesteps=timeframe_range,
        locs=loc_ids,
        species=functional_groups,
        scenarios=[1]
    )

    env_md = EnvLayer(
        fn_path,
        gpkg_path,
        reef_id_col,
        cluster_id_col,
        "",
        "",
        "",
        "",
        timeframe_range
    )

    criteria_weights::Vector{Union{DecisionWeights,DecisionThresholds}} = [
        SeedCriteriaWeights(),
        FogCriteriaWeights(),
        DepthThresholds()
    ]

    model::Model = Model((
        EnvironmentalLayer(dhw_scens, wave_scens, cyc_scens),
        Intervention(),
        criteria_weights...,
        Coral(),
        GrowthAcceleration()
    ))

    return adjust_sampling_bounds(
        RMEDomain(
            "ReefMod Engine",
            RCP,
            env_md,
            "",
            conn_data,
            spatial_data,
            reef_id_col,
            cluster_id_col,
            init_coral_cover,
            CoralGrowth(nrow(spatial_data)),
            reef_ids,
            dhw_scens,
            wave_scens,
            cyc_scens,
            reef_ids,
            reef_ids,
            reef_ids,
            model,
            SimConstants()
        )
    )
end

function adjust_sampling_bounds(dom::Domain)::Domain
    _timesteps = timesteps(dom)
    timeframe = _timesteps[end] - _timesteps[1]

    _year_start_factors = year_start_factors(dom)
    to_update_mask = getindex.(_year_start_factors[:, :dist_params], 2) .> timeframe
    to_update = _year_start_factors[to_update_mask, :]
    to_update_fieldnames = Tuple(to_update[:, :fieldname])

    @info "Adjusting sampling bounds for factors $to_update_fieldnames to fit within domain timeframe."
    new_dist_params = NamedTuple{to_update_fieldnames}(
        (d[1], timeframe) for d in to_update.dist_params
    )
    set_factor_bounds!(dom; new_dist_params...)

    return dom
end

"""
    _get_relevant_files(fn_path, ident)

Filter files found in given path down to those that have the provided identifier.

# Arguments
- `fn_path` : directory of files
- `ident` : keep files with `ident` in their filename.
"""
function _get_relevant_files(fn_path::String, ident::String)
    valid_files = filter(isfile, readdir(fn_path; join=true))
    return filter(x -> occursin(ident, x), valid_files)
end

function _data_folder_path(data_path, target_folder_name; file_extension="csv")
    data_folder_names = readdir(data_path)
    target_folder_mask = [
        occursin(Regex(target_folder_name), file) for file in data_folder_names
    ]
    target_folder_names = data_folder_names[target_folder_mask]
    target_folder_name = if length(target_folder_names) == 1
        target_folder_names[1]
    else
        dhw_csv_mask = [
            occursin(Regex(file_extension), csv_file) for csv_file in target_folder_names
        ]
        target_folder_names[dhw_csv_mask][1]
    end

    return joinpath(data_path, target_folder_name)
end

"""
    load_DHW(::Type{RMEDomain}, data_path::String, rcp::String, timeframe=(2022, 2100))::YAXArray

Loads ReefMod DHW data as a datacube.

# Arguments
- `RMEDomain`
- `data_path` : path to ReefMod data
- `rcp` : RCP identifier
- `timeframe` : range of years to represent.

# Returns
YAXArray[timesteps, locs, scenarios]
"""
function load_DHW(
    ::Type{RMEDomain}, data_path::String, rcp::String;
    timeframe::Union{Nothing,Tuple{Int,Int}}=nothing
)::YAXArray
    dhw_path = _data_folder_path(data_path, "dhw")
    rcp_files = _get_relevant_files(dhw_path, rcp)
    rcp_files = filter(x -> occursin("SSP", x), rcp_files)
    if isempty(rcp_files)
        ArgumentError("No DHW data files found in: $(dhw_path)")
    end

    first_file = CSV.read(rcp_files[1], DataFrame)
    loc_ids = String.(first_file[:, 1])

    data_tf = parse.(Int64, names(first_file[:, 2:end]))

    _timeframe = isnothing(timeframe) ? extrema(data_tf) : timeframe

    tf_start = findall(_timeframe[1] .∈ data_tf)[1]
    tf_end = findall(_timeframe[2] .∈ data_tf)[1]

    d1 = first_file[:, (tf_start + 1):(tf_end + 1)]
    data_shape = reverse(size(d1))
    data_cube = zeros(data_shape..., length(rcp_files))
    data_cube[:, :, 1] .= Matrix(d1)'

    local tf_start
    local tf_end
    keep_ds = fill(true, length(rcp_files))
    for (i, rcp_data_fn) in enumerate(rcp_files[2:end])
        d = CSV.read(rcp_data_fn, DataFrame)[:, 2:end]

        if size(d, 1) == 0
            @info "Empty file?" rcp_data_fn
            continue
        end

        data_tf = parse.(Int64, names(d))
        try
            tf_start = findfirst(_timeframe[1] .∈ data_tf)
            tf_end = findfirst(_timeframe[2] .∈ data_tf)
        catch err
            if !(err isa BoundsError)
                rethrow(err)
            end

            @warn "Building DHW: Could not find matching time frame, skipping $rcp_data_fn"
            keep_ds[i] = false  # mark scenario for removal
            continue
        end

        data_cube[:, :, i + 1] .= Matrix(d[:, tf_start:tf_end])'
    end

    # Only return valid scenarios
    return DataCube(
        data_cube[:, :, keep_ds];
        timesteps=_timeframe[1]:_timeframe[2],
        locs=loc_ids,
        scenarios=rcp_files[keep_ds]
    )
end

"""
    load_connectivity(::Type{RMEDomain}, data_path::String, loc_ids::Vector{String})::YAXArray

Loads the average connectivity from  binary connectivity files and returns them as a
YAXArray. To load from the `csv` connectivity files, refer to`load_connectivity_csv`.

# Arguments
- `RMEDomain`
- `data_path` : path to ReefMod data
- `loc_ids` : location ids

# Returns
YAXArray with dimensions (Source ⋅ Sink) and size (`n_locations` ⋅ `n_locations`).
"""
function load_connectivity(
    ::Type{RMEDomain}, data_path::String, loc_ids::Vector{String}
)::YAXArray
    conn_path = joinpath(data_path, "con_bin")
    conn_files = _get_relevant_files(conn_path, "CONNECT_ACRO")
    if isempty(conn_files)
        ArgumentError("No CONNECT_ACRO data files found in: $(conn_path)")
    end

    n_locs = length(loc_ids)
    tmp_mat = zeros(n_locs, n_locs, length(conn_files))
    for (i, fn) in enumerate(conn_files)
        # File pattern used is "CONNECT_ACRO_[YEAR]_[DAY].bin"
        # We use a clunky regex approach to identify the year.
        # tmp = replace(split(fn, r"(?=CONNECT_ACRO_[0-9,4]+)")[2], "CONNECT_ACRO_"=>"")
        # year_id = split(tmp, "_")[1]

        # Turns out, there's only data for each year so just read in directly
        # Have to read in binary data - read first two values as Int32, and the rest
        # as Float32. Then reshape into a square (n_locs * n_locs) matrix.
        data = IOBuffer(read(fn))
        x = read(data, Int32)
        y = read(data, Int32)

        ds = Vector{Float32}(undef, x * y)
        tmp_mat[:, :, i] .= reshape(read!(data, ds), (n_locs, n_locs))
    end

    # Mean over all years
    conn_data::Matrix{Float64} = dropdims(mean(tmp_mat; dims=3); dims=3)
    return DataCube(
        conn_data;
        Source=loc_ids,
        Sink=loc_ids
    )
end

"""
    load_connectivity_csv(::Type{RMEDomain}, data_path::String, loc_ids::Vector{String})::YAXArray

Loads the average connectivity from connectivity `csv` files. To load the binary version
of the connectivity files refer to `load_connectivity`.

# Arguments
- `RMEDomain`
- `data_path` : path to ReefMod data
- `loc_ids` : location ids

# Returns
YAXArray with dimensions (Source ⋅ Sink) and size (`n_locations` ⋅ `n_locations`).
"""
function load_connectivity_csv(
    ::Type{RMEDomain}, data_path::String, loc_ids::Vector{String};
    force_single_reef::Union{Nothing,Vector{Int64}}=nothing
)::YAXArray
    conn_path = _data_folder_path(data_path, "con")
    conn_files = _get_relevant_files(conn_path, "CONNECT_ACRO")
    if isempty(conn_files)
        ArgumentError("No CONNECT_ACRO data files found in: $(conn_path)")
    end

    # Ensure compatibility when `force_single_reef == true`
    locs_selector =
        isnothing(force_single_reef) ? (:, :) : (force_single_reef, force_single_reef)

    n_locs = length(loc_ids)
    tmp_mat = zeros(n_locs, n_locs, length(conn_files))
    for (i, fpath) in enumerate(conn_files)
        # TODO Allow some kind of customization when `force_single_reef == true`
        tmp_mat[:, :, i] .= Matrix(CSV.read(fpath, DataFrame; header=false))[locs_selector...]
    end

    # Mean over all years
    conn_data::Matrix{Float64} = dropdims(mean(tmp_mat; dims=3); dims=3)
    return DataCube(
        conn_data;
        Source=loc_ids,
        Sink=loc_ids
    )
end

"""
    load_cyclones(::Type{RMEDomain}, data_path::String, loc_ids::Vector{String}, tf::Tuple{Int64, Int64})::YAXArray

# Arguments
- `RMEDomain`
- `data_path` : path to ReefMod data
- `loc_ids` : location ids
- `tf` : the time frame represented by the dataset.
    Subsets the data in cases where `tf[2] < length(data)` such that only
    `1:(tf[2] - tf[1])+1` is read in.

# Returns
YAXArray[timesteps, locs, scenarios]
"""
function load_cyclones(
    ::Type{RMEDomain},
    data_path::String,
    loc_ids::Vector{String},
    tf::Tuple{Int64,Int64}
)::YAXArray
    # NOTE: This reads from the provided CSV files
    #       Replace with approach that reads directly from binary files
    #       Currently cannot get values in binary files to match
    cyc_path = joinpath(data_path, "cyc_csv")
    cyc_files = _get_relevant_files(cyc_path, "Cyclone_simulation_")
    if length(cyc_path) == 0
        ArgumentError("No cyclone data files found in: $(cyc_path)")
    end

    num_years = 100
    cyc_data = zeros(length(loc_ids), num_years, length(cyc_files))
    for (i, fn) in enumerate(cyc_files)
        # Read each cyclone trajectory
        cyc_data[:, :, i] = Matrix(
            CSV.read(fn, DataFrame; drop=[1], header=false, comment="#")
        )
    end

    # Cut down to the given time frame assuming the first entry represents the first index
    cyc_data = permutedims(cyc_data, (2, 1, 3))[1:((tf[2] - tf[1]) + 1), :, :]
    return DataCube(
        cyc_data;
        timesteps=tf[1]:tf[2],
        locs=loc_ids,
        scenarios=1:length(cyc_files)
    )
end

"""
    load_initial_cover(::Type{RMEDomain}, data_path::String, loc_ids::Vector{String}, loc_data::DataFrame)::YAXArray

# Arguments
- `RMEDomain`
- `data_path` : path to ReefMod data
- `loc_ids` : location ids

# Returns
YAXArray[locs, species]
"""
function load_initial_cover(
    ::Type{RMEDomain}, data_path::String, loc_ids::Vector{String}, loc_data::DataFrame
)::YAXArray
    icc_path = _data_folder_path(data_path, "initial")
    icc_files = _get_relevant_files(icc_path, "coral_")

    # ADRIA no longer models Arborescent Acropora
    file_mask = [!contains(filename, "coral_sp1") for filename in icc_files]
    icc_files = icc_files[file_mask]

    if isempty(icc_files)
        ArgumentError("No coral cover data files found in: $(icc_path)")
    end

    # Shape is locations, scenarios, species
    icc_data = zeros(length(loc_ids), 20, length(icc_files))

    for (i, fn) in enumerate(icc_files)
        # This ensures that this works when `force_single_reef = true`
        icc_loc_ids = dropdims(
            Matrix(
                CSV.read(fn, DataFrame; select=[1], header=false, comment="#")
            ); dims=2)

        icc_data[:, :, i] = Matrix(
            CSV.read(fn, DataFrame; drop=[1], header=false, comment="#")
        )[
            icc_loc_ids .== loc_ids, :
        ]
    end

    # Use ReefMod distribution for coral size class population (shape parameters have units log(cm^2))
    # as suggested by YM (pers comm. 2023-08-08 12:55pm AEST). Distribution is used to split ReefMod initial
    # species covers into ADRIA's size classes by weighting with the cdf.
    reef_mod_area_dist = LogNormal(log(700), log(4))
    bin_edges_area = colony_mean_area(bin_edges(; unit=:cm))

    # Find integral density between bounds of each size class areas to create weights for each size class.
    cdf_integral = cdf.(reef_mod_area_dist, bin_edges_area)
    size_class_weights = (cdf_integral[:, 2:end] .- cdf_integral[:, 1:(end - 1)])
    size_class_weights = size_class_weights ./ sum(size_class_weights; dims=2)

    # Take the mean over repeats, as suggested by YM (pers comm. 2023-02-27 12:40pm AEDT).
    # Convert from percent to relative values.
    icc_data = ((dropdims(mean(icc_data; dims=2); dims=2)) ./ 100.0)

    # Repeat species over each size class and reshape to give ADRIA compatible size
    # [(n_groups × n_sizes) ⋅ n_locs]
    # Multiply by size class weights to give initial cover distribution over each size class.
    n_sizes = size(size_class_weights, 2)
    n_groups = length(icc_files)
    n_locs = size(icc_data, 1)
    icc_data = icc_data .* reshape(size_class_weights, (1, n_groups, n_sizes))
    # Reshape and Permute icc_data from
    # [n_locs ⋅ n_groups ⋅ n_sizes] to [(n_groups × n_sizes) ⋅ n_locs]
    icc_data = reshape(
        permutedims(icc_data, (3, 2, 1)), (n_sizes * n_groups, n_locs)
    )

    # Convert values relative to absolute area to values relative to k area
    icc_data = _convert_abs_to_k(icc_data, loc_data)

    return DataCube(
        icc_data;
        species=1:(length(icc_files) * n_sizes),
        locs=loc_ids
    )
end

"""
    switch_RCPs!(d::RMEDomain, RCP::String)::Domain

Switch environmental datasets to represent the given RCP.
"""
function switch_RCPs!(d::RMEDomain, RCP::String)::RMEDomain
    @set! d.RCP = RCP
    data_files = joinpath(d.env_layer_md.dpkg_path, "data_files")

    timeframe::Tuple{Int64,Int64} = (
        d.env_layer_md.timeframe[1],
        d.env_layer_md.timeframe[end]
    )
    @set! d.dhw_scens = load_DHW(RMEDomain, data_files, RCP; timeframe=timeframe)

    # Cyclones are not RCP-specific?
    # @set! d.wave_scens = load_cyclones(RMEDomain, data_files, loc_ids)

    return d
end

function Base.show(io::IO, mime::MIME"text/plain", d::RMEDomain)::Nothing
    println("""
        ReefMod Engine Domain: $(d.name)

        Number of locations: $(n_locations(d))
        Timeframe: $(d.env_layer_md.timeframe[1]) - $(d.env_layer_md.timeframe[end])
    """)

    return nothing
end
