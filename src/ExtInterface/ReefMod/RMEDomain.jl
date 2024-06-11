using Distributions
using Statistics

using CSV, DataFrames, ModelParameters

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
    const site_data::DataFrames.DataFrame
    const site_id_col::String
    const cluster_id_col::String
    init_coral_cover::YAXArray{Float64}
    const coral_growth::CoralGrowth
    const site_ids::Vector{String}
    dhw_scens::YAXArray{Float64}

    # `wave_scens` holds empty wave data to maintain compatibility with
    # ADRIA's dMCDA methods
    wave_scens::YAXArray{Float64}

    # `cyclone_mortality_scens` holds dummy cyclone mortality data to maintain compatibility
    # with ADRIA's dMCDA methods
    cyclone_mortality_scens::YAXArray{Float64}

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
        rename!(spatial_data, Dict("LOC_NAME_S"=>"cluster_id"))
    catch err
        if !(err isa ArgumentError)
            rethrow(err)
        end

        rename!(spatial_data, Dict("reef_name"=>"cluster_id"))
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
        id_order = [first(findall(x .== spatial_data.LABEL_ID)) for x in string.(id_list[:, 1])]
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
    load_domain(::Type{RMEDomain}, fn_path::String, RCP::String)::RMEDomain

Load a ReefMod Engine dataset.

# Arguments
- `RMEDomain` : DataType
- `fn_path` : path to ReefMod Engine dataset
- `RCP` : Representative Concentration Pathway scenario ID

# Returns
RMEDomain
"""
function load_domain(::Type{RMEDomain}, fn_path::String, RCP::String)::RMEDomain
    isdir(fn_path) ? true : error("Path does not exist or is not a directory.")

    data_files = joinpath(fn_path, "data_files")
    dhw_scens::YAXArray{Float64} = load_DHW(RMEDomain, data_files, RCP)
    loc_ids::Vector{String} = collect(dhw_scens.locs)

    # Load the geopackage file
    gpkg_path = joinpath(data_files, "region", "reefmod_gbr.gpkg")
    spatial_data = GDF.read(gpkg_path)

    # Standardize IDs to use for site/reef and cluster
    _standardize_cluster_ids!(spatial_data)

    reef_id_col = "UNIQUE_ID"
    cluster_id_col = "UNIQUE_ID"

    reef_ids::Vector{String} = spatial_data[:, reef_id_col]

    # Load accompanying ID list
    # TODO: Create canonical geopackage file that aligns all IDs.
    #       Doing this removes the need for the manual correction below and removes the
    #       dependency on this file.
    id_list = CSV.read(
        joinpath(data_files, "id", "id_list_2023_03_30.csv"),
        DataFrame;
        header=false,
        comment="#"
    )

    _manual_id_corrections!(spatial_data, id_list)

    # Check that the two lists of location ids are identical
    try
        @assert isempty(findall(spatial_data.LABEL_ID .!= id_list[:, 1]))
    catch
        @assert isempty(findall(spatial_data.RME_GBRMPA_ID .!= id_list[:, 1]))
    end

    # Convert area in km² to m²
    spatial_data[:, :area] .= id_list[:, 2] .* 1e6

    # Calculate `k` area (1.0 - "ungrazable" area)
    spatial_data[:, :k] .= 1.0 .- id_list[:, 3]

    # Need to load initial coral cover after we know `k` area.
    init_coral_cover::YAXArray{Float64} = load_initial_cover(RMEDomain, data_files, loc_ids, spatial_data)

    conn_data::YAXArray{Float64} = load_connectivity(RMEDomain, data_files, loc_ids)

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
    spatial_data[:, :zone_type] .= zones

    # This loads cyclone categories, not mortalities, so ignoring for now.
    # cyc_scens = load_cyclones(RMEDomain, data_files, loc_ids, timeframe)

    timeframe = (2022, 2100)
    timeframe_range = timeframe[1]:timeframe[2]

    wave_scens::YAXArray{Float64} = ZeroDataCube(;
        T=Float64, timesteps=timeframe_range, locs=loc_ids, scenarios=[1]
    )

    functional_groups = functional_group_names()
    cyc_scens::YAXArray{Float64} = ZeroDataCube(;
        T=Float64, timesteps=timeframe_range, locs=loc_ids, species=functional_groups, scenarios=[1]
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
        timeframe_range,
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
        Coral()
    ))

    return RMEDomain(
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
        model,
        SimConstants(),
    )
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
    ::Type{RMEDomain}, data_path::String, rcp::String, timeframe=(2022, 2100)
)::YAXArray
    dhw_path = joinpath(data_path, "dhw")
    rcp_files = _get_relevant_files(dhw_path, rcp)
    rcp_files = filter(x -> occursin("SSP", x), rcp_files)
    if isempty(rcp_files)
        ArgumentError("No DHW data files found in: $(dhw_path)")
    end

    first_file = CSV.read(rcp_files[1], DataFrame)
    loc_ids = String.(first_file[:, 1])

    data_tf = parse.(Int64, names(first_file[:, 2:end]))
    tf_start = findall(timeframe[1] .∈ data_tf)[1]
    tf_end = findall(timeframe[2] .∈ data_tf)[1]

    d1 = first_file[:, (tf_start+1):(tf_end+1)]
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
            tf_start = findall(timeframe[1] .∈ data_tf)[1]
            tf_end = findall(timeframe[2] .∈ data_tf)[1]
        catch err
            if !(err isa BoundsError)
                rethrow(err)
            end

            @warn "Building DHW: Could not find matching time frame, skipping $rcp_data_fn"
            keep_ds[i] = false  # mark scenario for removal
            continue
        end

        data_cube[:, :, i+1] .= Matrix(d[:, tf_start:tf_end])'
    end

    # Only return valid scenarios
    return DataCube(
        data_cube[:, :, keep_ds];
        timesteps=timeframe[1]:timeframe[2],
        locs=loc_ids,
        scenarios=rcp_files[keep_ds],
    )
end

"""
    load_connectivity(::Type{RMEDomain}, data_path::String, loc_ids::Vector{String})::YAXArray

Loads the average connectivity matrix.

# Arguments
- `RMEDomain`
- `data_path` : path to ReefMod data
- `loc_ids` : location ids

# Returns
YAXArray[source, sinks]
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
        Sink=loc_ids,
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
    tf::Tuple{Int64,Int64},
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
    cyc_data = permutedims(cyc_data, (2, 1, 3))[1:((tf[2]-tf[1])+1), :, :]
    return DataCube(
        cyc_data;
        timesteps=tf[1]:tf[2],
        locs=loc_ids,
        scenarios=1:length(cyc_files),
    )
end

"""
    load_initial_cover(::Type{RMEDomain}, data_path::String, loc_ids::Vector{String}, site_data::DataFrame)::YAXArray

# Arguments
- `RMEDomain`
- `data_path` : path to ReefMod data
- `loc_ids` : location ids

# Returns
YAXArray[locs, species]
"""
function load_initial_cover(
    ::Type{RMEDomain}, data_path::String, loc_ids::Vector{String}, site_data::DataFrame
)::YAXArray
    icc_path = joinpath(data_path, "initial")
    icc_files = _get_relevant_files(icc_path, "coral_")
    if isempty(icc_files)
        ArgumentError("No coral cover data files found in: $(icc_path)")
    end

    # Shape is locations, scenarios, species
    icc_data = zeros(length(loc_ids), 20, length(icc_files))
    for (i, fn) in enumerate(icc_files)
        icc_data[:, :, i] = Matrix(
            CSV.read(fn, DataFrame; drop=[1], header=false, comment="#")
        )
    end

    # Use ReefMod distribution for coral size class population (shape parameters have units log(cm^2))
    # as suggested by YM (pers comm. 2023-08-08 12:55pm AEST). Distribution is used to split ReefMod initial
    # species covers into ADRIA's 6 size classes by weighting with the cdf.
    reef_mod_area_dist = LogNormal(log(700), log(4))
    bin_edges_area = colony_mean_area(Float64[0, 2, 5, 10, 20, 40, 80])

    # Find integral density between bounds of each size class areas to create weights for each size class.
    cdf_integral = cdf.(reef_mod_area_dist, bin_edges_area)
    size_class_weights = (cdf_integral[2:end] .- cdf_integral[1:(end-1)])
    size_class_weights = size_class_weights ./ sum(size_class_weights)

    # Take the mean over repeats, as suggested by YM (pers comm. 2023-02-27 12:40pm AEDT).
    # Convert from percent to relative values.
    icc_data = ((dropdims(mean(icc_data; dims=2); dims=2)) ./ 100.0)

    # Repeat species over each size class and reshape to give ADRIA compatible size (36 * n_locs).
    # Multiply by size class weights to give initial cover distribution over each size class.
    icc_data = Matrix(hcat(reduce.(vcat, eachrow(icc_data .* [size_class_weights]))...))

    # Convert values relative to absolute area to values relative to k area
    icc_data = _convert_abs_to_k(icc_data, site_data)

    return DataCube(
        icc_data;
        species=1:(length(icc_files)*6),
        locs=loc_ids,
    )
end

"""
    switch_RCPs!(d::RMEDomain, RCP::String)::Domain

Switch environmental datasets to represent the given RCP.
"""
function switch_RCPs!(d::RMEDomain, RCP::String)::RMEDomain
    @set! d.RCP = RCP
    data_files = joinpath(d.env_layer_md.dpkg_path, "data_files")
    @set! d.dhw_scens = load_DHW(RMEDomain, data_files, RCP)

    # Cyclones are not RCP-specific?
    # @set! d.wave_scens = load_cyclones(RMEDomain, data_files, loc_ids)

    return d
end

function Base.show(io::IO, mime::MIME"text/plain", d::RMEDomain)::Nothing
    println("""
        ReefMod Engine Domain: $(d.name)

        Number of locations: $(n_locations(d))
    """)

    return nothing
end
