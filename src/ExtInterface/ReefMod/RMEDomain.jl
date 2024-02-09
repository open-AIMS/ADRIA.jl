using ADRIA: Domain, GDF, SimConstants
using CSV
using DataFrames
using Distributions
using ModelParameters
using Statistics

mutable struct RMEDomain <: Domain
    const name::String
    RCP::String
    env_layer_md::ADRIA.EnvLayer
    scenario_invoke_time::String  # time latest set of scenarios were run
    const conn::YAXArray{Float64}
    const in_conn::Vector{Float64}
    const out_conn::Vector{Float64}
    const strong_pred::Vector{Int64}
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
    data_files = joinpath(fn_path, "data_files")
    dhw_scens = load_DHW(RMEDomain, data_files, RCP)
    loc_ids = axiskeys(dhw_scens)[2]

    site_data_path = joinpath(data_files, "region", "reefmod_gbr.gpkg")
    site_data = GDF.read(site_data_path)
    site_id_col = "LOC_NAME_S"
    cluster_id_col = "LOC_NAME_S"
    site_ids::Vector{String} = site_data[:, site_id_col]

    id_list = CSV.read(
        joinpath(data_files, "id", "id_list_2023_03_30.csv"),
        DataFrame;
        header=false,
        comment="#",
    )

    # Re-order spatial data to match RME dataset
    # MANUAL CORRECTION
    site_data[site_data.LABEL_ID.=="20198", :LABEL_ID] .= "20-198"
    id_order::Vector{Int64} = [first(findall(x .== site_data.LABEL_ID)) for x in string.(id_list[:, 1])]
    site_data = site_data[id_order, :]

    # Check that the two lists of location ids are identical
    @assert isempty(findall(site_data.LABEL_ID .!= id_list[:, 1]))

    # Convert area in km² to m²
    site_data[:, :area] .= id_list[:, 2] .* 1e6

    # Calculate `k` area (1.0 - "ungrazable" area)
    site_data[:, :k] .= 1.0 .- id_list[:, 3]

    # Need to load initial coral cover after we know `k` area.
    init_coral_cover = load_initial_cover(RMEDomain, data_files, loc_ids, site_data)

    conn_data = load_connectivity(RMEDomain, data_files, loc_ids)
    in_conn, out_conn, strong_pred = ADRIA.connectivity_strength(
        conn_data, vec(site_data.area .* site_data.k), similar(conn_data)
    )

    # Set all site depths to 6m below sea level
    # (ReefMod does not account for depth)
    # Ensure column is of float type
    site_data[:, :depth_med] .= 6.0
    site_data[!, :depth_med] = convert.(Float64, site_data[!, :depth_med])

    # Add GBRMPA zone type info as well
    gbr_zt_path = joinpath(data_files, "region", "gbrmpa_zone_type.csv")
    gbr_zone_types = CSV.read(gbr_zt_path, DataFrame; types=String)
    missing_rows = ismissing.(gbr_zone_types[:, "GBRMPA Zone Types"])
    gbr_zone_types[missing_rows, "GBRMPA Zone Types"] .= ""
    zones = gbr_zone_types[:, "GBRMPA Zone Types"]
    zones = replace.(zones, "Zone" => "", " " => "")
    site_data[:, :zone_type] .= zones

    # This loads cyclone categories, not mortalities, so ignoring for now.
    # cyc_scens = load_cyclones(RMEDomain, data_files, loc_ids, timeframe)

    timeframe = (2022, 2100)
    timeframe_range = timeframe[1]:timeframe[2]

    dim_timesteps = Dim{:timesteps}(timeframe_range)
    dim_locs = Dim{:locs}(loc_ids)
    dim_species = Dim{:species}(1:6)
    dim_scenarios = Dim{:scenarios}([1])

    wave_axlist = (dim_timesteps, dim_locs, dim_scenarios)
    wave_data = zeros(length(timeframe_range), nrow(site_data), 1)
    wave_scens::YAXArray{Float64} = YAXArray(wave_axlist, wave_data)

    cyc_axlist = (dim_timesteps, dim_locs, dim_species, dim_scenarios)
    cyc_data = zeros(length(timeframe_range), nrow(site_data), 6, 1)
    cyc_scens::YAXArray{Float64} = YAXArray(cyc_axlist, cyc_data)

    env_md = EnvLayer(
        fn_path,
        site_data_path,
        site_id_col,
        cluster_id_col,
        "",
        "",
        "",
        "",
        timeframe_range,
    )

    model::Model = Model((
        EnvironmentalLayer(dhw_scens, wave_scens, cyc_scens),
        Intervention(),
        CriteriaWeights(),
        Coral()
    ))

    return RMEDomain(
        "ReefMod Engine",
        RCP,
        env_md,
        "",
        conn_data,
        in_conn,
        out_conn,
        strong_pred,
        site_data,
        site_id_col,
        cluster_id_col,
        init_coral_cover,
        CoralGrowth(nrow(site_data)),
        site_ids,
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
    load_DHW(::Type{RMEDomain}, data_path::String, rcp::String, timeframe=(2022, 2100))

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
)::NamedDimsArray
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

    axlist = (
        Dim{:timesteps}(timeframe[1]:timeframe[2]),
        Dim{:locs}(loc_ids),
        Dim{:scenarios}(rcp_files[keep_ds]),
    )
    # Only return valid scenarios
    return YAXArray(axlist, data_cube[:, :, keep_ds])
end

"""
    load_connectivity(::Type{RMEDomain}, data_path::String, loc_ids::Vector{String})

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
)::NamedDimsArray
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
    return YAXArray((Dim{:Source}(loc_ids), Dim{:Sink}(loc_ids)), conn_data)
end

"""
    load_cyclones(::Type{RMEDomain}, data_path::String, loc_ids::Vector{String}, tf::Tuple{Int64, Int64})::NamedDimsArray

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

    axlist = (Dim{:timesteps}(tf[1]:tf[2]), Dim{:locs}(loc_ids), Dim{:scenarios}(1:length(cyc_files)))
    return YAXArray(axlist, cyc_data)
end

"""
    load_initial_cover(::Type{RMEDomain}, data_path::String, loc_ids::Vector{String}, site_data::DataFrame)::NamedDimsArray

# Arguments
- `RMEDomain`
- `data_path` : path to ReefMod data
- `loc_ids` : location ids

# Returns
YAXArray[locs, species]
"""
function load_initial_cover(
    ::Type{RMEDomain}, data_path::String, loc_ids::Vector{String}, site_data::DataFrame
)::NamedDimsArray
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

    axlist = (Dim{:species}(1:(length(icc_files)*6)), Dim{:locs}(loc_ids))
    return YAXArray(axlist, icc_data)
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
    # loc_ids = axiskeys(d.dhw_scens)[2]
    # @set! d.wave_scens = load_cyclones(RMEDomain, data_files, loc_ids)

    return d
end

function Base.show(io::IO, mime::MIME"text/plain", d::RMEDomain)::Nothing
    # df = model_spec(d)
    println("""
        ReefMod Engine Domain: $(d.name)

        Number of locations: $(n_locations(d))
    """)
    return nothing
end
