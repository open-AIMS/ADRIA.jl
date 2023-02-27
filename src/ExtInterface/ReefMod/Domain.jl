using NamedDims, AxisKeys, CSV
using CSV, DataFrames, Statistics
import GeoDataFrames as GDF

using ModelParameters
using ADRIA: SimConstants, Domain, site_distances


mutable struct ReefModDomain <: Domain
    const name::String
    RCP::String
    const TP_data
    const in_conn
    const out_conn
    const strong_pred
    const site_data
    const site_distances
    const median_site_distance
    const site_id_col
    const unique_site_id_col
    init_coral_cover
    const site_ids
    dhw_scens

    # `wave_scens` Actually holds cyclones, but to maintain compatibility with
    # ADRIA's dMCDA methods
    wave_scens

    model
    sim_constants::SimConstants
end


"""
    _get_relevant_files(fn_path, ident)

Filter files found in given path down to those that have the provided identifier.

# Arguments
- `fn_path` : directory of files
- `ident` : filter out files which do not have this in their filename.
"""
function _get_relevant_files(fn_path::String, ident::String)
    valid_files = filter(isfile, readdir(fn_path; join=true))
    return filter(x -> occursin(ident, x), valid_files)
end

"""
    load_DHW(::Type{ReefModDomain}, data_path::String, rcp::String, timeframe=(2022, 2100))

Loads ReefMod DHW data as a datacube.

# Arguments
- `ReefModDomain`
- `data_path` : path to ReefMod data
- `rcp` : RCP identifier
- `timeframe` : range of years to represent.
"""
function load_DHW(::Type{ReefModDomain}, data_path::String, rcp::String, timeframe=(2022, 2100))::NamedDimsArray
    dhw_path = joinpath(data_path, "dhw")
    rcp_files = _get_relevant_files(dhw_path, rcp)
    if length(rcp_files) == 0
        ArgumentError("No DHW data files found in: $(dhw_path)")
    end

    first_file = CSV.read(rcp_files[1], DataFrame)
    loc_ids = String.(first_file[:, 1])

    data_tf = parse.(Int64, names(first_file[:, 2:end]))
    tf_start = findall(timeframe[1] .∈ data_tf)[1]
    tf_end = findall(timeframe[2] .∈ data_tf)[1]

    d1 = first_file[:, tf_start+1:tf_end+1]
    data_shape = reverse(size(d1))
    data_cube = zeros(data_shape..., length(rcp_files))
    data_cube[:, :, 1] .= Matrix(d1)'

    local tf_start
    local tf_end
    keep_ds = fill(true, length(rcp_files))
    for (i, rcp_data_fn) in enumerate(rcp_files[2:end])
        d = CSV.read(rcp_data_fn, DataFrame)[:, 2:end]

        if size(d, 1) == 0
            @info "empty file?" rcp_data_fn
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
            keep_ds[i] = false  # mark member for removal
            continue
        end

        data_cube[:, :, i+1] .= Matrix(d[:, tf_start:tf_end])'
    end

    # Only return valid members
    return NamedDimsArray(data_cube[:, :, keep_ds], timesteps=timeframe[1]:timeframe[2], locs=loc_ids, member=rcp_files[keep_ds])
end

"""
    load_connectivity(::Type{ReefModDomain}, data_path::String, loc_ids::Vector{String})

Loads the average connectivity matrix.

# Arguments
- `ReefModDomain`
- `data_path` : path to ReefMod data
- `loc_ids` : location ids
"""
function load_connectivity(::Type{ReefModDomain}, data_path::String, loc_ids::Vector{String})::NamedDimsArray
    conn_path = joinpath(data_path, "con_bin")
    conn_files = _get_relevant_files(conn_path, "CONNECT_ACRO")
    if length(conn_files) == 0
        ArgumentError("No DHW data files found in: $(conn_path)")
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
    conn_data = dropdims(mean(tmp_mat, dims=3), dims=3)
    return NamedDimsArray(conn_data, source=loc_ids, sinks=loc_ids)
end

"""
    load_cyclones(::Type{ReefModDomain}, data_path::String, loc_ids::Vector{String})::NamedDimsArray

# Arguments
- `ReefModDomain`
- `data_path` : path to ReefMod data
- `loc_ids` : location ids

# Returns
NamedDimsArray[locations, years, members]
"""
function load_cyclones(::Type{ReefModDomain}, data_path::String, loc_ids::Vector{String})::NamedDimsArray
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
        cyc_data[:, :, i] = Matrix(CSV.read(fn, DataFrame; drop=[1], header=false))
    end

    # Mean over all years
    return NamedDimsArray(cyc_data, locs=loc_ids, years=1:num_years, members=1:length(cyc_files))
end

"""
    load_initial_cover(::Type{ReefModDomain}, data_path::String, loc_ids::Vector{String})::NamedDimsArray

# Arguments
- `ReefModDomain`
- `data_path` : path to ReefMod data
- `loc_ids` : location ids

# Returns
NamedDimsArray[locations, species, members]
"""
function load_initial_cover(::Type{ReefModDomain}, data_path::String, loc_ids::Vector{String})::NamedDimsArray
    icc_path = joinpath(data_path, "initial")
    icc_files = _get_relevant_files(icc_path, "coral_")
    if length(icc_files) == 0
        ArgumentError("No cyclone data files found in: $(icc_path)")
    end

    # Shape is locations, members, species
    icc_data = zeros(length(loc_ids), 20, length(icc_files))
    for (i, fn) in enumerate(icc_files)
        icc_data[:, :, i] = Matrix(CSV.read(fn, DataFrame; drop=[1], header=false))
    end

    # Reorder dims to: locations, species, members
    return NamedDimsArray(permutedims(icc_data, (1, 3, 2)), locs=loc_ids, species=1:length(icc_files), members=1:20)
end


"""
    load_domain(::Type{ReefModDomain}, fn_path, RCP)::ReefModDomain

Load a Domain for use with ReefMod.

# Arguments
- `ReefModDomain`
- `fn_path`
- `RCP`
"""
function load_domain(::Type{ReefModDomain}, fn_path::String, RCP::String)::ReefModDomain
    data_files = joinpath(fn_path, "data_files")
    dhw_scens = load_DHW(ReefModDomain, data_files, RCP)
    loc_ids = axiskeys(dhw_scens)[2]

    conn_data = load_connectivity(ReefModDomain, data_files, loc_ids)
    in_conn, out_conn, strong_pred = ADRIA.connectivity_strength(conn_data)

    site_data = GDF.read(joinpath(data_files, "region", "reefmod_gbr.gpkg"))
    site_dist, med_site_dist = ADRIA.site_distances(site_data)
    site_id_col = "LOC_NAME_S"
    unique_site_id_col = "LOC_NAME_S"
    init_coral_cover = load_initial_cover(ReefModDomain, data_files, loc_ids)
    site_ids = site_data[:, unique_site_id_col]

    cyc_scens = load_cyclones(ReefModDomain, data_files, loc_ids)

    model::Model = Model((EnvironmentalLayer(dhw_scens, cyc_scens), Intervention(), Criteria()))

    return ReefModDomain(
        "ReefMod", RCP,
        conn_data, in_conn, out_conn, strong_pred,
        site_data, site_dist, med_site_dist,
        site_id_col, unique_site_id_col,
        init_coral_cover,
        site_ids,
        dhw_scens, cyc_scens,
        model, SimConstants())
end


function Base.show(io::IO, mime::MIME"text/plain", d::ReefModDomain)

    # df = model_spec(d)
    println("""
    ReefMod Domain: $(d.name)

    Number of sites: $(n_locations(d))
    """)
end

