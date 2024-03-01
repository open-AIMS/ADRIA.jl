using ADRIA: EnvLayer, GDF, ZeroDataCube, DataCube
using ArchGDAL: centroid
using DataFrames
using YAXArrays
using CSV

struct CScapeResultSet <: ResultSet
    name::String
    RCP::String

    site_ids
    site_area::Vector{Float64}
    site_max_coral_cover::Vector{Float64}
    site_centroids
    env_layer_md::EnvLayer
    connectivity_data
    site_data

    sim_constants

    # raw::AbstractArray
    outcomes
end


function load_results(::Type{CScapeResultSet}, result_loc::String)
    !isdir(result_loc) ? error("Expected a directory but received $(result_loc)") : nothing
    
    result_path = _get_output_path(result_loc)
    raw_set = open_dataset(result_path)

    res_name::String = _get_result_name(raw_set)
    res_rcp::String = _get_rcp(raw_set)
    
    !haskey(raw_set.cubes, :reef_siteid) ? error("Unable to find location ids.") : nothing
    location_ids = _get_reefids(raw_set.cubes[:reef_siteid])
    # "Moore_MR_C_1" -> "MR_C_1"
    short_loc_ids = [id[(findfirst('_', id)+1):end] for id in location_ids]

    gpkg_path = _get_gpkg_path(result_loc)
    geodata = GDF.read(gpkg_path)

    connectivity_path = joinpath(result_loc, "connectivity", "connectivity.csv")
    connectivity = CSV.read(connectivity_path, DataFrame, comment="#", header=false)
    
    # There is missing location data in site data. Use intersection
    gpkg_mask = BitVector([loc_name in short_loc_ids for loc_name in geodata.site_id])
    res_mask  = BitVector([loc_name in geodata.site_id for loc_name in short_loc_ids])
    conn_mask = BitVector(
        [!ismissing(loc_name) &&
         (loc_name in location_ids) && 
         (loc_name[findfirst('_', loc_name)+1:end] in geodata.site_id)
        for loc_name in connectivity[1, :]]
    )

    geodata = geodata[gpkg_mask, :]
    conn_sites = connectivity[1, :][conn_mask] 
    connectivity = connectivity[conn_mask, conn_mask]

    geo_id_order = [first(findall(x .== geodata.site_id)) for x in short_loc_ids[res_mask]]
    conn_id_order = [first(findall(x .== Vector(conn_sites))) for x in location_ids[res_mask]]
    
    geodata = geodata[geo_id_order, :]
    connectivity = connectivity[conn_id_order, conn_id_order]

    timeframe = 2009:2099
    if !haskey(raw_set.properties, "temporal_range")
        @warn "Unable to find timeframe defaulting to $(timeframe[1]):$(timeframe[end])"
    else
        tf_str = split(raw_set.properties["temporal_range"], ":")
        timeframe = parse(Int, tf_str[1]):parse(Int, tf_str[2])
    end

    env_layer_md::EnvLayer = EnvLayer(
        result_loc,
        gpkg_path,
        "site_id",
        "Reef",
        "",
        connectivity_path,
        "",
        "",
        timeframe
    )

    location_max_coral_cover = 1 .- geodata.k ./ 100
    location_centroids = [centroid(multipoly) for multipoly ∈ geodata.geom]

    outcome_keys = [
        :larvae, 
        :external_larvae, 
        :internal_received_larvae, 
        :settlers, 
        :eggs, 
        :size_classes,
        :cover
    ]
    outcomes = Dict{Symbol, YAXArray}()
    for outcome_key ∈ outcome_keys
        if !haskey(raw_set.cubes, outcome_key)
            @warn "Unable to find $(string(outcome_key). Skipping)"
            continue
        end
        if outcome_key == :cover
            outcomes[:relative_cover] = raw_set.cubes[outcome_key]
            continue
        end
        outcomes[outcome_key] = raw_set.cubes[outcome_key]
    end

    return CScapeResultSet(
        res_name,
        res_rcp,
        geodata.site_id,
        geodata.area,
        location_max_coral_cover,
        location_centroids,
        env_layer_md,
        connectivity,
        geodata,
        SimConstants(),
        outcomes,
    )
end

"""
    _get_output_path(result_loc::String)::String

Get the filename of the netcdf file contained in the result subdirectory. Must have a NetCDF
suffix and begin with `NetCDF_Scn`.
"""
function _get_output_path(result_loc::String)::String
    res_subdir = joinpath(result_loc, "results")
    possible_files = filter(isfile, readdir(res_subdir, join=true))
    possible_files = filter(x -> occursin(r"NetCDF_Scn.*.nc", x), possible_files)
    if length(possible_files) == 0
        error("Unable to find result netcdf file in subdirectory $(res_subdir)")
    elseif length(possible_files) > 1
        @warn "Found multiple netcdf files, using first."
    end
    return possible_files[1]
end

"""
    _get_gpkg_path(result_loc::String)

Get the path to the gpkg file contained in the site_data subdirectory. 
"""
function _get_gpkg_path(result_loc::String)
    gpkg_dir = joinpath(result_loc, "site_data")
    possible_files = filter(isfile, readdir(gpkg_dir, join=true))
    possible_files = filter(x -> occursin(".gpkg", x), possible_files)
    if length(possible_files) == 0
        error("Unable to find site data in $(gpkg_dir)")
    elseif length(possible_files) > 1
        @warn "Found multiple gpkg files, using first."
    end
    return possible_files[1]
end
"""
    _get_result_name()::String

Get the name of the data set from the properties of the dataset.
"""
function _get_result_name(ds::Dataset)::String
    name = "CScape Results"
    if !haskey(ds.properties, "title")
        @warn "Unable to find key `title` in dataset properties, \
               defaulting to `CScape Results`"
    else
        name = ds.properties["title"]
    end
    return name
end

"""
    _get_rcp(ds::Datset)::String

Get the RCP from the ssp field in dataset properties. Assume RCP is last two numbers.
"""
function _get_rcp(ds::Dataset)::String
    rcp = ""
    if !haskey(ds.properties, "ssp")
        @warn "Unable to find key `ssp` in dataset properties."
    else
        rcp = ds.properties["ssp"][end-1:end]
    end
    return rcp
end

"""
    _get_reefids(ds::YAXArray)::String

Reef IDs are stored in a space seperated list in the properties of reef_siteids cube.
"""
function _get_reefids(reef_cube::YAXArray)::Vector{String}
    # Site IDs are necessary to extract the correct data from the geopackage
    if !haskey(reef_cube.properties, "flag_meanings")
        error("Unable to find key `flag_meanings` in Cube properties.")
    end

    reef_ids = split(reef_cube.properties["flag_meanings"], " ")
    if (reef_ids[1] == reef_ids[2])
        return reef_ids[2:end] # Possible first element duplication
    end 
    return reef_ids
end


function Base.show(io::IO, mime::MIME"text/plain", rs::CScapeResultSet)
    rcps = rs.RCP
    scens = length(rs.outcomes[:relative_cover].draws)
    sites = length(rs.outcomes[:relative_cover].reef_sites)
    tf = rs.env_layer_md.timeframe

    println("""
    Name: $(rs.name)

    Results stored at: $(rs.env_layer_md.dpkg_path)

    RCP(s) represented: $(rcps)
    Scenarios run: $(scens)
    Number of sites: $(sites)
    Timesteps: $(tf)
    """)
end
