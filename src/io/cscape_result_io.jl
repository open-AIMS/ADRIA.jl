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

function load_results(::Type{CScapeResultSet}, result_loc::String)::CScapeResultSet
    !isdir(result_loc) ? error("Expected a directory but received $(result_loc)") : nothing
    
    result_path = _get_output_path(result_loc)
    raw_set = open_dataset(result_path)

    res_name::String = _get_result_name(raw_set)
    res_rcp::String = _get_rcp(raw_set)
    
    !haskey(raw_set.cubes, :reef_siteid) ? error("Unable to find location ids.") : nothing
    location_ids = _get_reefids(raw_set.cubes[:reef_siteid])

    gpkg_path = _get_gpkg_path(result_loc)
    geodata = GDF.read(gpkg_path)

    connectivity_path = joinpath(result_loc, "connectivity.csv")
    connectivity = CSV.read(connectivity_path, DataFrame, comment="#", header=true)
    
    # There is missing location data in site data. Use intersection
    gpkg_mask = BitVector([loc_name in location_ids for loc_name in geodata.reef_siteid])
    res_mask  = BitVector([loc_name in geodata.reef_siteid for loc_name in location_ids])
    conn_mask = BitVector([
        (loc_name in location_ids) && (loc_name in geodata.reef_siteid)
        for loc_name in names(connectivity)
    ])

    geodata = geodata[gpkg_mask, :]
    conn_sites = names(connectivity)[conn_mask] 
    connectivity = connectivity[conn_mask[2:end], conn_mask]
    
    # Re order data to match location ordering
    geo_id_order = [first(findall(x .== geodata.reef_siteid)) for x in location_ids[res_mask]]
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
    if !haskey(raw_set.cubes, :cover)
        @warn "Unable to find cover. Skipping)"
    else
        outcomes[:relative_taxa_cover] = create_relative_taxa_cover(raw_set.cover[reef_site=res_mask]) ./100
        outcomes[:relative_cover] = create_relative_cover(raw_set.cover[reef_site=res_mask]) ./100
        outcomes[:size_classes] = create_size_classes(raw_set.cover[reef_site=res_mask])
    end

    if !haskey(raw_set.cubes, :size_classes)
        @warn "Unable to find size_classes. Skipping"
    else
        outcomes[:size_classes]
    end
    #for outcome_key ∈ outcome_keys
    #    if !haskey(raw_set.cubes, outcome_key)
    #        @warn "Unable to find $(string(outcome_key). Skipping)"
    #        continue
    #    end
    #    if outcome_key == :cover
    #        continue
    #    end
    #    outcomes[outcome_key] = raw_set.cubes[outcome_key][reef_sites=res_mask]
    #end

    return CScapeResultSet(
        res_name,
        res_rcp,
        location_ids[res_mask],
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
    create_size_classes(size_class_cube::YAXArray)::YAXArray

Create outcome array describing the distribution of size classes.
"""
function create_size_classes(size_class_cube::YAXArray)::YAXArray
    return dropdims(sum(size_class_cube, dims=:intervened), dims=:intervened)
end

function create_relative_taxa_cover(cover_cube::YAXArray)::YAXArray
    cover_cube = rename_dims(cover_cube)
    rel_taxa_cover::YAXArray = dropdims(sum(
        cover_cube, dims=(:thermal_tolerance, :intervened)
    ), dims=(:thermal_tolerance, :intervened))
    return permutedims(rel_taxa_cover, [2, 3, 4, 1])
end

function create_relative_cover(cover_cube::YAXArray)::YAXArray
    rel_taxa_cover = create_relative_taxa_cover(cover_cube)
    return dropdims(sum(rel_taxa_cover, dims=:taxa), dims=:taxa)
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

Reef IDs are stored in a space seperated list in the properties of reef_siteid cube.
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

function short_habitat_to_long(encoding::SubString{String})::String
    if (encoding == "C")
        return "Crest"
    elseif (encoding == "OF")
        return "Outer Flat" 
    elseif (encoding == "S")
        return "Slope"
    elseif (encoding == "SS")
        return "Sheltered Slope"
    else 
        return ""
    end
end

"""
    rename_dims(cscape_cube::YAXArray)::Nothing

Rename the names of the dimensions to align with ADRIA's expected dimension names.
"""
function rename_dims(cscape_cube::YAXArray)::YAXArray
    dim_names = name.(cscape_cube.axes)
    if :draws in dim_names
        cscape_cube = renameaxis!(cscape_cube, :draws=>:scenarios)
    end
    if :year in dim_names
        cscape_cube = renameaxis!(cscape_cube, :year=>:timesteps)
    end
    if :reef_sites in dim_names
        cscape_cube = renameaxis!(cscape_cube, :reef_sites=>:sites)
    end
    if :ft in dim_names
        cscape_cube = renameaxis!(cscape_cube, :ft=>:taxa)
    end
    return cscape_cube
end

"""
    equivalient_id(short::String, long::String)::Bool

Location IDs stored in connectivity data files have a different more verbose format.
"""
function equivalent_id(short::String, long::String)::Bool
    short_frags = split(short, '_')
    long_frags = split(long, '_')
    length(short_frags) != length(long_frags) ? error("Unexpected Reef ID format") : nothing
    if short_frags[1] != long_frags[1]
        return false
    elseif short_habitat_to_long(short_frags[3]) != long_frags[3]
        return false
    elseif short_frags[4] != long_frags[4] && short_frags[4] != long_frags[4][1:end-1]
        return false
    end
    return true
end

"""
    conn_ids_format(loc_ids::Vector{String})::Vector{String}

Reef IDs in CScape Connectivity files are named differently, so we need to create an equivalent reef ID mask for data extraction.
"""
function conn_ids_format(connectivity_path::String, loc_ids::Vector{String})::Vector{String}
    fn = filter(x -> occursin("csv", x), readdir(connectivity_path))[1]
    conn_loc_ids = names(CSV.read(fn, DataFrame, comment="#", header=true))[2:end]
    loc_ordering = [findfirst(equivalent_id.(loc_ids, c_id)) for c_id in conn_loc_ids]
end

function Base.show(io::IO, mime::MIME"text/plain", rs::CScapeResultSet)
    rcps = rs.RCP
    scens = length(rs.outcomes[:relative_cover].scenarios)
    sites = length(rs.outcomes[:relative_cover].sites)
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

function reformat_temporal_inputs(temporal_input::DataFrame, RCP::String="1")::YAXArray
    rcp::Int = parse(Int, RCP)
    masked_temp_input = temporal_input[temporal_input.RCP == rcp, :]
    year_dim = unique(masked_temp_input.year)
    loc_dim = unique(masked_temp_input.reef_siteid)
end
