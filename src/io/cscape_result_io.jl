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
    scenario_groups

    sim_constants

    # raw::AbstractArray
    outcomes
    # Cscape uses different size classes
    coral_size_diameter::YAXArray
end

"""
    load_results(::Type{CScapeResultSet}, result_loc::String)::CScapeResultSet

Interface for loading CScape model outputs
"""
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

    connectivity_path = joinpath(result_loc, "connectivity/connectivity.csv")
    connectivity = CSV.read(connectivity_path, DataFrame, comment="#", header=true)
    
    # There is missing location data in site data. Use intersection
    gpkg_mask = BitVector([loc_name in location_ids for loc_name in geodata.reef_siteid])
    res_mask  = BitVector([loc_name in geodata.reef_siteid for loc_name in location_ids])
    conn_mask = BitVector([
        (loc_name in location_ids) && (loc_name in geodata.reef_siteid)
        for loc_name in names(connectivity)
    ])

    raw_set.properties["location_mask"] = res_mask

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
    for key in outcome_keys
        if !haskey(raw_set.cubes, key)
            @warn "Unable to find "*string(key)*". Skipping."
            continue
        end
        outcomes[key] = reformat_cube(raw_set.cubes[key], res_mask)
    end
    # add precomputed metrics for comptability
    outcomes[:relative_juveniles] = _cscape_relative_juvenile(raw_set)
    outcomes[:relative_cover] = _cscape_relative_cover(raw_set)
    outcomes[:relative_taxa_cover] = _cscape_relative_taxa_cover(raw_set, geodata.area)
    scen_groups = Dict(:counterfactual=>BitVector(true for _ in raw_set.draws))

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
        scen_groups,
        SimConstants(),
        outcomes,
        reformat_cube(raw_set.coral_size_diameter)
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

"""
    reformat_cube(cscape_cube::YAXArray)::YAXArray

Rename reorder the names of the dimensions to align with ADRIA's expected dimension names.
"""
function reformat_cube(cscape_cube::YAXArray, loc_mask::BitVector)
    cscape_cube = reformat_cube(cscape_cube)
    if :sites ∈ name.(cscape_cube.axes)
        cscape_cube = cscape_cube[sites=loc_mask]
    end
    return cscape_cube
end
function reformat_cube(cscape_cube::YAXArray)::YAXArray
    dim_names = name.(cscape_cube.axes)
    cscape_names = [:year, :reef_sites, :ft, :draws]
    adria_names = [:timesteps, :sites, :taxa, :scenarios]
    final_ordering::Vector{Int} = Vector{Int}(undef, length(dim_names))
    current_index = 1
    # rename expected dimensions and update the permutation vector
    for (cscape_nm, adria_nm) in zip(cscape_names, adria_names)
        if cscape_nm ∉ dim_names
            continue
        end
        cscape_cube = renameaxis!(cscape_cube, cscape_nm=>adria_nm)
        final_ordering[current_index] = findfirst(x->x==cscape_nm, dim_names)
        current_index += 1
    end
     # append the indices of non-adria axis to end of permutation vector
    for dim_name in dim_names
        if dim_name in cscape_names
            continue
        end
        final_ordering[current_index] = findfirst(x->x==dim_name, dim_names)
        current_index += 1;
    end
    if haskey(cscape_cube.properties, "units")
        if cscape_cube.properties["units"] == "percent"
            cscape_cube = cscape_cube ./ 100
            cscape_cube.properties["units"] == "proportion [0, 1]"
        end
    end
    cscape_cube = permutedims(cscape_cube, final_ordering)
    return cscape_cube
end

"""
    _cscape_relative_cover(dataset::Dataset)::YAXArray

Calculate relative cover metric for cscape data.
"""
function _cscape_relative_cover(dataset::Dataset)::YAXArray
    loc_mask = dataset.properties["location_mask"]
    cube = reformat_cube(dataset.cover, loc_mask)
    dim_sum = (:thermal_tolerance, :intervened, :taxa)
    return dropdims(sum(cube, dims=dim_sum), dims=dim_sum)
end

"""
    _cscape_relative_taxa_cover(dataset::Dataset)::YAXArray

Calculate relative taxa cover metric for cscape data.

Note: It would be possible to instead write this function as a wrapper of the already s
implemented metric. However the metric computes the metrics for individual scenarios one 
at a time so it would be slower.
"""
function _cscape_relative_taxa_cover(
    dataset::Dataset, location_area::AbstractArray
)::YAXArray
    loc_mask = dataset.properties["location_mask"]
    cube = reformat_cube(dataset.cover, loc_mask)
    dim_sum = (:thermal_tolerance, :intervened)
    cube = dropdims(sum(cube, dims=dim_sum), dims=dim_sum)
    # Reshape to allow vectorised multiplication
    location_area = reshape(location_area, (1, length(cube.sites), 1, 1))
    cube = dropdims(sum(cube .* location_area, dims=:sites), dims=:sites)
    return cube ./ sum(location_area)
end

function _cscape_relative_juvenile(dataset::Dataset)::YAXArray
    loc_mask = dataset.properties["location_mask"]
    cover = reformat_cube(dataset.cover, loc_mask)
    size_counts = reformat_cube(dataset.size_classes, loc_mask)
    size_desc = reformat_cube(dataset.coral_size_diameter)
    return _relative_size_class_cover(cover, size_counts, size_desc, 0.0, 5.0)
end

function _safeFindFirst(condition, iterable)
    ind = findfirst(condition, iterable)
    isnothing(ind) && return length(iterable)
    return ind
end

"""
    _create_juvenile_index_threshold(size_classes::YAXArray, threshold::Float64)::Vector{Int64}
"""
function _create_juvenile_index_threshold(size_classes::YAXArray, threshold::Float64)::Vector{Int64}
    return [_safeFindFirst(x -> x > threshold, size_classes[taxa=At(i)]) for i in size_classes.taxa]
end

function _relative_size_class_cover(
    cover::AbstractArray{<:Real},
    count::AbstractArray{<:Real},
    size_desc::AbstractArray{<:Real},
    lower_bound::Float64,
    upper_bound::Float64
)::YAXArray 
    lb_indices::Vector{Int} = _create_juvenile_index_threshold(size_desc, lower_bound)
    ub_indices::Vector{Int} = _create_juvenile_index_threshold(size_desc, upper_bound)
    count = dropdims(sum(count, dims=(:thermal_tolerance, :intervened)), dims=(:thermal_tolerance, :intervened))
    # convert coral diameter to areas and convert cm^2 to m^2
    coral_areas = (size_desc .* size_desc .* π) ./ 4e-4
    # Force reshape and disk array to load data and not compute reshape lazily
    coral_areas = reshape(
        coral_areas, (1, 1, length(coral_areas.taxa), 1, length(coral_areas.size_bins))
    )[:, :, :, :, :]
    relative_cover = ZeroDataCube(
        ;T=Float64,
        timesteps=collect(count.timesteps),
        sites=collect(count.sites),
        taxa=collect(count.taxa),
        scenarios=collect(count.scenarios)
    )
    count = count .* coral_areas
    for coral_class in relative_cover.taxa
        all_relative_cover = dropdims(sum(
            count[taxa=At(coral_class)], 
            dims=:size_bins
        ), dims=:size_bins)
        all_relative_cover[all_relative_cover .== 0] .= 1.0
        relative_cover[taxa=At(coral_class)] = dropdims(sum(
            count[taxa=At(coral_class), size_bins=At(lb_indices[coral_class]:ub_indices[coral_class])], 
            dims=:size_bins
            ), dims=:size_bins) ./ all_relative_cover
    end
    cover = dropdims(sum(cover, dims=(:thermal_tolerance, :intervened)), dims=(:thermal_tolerance, :intervened))
    # Force load and evaluatiun of reshape to speed up succeeding multiplication
    return dropdims(sum(relative_cover .* cover, dims=:taxa), dims=:taxa) end

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
