using ADRIA: EnvLayer, GDF, ZeroDataCube, DataCube
using ArchGDAL: centroid
using CSV
using DataFrames
using YAXArrays

struct CScapeResultSet <: ResultSet
    name::String
    RCP::String

    loc_ids
    loc_area::Vector{Float64}
    loc_max_coral_cover::Vector{Float64}
    loc_centroids
    env_layer_md::EnvLayer
    connectivity_data
    loc_data
    scenario_groups
    raw_data::Vector{NcFile}

    inputs
    sim_constants
    model_spec::DataFrame

    # raw::AbstractArray
    outcomes
    # Cscape uses different size classes
    coral_size_diameter::YAXArray
end

"""
    load_results(::Type{CScapeResultSet}, data_dir::String)::CScapeResultSet
    load_results(::Type{CscapeResultSet}, data_dir::String, result_dir::String)::CScapeResultSet
    load_results(::Type{CScapeResultSet}, data_dir::String, result_files::Vector{String})::CScapeResultSet

Interface for loading C~scape model outputs.

See the [Loading C~scape Results](@ref) section for details on expected directory structure.

# Examples
```julia
rs = ADRIA.load_results(CScapeResultSet, "a C~scape dataset of interest")
```
"""
function load_results(::Type{CScapeResultSet}, data_dir::String)::CScapeResultSet
    return load_results(CScapeResultSet, data_dir, joinpath(data_dir, "results"))
end
function load_results(
    ::Type{CScapeResultSet}, data_dir::String, result_dir::String
)::CScapeResultSet
    return load_results(CScapeResultSet, data_dir, _get_result_paths(result_dir))
end
function load_results(
    ::Type{CScapeResultSet}, data_dir::String, result_files::Vector{String}
)::CScapeResultSet
    !isdir(data_dir) ? error("Expected a directory but received $(data_dir)") : nothing

    scenario_spec_path::String = joinpath(data_dir, "ScenarioID.csv")
    scenario_spec::DataFrame = DataFrame(CSV.File(scenario_spec_path))

    # 100x faster then YAXArrays
    datasets::Vector{NcFile} = NetCDF.open.(result_files)
    inputs::DataFrame = _recreate_inputs_dataframe(datasets, scenario_spec)
    model_spec::DataFrame = _create_model_spec(CScapeResultSet, inputs)

    # NetCDF auto closes when the reference to raw_set is lost
    raw_set = datasets[1]

    res_name::String = _get_result_name(raw_set)
    res_rcp::String = _get_rcp(raw_set)

    init_cover_path = joinpath(data_dir, "initial_cover", "initial_cover.csv")
    init_data::DataFrame = CSV.read(init_cover_path, DataFrame; header=true)

    !haskey(raw_set.vars, "reef_siteid") ? error("Unable to find location ids.") : nothing
    location_ids = init_data.reef_siteid

    gpkg_path = _get_gpkg_path(data_dir)
    geodata = GDF.read(gpkg_path)

    connectivity_path = joinpath(data_dir, "connectivity/connectivity.csv")
    connectivity = CSV.read(connectivity_path, DataFrame; comment="#", header=true)

    # There is missing location data in site data. Use intersection
    gpkg_mask = BitVector([loc_name in location_ids for loc_name in geodata.reef_siteid])
    conn_mask = BitVector([
        (loc_name in location_ids) && (loc_name in geodata.reef_siteid)
        for loc_name in names(connectivity)
    ])

    geodata = geodata[gpkg_mask, :]
    conn_sites = names(connectivity)[conn_mask]
    connectivity = connectivity[conn_mask[2:end], conn_mask]

    # Re order data to match location ordering
    geo_id_order = [first(findall(x .== geodata.reef_siteid)) for x in location_ids]
    conn_id_order = [first(findall(x .== Vector(conn_sites))) for x in location_ids]

    geodata = geodata[geo_id_order, :]
    connectivity = connectivity[conn_id_order, conn_id_order]

    timeframe = 2007:2099
    if !haskey(raw_set.gatts, "temporal_range")
        @warn "Unable to find timeframe defaulting to $(timeframe[1]):$(timeframe[end])"
    else
        tf_str = split(raw_set.gatts["temporal_range"], ":")
        if tf_str[1] != "Inf"
            timeframe = parse(Int, tf_str[1]):parse(Int, tf_str[2])
        else
            timeframe = raw_set.vars["year"][1]:raw_set.vars["year"][end]
        end
    end

    env_layer_md::EnvLayer = EnvLayer(
        data_dir,
        gpkg_path,
        "site_id",
        "Reef",
        "",
        connectivity_path,
        "",
        "",
        timeframe
    )

    location_max_coral_cover = geodata.k ./ 100
    location_centroids = [centroid(multipoly) for multipoly ∈ geodata.geom]

    outcomes = Dict{Symbol,YAXArray}()
    # add precomputed metrics for comptability
    outcomes[:relative_cover] = _cscape_relative_cover(datasets)

    scen_groups = Dict(
        :counterfactual => BitVector(true for _ in outcomes[:relative_cover].scenarios)
    )

    return CScapeResultSet(
        res_name,
        res_rcp,
        location_ids,
        geodata.area,
        location_max_coral_cover,
        location_centroids,
        env_layer_md,
        connectivity,
        geodata,
        scen_groups,
        datasets,
        inputs,
        SimConstants(),
        model_spec,
        outcomes,
        _construct_coral_sizes(raw_set)
    )
end

"""
    _ncvar_dim_names(ncvar::NcVar)::Vector{Symbol}

Get the dimension names from a netcdf variable
"""
function _ncvar_dim_names(ncvar::NcVar)::Vector{Symbol}
    return [Symbol(d.name) for d in ncvar.dim]
end

"""
    _construct_dim(nc_dim::NcDim, nc_handle::NcFile)::Dim

Construct a dimension from from the netcdf file. If the dimension has no variable array
associated with it, default to standard index range.
"""
function _construct_dim(nc_dim::NcDim, nc_handle::NcFile)::Dim
    dim_name::String = nc_dim.name

    # If there is variable for the dimension return the standard index range, 1:n
    if !haskey(nc_handle.vars, dim_name)
        return Dim{Symbol(dim_name)}(1:Int64(nc_dim.dimlen))
    end

    return Dim{Symbol(dim_name)}(collect(nc_handle[dim_name]))
end

"""
    _construct_axlist(nc_var::NcVar, nc_handle::NcFile)::Tuple

Construct the axis list required for the construction of a YAXrrays.
"""
function _construct_axlist(nc_var::NcVar, nc_handle::NcFile)::Tuple
    return axlist = Tuple([
        _construct_dim(d, nc_handle) for d in nc_var.dim
    ])
end

"""
    _construct_coral_sizes(nc_handle::NcFile)::YAXArray

Construct the coral size diameter YAXArray for the ResultSet.
"""
function _construct_coral_sizes(nc_handle::NcFile)::YAXArray
    coral_size_var = nc_handle["coral_size_diameter"]
    coral_sizes::Matrix{Float64} = NetCDF.readvar(coral_size_var)
    dims::Tuple = _construct_axlist(coral_size_var, nc_handle)
    return YAXArray(dims, coral_sizes, coral_size_var.atts)
end

"""
    _get_scenario_id(datasets::Dataset)::Int

Get the scenario id contained in the meta data of the netcdf. Transform to form expected in
scenario id expected in scenario spec.

1   -> 100001
701 -> 100701
"""
function _get_scenario_id(nc_handle::NcFile)::Int
    scenario_id = parse(Int, nc_handle.gatts["scenario_ID"])
    if scenario_id > 100000
        return scenario_id
    end
    return scenario_id + 100000
end

"""
    _get_functional_types(nc_handle::NcFile)::Vector{String}
"""
function _get_functional_types(nc_handle::NcFile)::Vector{String}
    raw_names::String = nc_handle["coral_size_diameter"].atts["column_names"]
    return string.(split(raw_names, " "))
end

function _default_missing(value, default)
    return ismissing(value) ? default : value
end

"""
    _recreate_inputs_dataframe(nc_handles::Vector{NcFile}, scenario_spec::DataFrame)::DataFrame

Construct the inputs datadrame from dataset properties and scenario table.
"""
function _recreate_inputs_dataframe(
    nc_handles::Vector{NcFile}, scenario_spec::DataFrame
)::DataFrame
    # Get rows from scenario spec dataframe corresponding to the scenarios
    scenario_idxs::Vector{Int} = [
        findfirst(x -> x == idx, scenario_spec.ID) for idx in _get_scenario_id.(nc_handles)
    ]
    scenario_rows::Vector{DataFrameRow} = [scenario_spec[idx, :] for idx in scenario_idxs]

    # Convert climate scenarios to factors
    rcps::Vector{Float64} =
        parse.(Float64, [nc_handle.gatts["ssp"][(end - 1):end] for nc_handle in nc_handles])

    # Convert climate models to factors
    input_scenarios::Vector{Tuple{String,Int}} = [
        (
            nc_handle.gatts["climate model"],
            nc_handle.gatts["climate_model_realisation_point"]
        ) for nc_handle in nc_handles
    ]
    input_inds::Vector{Int} = [
        findfirst(x -> x == scen, unique(input_scenarios)) for scen in input_scenarios
    ]

    fragmented_spec::Vector{DataFrame} =
        _create_inputs_dataframe.(
            nc_handles, scenario_rows, rcps, input_inds
        )
    scenarios::DataFrame = reduce(vcat, fragmented_spec)
    return scenarios
end
"""
    _create_inputs_dataframe(nc_handle::NcFile, scenario_spec::DataFrameRow, rcp::Float64, input_index::Int)::DataFrame

Construct inputs dataframe for a singular netcdf file.
"""
function _create_inputs_dataframe(
    nc_handle::NcFile,
    scenario_spec::DataFrameRow,
    rcp::Float64,
    input_index::Int
)::DataFrame
    functional_types::Vector{String} = _get_functional_types(nc_handle)
    n_draws::Int =
        "draws" in keys(nc_handle.dim) ? length(get(nc_handle.dim, "draws", [1])) : 1

    dhws::Vector{Float64} = Float64.(repeat([input_index], n_draws))
    cyclones::Vector{Float64} = Float64.(repeat([input_index], n_draws))

    # Thermal tolerance bins are stored as lb_interval_ub
    lb_t, int_t1, int_t2, ub_t =
        parse.(
            Float64, split(scenario_spec.HeatToleranceGroups, "_")
        )
    init_heat_tol_mean, init_heat_tol_std =
        parse.(
            Float64, split(scenario_spec.HeatToleranceInit, "_")
        )
    heritability1, heritability2 = parse.(
        Float64, split(scenario_spec.Heritability, "_")
    )
    plasticity::Int64 = scenario_spec.Plasticity

    settle_probability = parse.(Float64, split(scenario_spec.settle_prob, "_"))

    settle_probs_kwargs = Dict(
        Symbol(ft * "_settle_probability") => prob
        for (ft, prob) in zip(functional_types, settle_probability)
    )

    # Intervention factors
    deployment_area = _default_missing(scenario_spec[Symbol("Deployment area")], 0.0)
    total_corals = _default_missing(scenario_spec.TotalCorals, 1.0)

    n_seeded = [0.0, 0.0, 0.0, 0.0, 0.0]
    taxa_deployed = if ismissing(scenario_spec.species)
        []
    else
        parse.(
            Int64, split(scenario_spec.species, '_')
        )
    end
    n_seeded[taxa_deployed] .= total_corals / length(taxa_deployed)
    corals_deployed = Dict(
        Symbol(ft * "_n_seeded") => n_corals
        for (ft, n_corals) in zip(functional_types, n_seeded)
    )

    enhancement_mean, enhancement_std =
        parse.(
            Float64, split(_default_missing(scenario_spec.Enhancement, "0_0"), '_')
        )

    intervention_start =
        _default_missing(
            scenario_spec.InterventionYears_start, nc_handle["year"][1]
        ) - nc_handle["year"][1]

    duration = _default_missing(scenario_spec.duration, 0.0)
    frequency = _default_missing(scenario_spec.frequency, 0.0)
    n_dep_locations = count(x -> x == '/', _default_missing(scenario_spec.Reef_siteids, ""))

    coral_dict = merge(settle_probs_kwargs, corals_deployed)

    return DataFrame(;
        dhw_scenario=dhws,
        cyc_scenario=cyclones,
        RCP=repeat([rcp], n_draws),
        thermal_tol_lb=repeat([lb_t], n_draws),
        thermal_tol_int1=repeat([int_t1], n_draws),
        thermal_tol_int2=repeat([int_t2], n_draws),
        thermal_tol_ub=repeat([ub_t], n_draws),
        init_heat_tol_mean=repeat([init_heat_tol_mean], n_draws),
        init_heat_tol_std=repeat([init_heat_tol_std], n_draws),
        heritability1=repeat([heritability1], n_draws),
        heritability2=repeat([heritability2], n_draws),
        plasticity=repeat([plasticity], n_draws),
        intervention_start=repeat([intervention_start], n_draws),
        intervention_duration=repeat([duration], n_draws),
        intervention_frequency=repeat([frequency], n_draws),
        deployment_area=repeat([deployment_area], n_draws),
        n_deployment_locations=repeat([n_dep_locations], n_draws),
        enhancement_mean=repeat([enhancement_mean], n_draws),
        enhancement_std=repeat([enhancement_std], n_draws),
        coral_dict...
    )
end

"""
    _create_model_spec(::Type{CScapeResultSet}, scenario_spec::DataFrame)::DataFrame

Create partial model specification from scenario specification.
"""
function _create_model_spec(::Type{CScapeResultSet}, scenario_spec::DataFrame)::DataFrame
    # Retrieve coral factors
    factor_names::Vector{String} = names(scenario_spec)
    # Settler Probability
    settle_names = filter(factor -> contains(factor, "settle_probability"), factor_names)
    settle_readable = human_readable_name.(settle_names)
    settle_names = Symbol.(settle_names)
    settle_ptype = repeat(["continuous"], length(settle_names))
    settle_lb = repeat([0.0], length(settle_names))
    settle_ub = repeat([1.0], length(settle_names))

    # Number of corals seeded
    seeded_names = filter(factor -> contains(factor, "n_seeded"), factor_names)
    seeded_readable = human_readable_name.(seeded_names)
    seeded_names = Symbol.(seeded_names)
    seeded_ptype = repeat(["continuous"], length(settle_names))
    seeded_lb = repeat([0.0], length(seeded_names))
    seeded_ub = repeat([5e7], length(seeded_names))

    # Construct default model spec
    fieldname::Vector{Symbol} = [
        :dhw_scenario, # Environment Layer
        :cyc_scenario,
        :thermal_tol_lb, # Corals
        :thermal_tol_int1,
        :thermal_tol_int2,
        :thermal_tol_ub,
        :init_heat_tol_mean,
        :init_heat_tol_std,
        :heritability1,
        :heritability2,
        :plasticity,
        :intervention_start, # Interventions
        :intervention_duration,
        :intervention_frequency,
        :deployment_area,
        :n_deployment_locations,
        :enhancement_mean,
        :enhancement_std,
        settle_names..., # Corals
        seeded_names... # Intervention
    ]
    descriptions::Vector{String} = [
        "DHW Scenario",
        "Cyclone Scenario",
        "Natural thermal tolerance lower bound",
        "Natural thermal tolerance int1",
        "Natural thermal tolerance int2",
        "Natural thermal tolerance upper bound",
        "Initial natural thermal tolerance mean",
        "Initial natural thermal tolerance standard deviation",
        "Heritability mean",
        "Heritability standard deviation",
        "Plasticity",
        "First Intervention Year",
        "Duration of interventions (years)",
        "Frequency of interventions after first year (years)",
        "Area of coral deployments (m^2)",
        "Number of locations intervened on",
        "Artificial thermal tolerance enhancement mean",
        "Artificial thermal tolerance enhancement standard deviation",
        settle_readable...,
        seeded_readable...
    ]
    human_names::Vector{String} = [
        "DHW scenario",
        "Cyclone scenario",
        "Thermal Tolerance Lower Bound",
        "Thermal Tolerance int1",
        "Thermal Tolerance int2",
        "Thermal Tolerance Upper Bound",
        "Initial Thermal Tolerance Mean",
        "Initial Thermal Tolerance Standard deviation",
        "Heritability Mean",
        "Heritability Standard Deviation",
        "Plasticity",
        "First Intervention Year",
        "Duration of Interventions",
        "Frequency of Interventions",
        "Area of Deployments",
        "Number of Deployment Locations",
        "Enhancement Mean",
        "Enhancement Standard Deviation",
        settle_readable...,
        seeded_readable...
    ]
    ptype::Vector{String} = [
        "unordered categorical",
        "unordered categorical",
        "continuous",
        "continuous",
        "continuous",
        "continuous",
        "continuous",
        "continuous",
        "continuous",
        "continuous",
        "continuous",
        "ordered categorical",
        "ordered categorical",
        "ordered categorical",
        "continuous",
        "ordered discrete",
        "continuous",
        "continuous",
        settle_ptype...,
        seeded_ptype...
    ]
    lower_bound::Vector{Float64} = [
        1.0,
        1.0,
        -6.0,
        0.0,
        0.0,
        8.0,
        -1.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        settle_lb...,
        seeded_lb...
    ]
    upper_bound::Vector{Float64} = [
        50.0,
        50.0,
        0.0,
        40.0,
        24.0,
        30.0,
        1.0,
        3.0,
        1.0,
        0.5,
        1.0,
        25.0,
        10.0,
        10.0,
        maximum(scenario_spec.deployment_area),
        100.0,
        30.0,
        2.0,
        settle_ub...,
        seeded_ub...
    ]
    dist_params::Vector{String} = [
        string((lb, ub)) for (lb, ub) in zip(lower_bound, upper_bound)
    ]
    component::Vector{String} = vcat(
        repeat(["EnvironmentalLayer"], 2),
        repeat(["Coral"], 9),
        repeat(["Intervention"], 7),
        repeat(["Coral"], length(settle_names)),
        repeat(["Intervention"], length(seeded_names))
    )
    return DataFrame(;
        component=component,
        fieldname=fieldname,
        description=descriptions,
        name=human_names,
        ptype=ptype,
        dist_params=dist_params,
        lower_bound=lower_bound,
        upper_bound=upper_bound
    )
end

"""
    _get_result_paths(result_dir::String)::Vector{String}

Get the names of all result netcdf result files in the given directory. Assumes filenames
match "NetCDF_Scn.*.nc.
"""
function _get_result_paths(result_dir::String)::Vector{String}
    possible_files = filter(isfile, readdir(result_dir; join=true))
    possible_files = filter(x -> occursin(r"NetCDF_Scn.*.nc", x), possible_files)
    if length(possible_files) == 0
        error("Unable to find result netcdf file in subdirectory $(res_subdir)")
    end
    return possible_files
end

"""
    _get_gpkg_path(data_dir::String)

Get the path to the gpkg file contained in the site_data subdirectory.
"""
function _get_gpkg_path(data_dir::String)
    gpkg_dir = joinpath(data_dir, "site_data")
    possible_files = filter(isfile, readdir(gpkg_dir; join=true))
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
function _get_result_name(ds::NcFile)::String
    name = "CScape Results"
    if !haskey(ds.gatts, "title")
        msg = "Unable to find key `title` in dataset properties, "
        msg *= "defaulting to `CScape Results`"
        @warn msg
    else
        name = ds.gatts["title"]
    end
    return name
end

"""
    _get_rcp(ds::NcFile)::String

Get the RCP from the ssp field in dataset properties. Assume RCP is last two numbers.
"""
function _get_rcp(ds::NcFile)::String
    rcp = ""
    if !haskey(ds.gatts, "ssp")
        @warn "Unable to find key `ssp` in dataset properties."
    else
        rcp = ds.gatts["ssp"][(end - 1):end]
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
        cscape_cube = renameaxis!(cscape_cube, cscape_nm => adria_nm)
        final_ordering[current_index] = findfirst(x -> x == cscape_nm, dim_names)
        current_index += 1
    end
    # append the indices of non-adria axis to end of permutation vector
    for dim_name in dim_names
        if dim_name in cscape_names
            continue
        end
        final_ordering[current_index] = findfirst(x -> x == dim_name, dim_names)
        current_index += 1
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

function _throw_missing_variable(nc_handle::NcFile, var_name::String)::Nothing
    msg = "NetCDF $(nc_handle.gatts["scenario_ID"]) does not "
    msg *= "contain $(String.(var_name)) variable"
    throw(ArgumentError(msg))

    return nothing
end

"""
    _drop_sum(cube::YAXArray, red_dims::Tuple)::YAXArray

Sum over given dimensions and drop the same given dimensions.
"""
function _drop_sum(cube::YAXArray, red_dims::Tuple)::YAXArray
    return dropdims(sum(cube; dims=red_dims); dims=red_dims)
end
function _drop_sum(cube::AbstractArray, red_dims::Tuple)::AbstractArray
    return dropdims(sum(cube; dims=red_dims); dims=red_dims)
end

"""
    _cscape_relative_cover(nc_handle::NcFile)::YAXArray
    _cscape_relative_cover(nc_handles::Vector{NcFiles})::YAXArray

Calculate relative cover metric for C~scape data.
"""
function _cscape_relative_cover(nc_handle::NcFile)::Array
    # Throw error and identify specific misformatted file.
    if !haskey(nc_handle.vars, "cover")
        _throw_missing_variable(nc_handle, :cover)
    end
    if !haskey(nc_handle.vars, "area")
        _throw_missing_variable(nc_handle, "area")
    end
    if !haskey(nc_handle.vars, "k")
        _throw_missing_variable(nc_handle, "k")
    end

    # thermal tolerance and functional group dimensions
    dim_sum = (3, 4)

    multi_scenario::Bool = :draws in keys(nc_handle.dim)

    n_scens = multi_scenario ? Int64(nc_handle.dim["draws"].dimlen) : 0
    n_tsteps::Int64 = Int64(nc_handle.dim["year"].dimlen)
    n_locs::Int64 = Int64(nc_handle.dim["reef_sites"].dimlen)

    if multi_scenario
        relative_cover = zeros(n_scens, n_tsteps, n_locs)
        reshape_tuple = (1, 1, n_locs)
    else
        relative_cover = zeros(n_tsteps, n_locs)
        reshape_tuple = (1, n_locs)
    end

    # Calculate relative cover from non-intervened areas
    # Force load with dimensions [intervened ⋅ reef_sites]

    # We want to read in data before any calculations are conducted to maintain
    # speed/performance. Using `read()` converts data into a base Matrix, so we
    # use a dummy selector to retain the nice YAXArray data type while also reading
    # the data into memory.
    n_dims = nc_handle["cover"].ndim

    area::Matrix{Float64} = NetCDF.readvar(nc_handle["area"])
    habitable_area::Vector{Float64} = NetCDF.readvar(nc_handle["k"])
    cover = NetCDF.readvar(nc_handle["cover"])
    relative_cover .=
        _drop_sum(
            cover[:, :, 1, :, :], dim_sum
        ) .* reshape(
            area[1, :],
            reshape_tuple
        )

    # Only calculate if there are area values > 0
    if !any(area[2, :] .> 0.0)
        relative_cover .+=
            _drop_sum(
                cover[:, :, 2, :, :], dim_sum
            ) .* reshape(
                area[2, :],
                reshape_tuple
            )
    end

    relative_cover ./=
        reshape(
            sum(area; dims=1) .* habitable_area' ./ 100, reshape_tuple
        ) .* 100

    # ADRIA assumes a shape of [timesteps ⋅ locations ⋅ scenarios]
    if multi_scenario
        return permutedims(relative_cover, (2, 3, 1))
    end

    return relative_cover
end
function _cscape_relative_cover(nc_handles::Vector{NcFile})::YAXArray
    # Assume same number of locations and timesteps for every dataset
    n_locs::Int64 = nc_handles[1].dim["reef_sites"].dimlen

    # Determine the number of entries for each file
    n_scens::Vector{Int64} = _n_scenarios.(nc_handles)

    relative_cover = ZeroDataCube(;
        T=Float64,
        timesteps=Vector(NetCDF.readvar(nc_handles[1]["year"])),
        sites=1:n_locs,
        scenarios=1:sum(n_scens)
    )

    # Could implement a simpler loading strategy in the case where all datasets hold a
    # single scenario. This would simplify loading...
    # all(diff(cumsum(n_scens)) .== 1)

    start_end_pos = _data_position(n_scens)
    @showprogress desc = "Loading datasets" for (ds_idx, nc_handle) in
                                                zip(start_end_pos, nc_handles)
        relative_cover[scenarios=ds_idx[1]:ds_idx[2]] = _cscape_relative_cover(nc_handle)
    end

    return relative_cover
end

"""
    _n_scenarios(nc_handle::NcFile)::Int64

Get the number of scenario or draws in a dataset.
"""
function _n_scenarios(nc_handle::NcFile)::Int64
    if "draws" in keys(nc_handle.dim)
        return nc_handle.dim["draws"].dimlen
    end

    return 1
end

"""
    _data_position(n_per_dataset::Vector{Int64})

Determine the index positions of a set of data if they were collated into a single dataset.

# Examples

```julia
# Here we have five datasets with hetrogenous number of data
n_scens = [1, 10, 5, 1, 12]
start_pos, end_pos = _data_position(n_scens)
# 5-element Vector{Tuple{Int64, Int64}}:
#  (1, 1)
#  (2, 11)
#  (12, 16)
#  (17, 17)
#  (18, 29)
```

# Arguments
- `n_per_dataset` : Indicates the number of scenarios per entry

# Returns
The start and end position of each entry
"""
function _data_position(n_per_dataset::Vector{Int64})::Vector{Tuple{Int64,Int64}}
    n = length(n_per_dataset)
    starts = Vector{Int64}(undef, n)
    ends = Vector{Int64}(undef, n)

    starts[1] = 1
    ends[1] = n_per_dataset[1]

    for i in 2:n
        starts[i] = ends[i - 1] + 1
        ends[i] = starts[i] + n_per_dataset[i] - 1
    end

    return collect(zip(starts, ends))
end

"""
    load_variable!(rs::CScapeResultSet, variable_name::Symbol; aggregate_ttol=sum)::YAXArray

Load given variable into a single YAXArray from all datasets in the result store.
If thermal tolerance is present remove the dimension using the aggregation function. Convert
dimensions and dimension ordering to standard ADRIA ordering.

Return the variable and write the variable to the outcomes dictionary.
"""
function load_variable!(
    rs::CScapeResultSet,
    variable_name::Symbol;
    aggregate_ttol=sum
)::YAXArray
    var_name::String = String(variable_name)
    if !haskey(rs.raw_data[1].vars, var_name)
        throw(ArgumentError("Unable to find variable $(var_name) in datasets"))
    end
    n_scens::Vector{Int64} = _n_scenarios.(rs.raw_data)

    dims = Tuple(
        _construct_dim(d, rs.raw_data[1]) for d in rs.raw_data[1].vars[var_name].dim
        if !(d.name in ["thermal_tolerance", "draws"])
    )
    dims = (
        Dim{:draws}(1:sum(n_scens)),
        dims...
    )
    shape = Tuple(length(d) for d in dims)
    output_variable = YAXArray(
        dims,
        zeros(eltype(rs.raw_data[1][var_name]), shape...)
    )

    # Only aggregate the thermal tolerance if the dimension is present
    agg_func = x -> NetCDF.readvar(x[var_name])
    if :thermal_tolerance in _ncvar_dim_names(rs.raw_data[1].vars[var_name])
        thermal_tol_dim_idx::Int64 = findfirst(
            x -> x == :thermal_tolerance,
            _ncvar_dim_names(rs.raw_data[1].vars[var_name])
        )
        agg_func =
            x -> dropdims(
                aggregate_ttol(
                    NetCDF.readvar(x[var_name]); dims=thermal_tol_dim_idx
                ); dims=thermal_tol_dim_idx)
    end

    cur_indx = 1
    @showprogress for (n_sc, nc_handle) in zip(n_scens, rs.raw_data)
        if n_sc == 1
            output_variable[draws=cur_indx] .= agg_func(nc_handle)
        else
            output_variable[draws=cur_indx:(cur_indx + n_sc - 1)] .= agg_func(nc_handle)
        end
        cur_indx += n_sc
    end
    # Convert cscape yaxarray to standard ADRIA formats
    output_variable = reformat_cube(output_variable)
    rs.outcomes[variable_name] = output_variable

    return output_variable
end

function Base.show(io::IO, mime::MIME"text/plain", rs::CScapeResultSet)
    rcps = rs.RCP
    scens = length(rs.outcomes[:relative_cover].scenarios)
    sites = length(rs.outcomes[:relative_cover].sites)
    tf = rs.env_layer_md.timeframe

    return println("""
        Name: $(rs.name)

        Results stored at: $(rs.env_layer_md.dpkg_path)

        RCP(s) represented: $(rcps)
        Scenarios run: $(scens)
        Number of sites: $(sites)
        Timesteps: $(tf)
    """)
end
