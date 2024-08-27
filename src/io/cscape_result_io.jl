using Core: Argument
using ADRIA: EnvLayer, GDF, ZeroDataCube, DataCube
using ArchGDAL: centroid
using CSV
using DataFrames
using YAXArrays

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
    raw_data::Vector{Dataset}

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

Interface for loading CScape model outputs.
"""
function load_results(::Type{CScapeResultSet}, data_dir::String)::CScapeResultSet
    return load_results(CScapeResultSet, data_dir, joinpath(data_dir, "results"))
end
function load_results(::Type{CScapeResultSet}, data_dir::String, result_dir::String)::CScapeResultSet
    return load_results(CScapeResultSet, data_dir, _get_result_paths(result_dir))
end
function load_results(::Type{CScapeResultSet}, data_dir::String, result_files::Vector{String})::CScapeResultSet
    !isdir(data_dir) ? error("Expected a directory but received $(data_dir)") : nothing

    scenario_spec_path::String = joinpath(data_dir, "ScenarioID.csv")
    scenario_spec::DataFrame = DataFrame(CSV.File(scenario_spec_path))

    datasets::Vector{Dataset} = open_dataset.(result_files)
    inputs::DataFrame = _recreate_inputs_dataframe(datasets, scenario_spec)
    model_spec::DataFrame = _create_model_spec(CScapeResultSet, inputs)

    # Assume all result set have the same locations
    raw_set = datasets[1]

    res_name::String = _get_result_name(raw_set)
    res_rcp::String = _get_rcp(raw_set)

    init_cover_path = joinpath(data_dir, "initial_cover", "initial_cover.csv")
    init_data::DataFrame = CSV.read(init_cover_path, DataFrame, header=true)

    !haskey(raw_set.cubes, :reef_siteid) ? error("Unable to find location ids.") : nothing
    location_ids = init_data.reef_siteid

    gpkg_path = _get_gpkg_path(data_dir)
    geodata = GDF.read(gpkg_path)

    geodata = _manual_site_additions(geodata, location_ids, raw_set)

    connectivity_path = joinpath(data_dir, "connectivity/connectivity.csv")
    connectivity = CSV.read(connectivity_path, DataFrame, comment="#", header=true)

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
    if !haskey(raw_set.properties, "temporal_range")
        @warn "Unable to find timeframe defaulting to $(timeframe[1]):$(timeframe[end])"
    else
        tf_str = split(raw_set.properties["temporal_range"], ":")
        if tf_str[1] != "Inf"
            timeframe = parse(Int, tf_str[1]):parse(Int, tf_str[2])
        else
            timeframe = raw_set.year[1]:raw_set.year[end]
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

    location_max_coral_cover = 1 .- geodata.k ./ 100
    location_centroids = [centroid(multipoly) for multipoly ∈ geodata.geom]

    outcomes = Dict{Symbol, YAXArray}()
    # add precomputed metrics for comptability
    outcomes[:relative_cover] = _cscape_relative_cover(datasets)
    # outcomes[:relative_taxa_cover] = _cscape_relative_taxa_cover(raw_set, geodata.area)

    scen_groups = Dict(:counterfactual=>BitVector(true for _ in outcomes[:relative_cover].scenarios))

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
        reformat_cube(raw_set.coral_size_diameter)
    )
end

"""
    _manual_site_additions(geodata::DataFrame, loc_ids, dataset::Dataset)::DataFrame

Add missing sites to geopackage
"""
function _manual_site_additions(geodata::DataFrame, loc_ids, dataset::Dataset)::DataFrame
    row_indx::Int64 = findfirst(x->x=="Moore_MR_S_39", geodata.reef_siteid)
    row_cpy = copy(geodata[row_indx, :])
    push!(geodata, row_cpy)
    data_indx = findfirst(x->x=="Moore_MR_S_40", loc_ids)
    geodata[end, :site_id] = "MR_S_40"
    geodata[end, :reef_siteid] = "Moore_MR_S_40"
    geodata[end, :area] = sum(dataset.area[reef_sites=data_indx], dims=:intervened)[1]
    geodata[end, :k] = dataset.k[reef_sites=data_indx][1]

    row_indx = findfirst(x->x=="Milln_MR_OF_3", geodata.reef_siteid)
    row_cpy = copy(geodata[row_indx, :])
    push!(geodata, row_cpy)
    data_indx = findfirst(x->x=="Milln_MR_OF_4", loc_ids)
    geodata[end, :site_id] = "MR_OF_4"
    geodata[end, :reef_siteid] = "Milln_MR_OF_4"
    geodata[end, :area] = sum(dataset.area[reef_sites=data_indx], dims=:intervened)[1]
    geodata[end, :k] = dataset.k[reef_sites=data_indx][1]

    row_indx = findfirst(x->x=="Elford_ER_S_50", geodata.reef_siteid)
    row_cpy = copy(geodata[row_indx, :])
    push!(geodata, row_cpy)
    data_indx = findfirst(x->x=="Elford_ER_S_51", loc_ids)
    geodata[end, :site_id] = "ER_S_51"
    geodata[end, :reef_siteid] = "Elford_ER_S_51"
    geodata[end, :area] = sum(dataset.area[reef_sites=data_indx], dims=:intervened)[1]
    geodata[end, :k] = dataset.k[reef_sites=data_indx][1]

    row_indx = findfirst(x->x=="Elford_ER_S_50", geodata.reef_siteid)
    row_cpy = copy(geodata[row_indx, :])
    push!(geodata, row_cpy)
    data_indx = findfirst(x->x=="Elford_ER_S_52", loc_ids)
    geodata[end, :site_id] = "ER_S_52"
    geodata[end, :reef_siteid] = "Elford_ER_S_52"
    geodata[end, :area] = sum(dataset.area[reef_sites=data_indx], dims=:intervened)[1]
    geodata[end, :k] = dataset.k[reef_sites=data_indx][1]
    return geodata
end

"""
    _get_scenario_id(datasets::Dataset)::Int

Get the scenario id contained in the meta data of the netcdf. Transform to form expected in
scenario id expected in scenario spec.

1   -> 100001
701 -> 100701
"""
function _get_scenario_id(dataset::Dataset)::Int
    scenario_id = parse(Int, dataset.properties["scenario_ID"])
    if scenario_id > 100000
        return scenario_id
    end
    return scenario_id + 100000
end

"""
    _get_functional_types(dataset::Dataset)::Vector{String}
"""
function _get_functional_types(dataset::Dataset)::Vector{String}
    if !haskey(dataset.coral_size_diameter.properties, "column_names")
        @warn "Unable to find names of functional types. Skipping."
        return Vector{String}(undef, 0)
    end
    raw_names::String = dataset.coral_size_diameter.properties["column_names"]
    return string.(split(raw_names, " "))
end

function _default_missing(value, default)
    return ismissing(value) ? default : value
end

"""
    _recreate_inputs_dataframe(datasets::Vector{Dataset}, scenario_spec::DataFrame)::DataFrame

Construct the inputs datadrame from dataset properties and scenario table.
"""
function _recreate_inputs_dataframe(
    datasets::Vector{Dataset}, scenario_spec::DataFrame
)::DataFrame
    # Get rows from scenario spec dataframe corresponding to the scenarios
    scenario_idxs::Vector{Int} = [
        findfirst(x->x==idx, scenario_spec.ID) for idx in _get_scenario_id.(datasets)
    ]
    scenario_rows::Vector{DataFrameRow} = [scenario_spec[idx, :] for idx in scenario_idxs]

    # Convert climate scenarios to factors
    rcps::Vector{Float64} =
        parse.(Float64, [dataset.properties["ssp"][end-1:end] for dataset in datasets])

    # Convert climate models to factors
    input_scenarios::Vector{Tuple{String, Int}} = [(
        dataset.properties["climate model"], dataset.properties["climate_model_realisation_point"]
    ) for dataset in datasets]
    input_inds::Vector{Int} =
        [findfirst(x -> x == scen, unique(input_scenarios)) for scen in input_scenarios]

    fragmented_spec::Vector{DataFrame} = _create_inputs_dataframe.(
        datasets, scenario_rows, rcps, input_inds
    )
    scenarios::DataFrame = reduce(vcat, fragmented_spec)
    return scenarios
end
"""
    _create_inputs_dataframe(dataset::Dataset, scenario_spec::DataFrameRow, rcp::Float64, input_index::Int)::DataFrame

Construct inputs dataframe for a singular netcdf file.
"""
function _create_inputs_dataframe(
    dataset::Dataset,
    scenario_spec::DataFrameRow,
    rcp::Float64,
    input_index::Int
)::DataFrame
    functional_types::Vector{String} = _get_functional_types(dataset)
    n_draws::Int = :draws in keys(dataset.axes) ? length(get(dataset.axes, :draws, [1])) : 1

    dhws::Vector{Float64} = Float64.(repeat([input_index], n_draws))
    cyclones::Vector{Float64} = Float64.(repeat([input_index], n_draws))

    # Thermal tolerance bins are stored as lb_interval_ub
    lb_t, int_t1, int_t2, ub_t = parse.(
        Float64, split(scenario_spec.HeatToleranceGroups, "_")
    )
    init_heat_tol_mean, init_heat_tol_std = parse.(
        Float64, split(scenario_spec.HeatToleranceInit, "_")
    )
    heritability1, heritability2 = parse.(
        Float64, split(scenario_spec.Heritability, "_")
    )
    plasticity::Int64 = scenario_spec.Plasticity

    settle_probability = parse.(Float64, split(scenario_spec.settle_prob, "_"))

    settle_probs_kwargs = Dict(
        Symbol(ft*"_settle_probability") => prob
        for (ft, prob) in zip(functional_types, settle_probability)
    )

    # Intervention factors
    deployment_area = _default_missing(scenario_spec[Symbol("Deployment area")], 0.0)
    total_corals = _default_missing(scenario_spec.TotalCorals, 1.0)

    n_seeded = [0.0, 0.0, 0.0, 0.0, 0.0]
    taxa_deployed = ismissing(scenario_spec.species) ? [] : parse.(
        Int64, split(scenario_spec.species, '_')
    )
    n_seeded[taxa_deployed] .= total_corals / length(taxa_deployed)
    corals_deployed = Dict(
        Symbol(ft*"_n_seeded") => n_corals
        for (ft, n_corals) in zip(functional_types, n_seeded)
    )

    enhancement_mean, enhancement_std = parse.(
        Float64, split(_default_missing(scenario_spec.Enhancement, "0_0"), '_')
    )

    intervention_start = _default_missing(
        scenario_spec.InterventionYears_start, dataset.year[1]
    ) - dataset.year[1]

    duration = _default_missing(scenario_spec.duration, 0.0)
    frequency = _default_missing(scenario_spec.frequency, 0.0)
    n_dep_locations = count(x->x=='/', _default_missing(scenario_spec.Reef_siteids, ""))

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
    return DataFrame(
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
    possible_files = filter(isfile, readdir(result_dir, join=true))
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
        msg = "Unable to find key `title` in dataset properties, "
        msg *= "defaulting to `CScape Results`"
        @warn msg
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

function _throw_missing_variable(dataset::Dataset, var_name::Symbol)::Nothing
    msg = "NetCDF $(dataset.properties["scenario_ID"]) does not "
    msg *= "contain $(String.(var_name)) variable"
    throw(ArgumentError(msg))

    return nothing
end

"""
    _drop_sum(cube::YAXArray, red_dims)::YAXArray

Sum over given dimensions and drop the same given dimensions.
"""
function _drop_sum(cube::YAXArray, red_dims)::YAXArray
    return dropdims(sum(cube, dims=red_dims), dims=red_dims)
end

"""
    _cscape_relative_cover(dataset::Dataset)::YAXArray
    _cscape_relative_cover(datasets::Vector{Dataset})::YAXArray

Calculate relative cover metric for C~scape data.
"""
function _cscape_relative_cover(dataset::Dataset)::Array
    # Throw error and identify specific misformatted file.
    if !haskey(dataset.cubes, :cover)
        _throw_missing_variable(dataset, :cover)
    end
    if !haskey(dataset.cubes, :area)
        _throw_missing_variable(dataset, :area)
    end
    if !haskey(dataset.cubes, :k)
        _throw_missing_variable(dataset, :k)
    end

    dim_sum = (:ft, :thermal_tolerance)

    multi_scenario::Bool = :draws in keys(dataset.axes)

    n_scens = multi_scenario ? length(dataset.draws) : 0
    n_tsteps::Int64 = length(dataset.year)
    n_locs::Int64 = length(dataset.reef_sites)

    if multi_scenario
        relative_cover = zeros(n_scens, n_tsteps, n_locs)
        reshape_tuple = (1, 1, n_locs)
    else
        relative_cover = zeros(n_tsteps, n_locs)
        reshape_tuple = (1, n_locs)
    end

    # Calculate relative cover from non-intervened areas
    # Force load with dimensions [intervened ⋅ reef_sites]
    area = dataset.area.data[:, :]
    relative_cover .= _drop_sum(
        dataset.cover[intervened=1], dim_sum
    ) .* reshape(
        area[1, :] .* dataset.k[:],
        reshape_tuple
    ) ./ 100.0

    relative_cover .+= _drop_sum(
        dataset.cover[intervened=2], dim_sum
    ) .* reshape(
        area[2, :] .* dataset.k[:],
        reshape_tuple
    ) ./ 100.0

    # ADRIA assumes a shape of [timesteps ⋅ locations ⋅ scenarios]
    if multi_scenario
        return permutedims(relative_cover, (2, 3, 1))
    end

    return relative_cover
end
function _cscape_relative_cover(datasets::Vector{Dataset})::YAXArray
    # Assume same number of locations and timesteps for every dataset
    n_locs::Int64 = length(datasets[1].reef_sites)

    # Determine the number of entries for each file
    n_scens::Vector{Int64} = _n_scenarios.(datasets)

    relative_cover = ZeroDataCube(
        T=Float64,
        timesteps=datasets[1].year.val,
        sites=1:n_locs,
        scenarios=1:sum(n_scens)
    )

    # Could implement a simpler loading strategy in the case where all datasets hold a
    # single scenario. This would simplify loading...
    # all(diff(cumsum(n_scens)) .== 1)

    start_end_pos = _data_position(n_scens)
    for (ds_idx, dataset) in zip(start_end_pos, datasets)
        relative_cover[scenarios=ds_idx[1]:ds_idx[2]] = _cscape_relative_cover(dataset)
    end

    return relative_cover
end

"""
    _n_scenarios(dataset::Dataset)::Int64

Get the number of scenario or draws in a dataset.
"""
function _n_scenarios(dataset::Dataset)::Int64
    if :draws in keys(dataset.axes)
        return length(dataset.draws)
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
function _data_position(n_per_dataset::Vector{Int64})::Vector{Tuple{Int64, Int64}}
    n = length(n_per_dataset)
    starts = Vector{Int64}(undef, n)
    ends = Vector{Int64}(undef, n)

    starts[1] = 1
    ends[1] = n_per_dataset[1]

    for i in 2:n
        starts[i] = ends[i-1] + 1
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
    if !haskey(rs.raw_data[1].cubes, variable_name)
        throw(ArgumentError("Unable to find variable $(String(variable_name)) in datasets"))
    end
    n_scens::Vector{Int64} = _n_scenarios.(rs.raw_data)

    dims = Tuple(
        d for d in rs.raw_data[1].cubes[variable_name].axes
        if !(name(d) in [:thermal_tolerance, :draws])
    )
    dims = (
        Dim{:draws}(1:sum(n_scens)),
        dims...
    )
    shape = Tuple(length(d) for d in dims)
    output_variable = YAXArray(
        dims,
        zeros(eltype(rs.raw_data[1].cubes[variable_name]), shape...)
    )

    # Only aggregate the thermal tolerance if the dimension is present
    agg_func = x -> x
    if :thermal_tolerance in name.(rs.raw_data[1].cubes[variable_name].axes)
        agg_func = x -> dropdims(aggregate_ttol(
            x.cubes[variable_name], dims=:thermal_tolerance
        ), dims=:thermal_tolerance)
    end

    cur_indx = 1
    @showprogress for (n_sc, dataset) in zip(n_scens, rs.raw_data)
        if n_sc == 1
            output_variable[draws=cur_indx] = agg_func(dataset)
        else
            output_variable[draws=cur_indx:cur_indx+n_sc-1] = agg_func(dataset)
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

    println("""
    Name: $(rs.name)

    Results stored at: $(rs.env_layer_md.dpkg_path)

    RCP(s) represented: $(rcps)
    Scenarios run: $(scens)
    Number of sites: $(sites)
    Timesteps: $(tf)
    """)
end
