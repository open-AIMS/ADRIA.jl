using ADRIA: EnvLayer, GDF, ZeroDataCube, DataCube
using ArchGDAL: centroid
using CSV
using DataFrames
using YAXArrays

"""
    CScapeResultSet <: ResultSet

Result set for C~scape model outputs, loaded via `load_results`.

C~scape writes one NetCDF per scenario. Result variables are stored with dimensions
`(year, reef_sites, intervened, ft, thermal_tolerance)`, plus a leading `draws` dimension
when a scenario carries multiple draws. `intervened` has length 2 (counterfactual,
intervened) and `ft` is the functional-type axis. On load these are renamed and reordered
to ADRIA's `(timesteps, groups, locations, scenarios)` convention by `_reformat_cube`;
`_load_variable!` covers how individual outcomes are aggregated.

Only `relative_cover` is populated eagerly. Other outcomes are computed on first request by
the C~scape `ADRIA.metrics.*` methods and cached in `outcomes`.

# Fields
- `name`: Dataset title, from the first NetCDF's `title` global attribute
- `RCP`: RCP string (last two digits of the `ssp` attribute)
- `loc_ids`: Location identifiers, ordered to match `initial_cover.csv`
- `loc_area`: Location areas (m²)
- `loc_max_coral_cover`: Carrying capacity `k` per location, as a proportion
- `loc_centroids`: Location centroid points
- `env_layer_md`: Environmental-layer / data-package paths and timeframe
- `connectivity_data`: Connectivity matrix (`DataFrame`), rows/cols ordered to `loc_ids`
- `loc_data`: Location spatial attributes from the GeoPackage
- `raw_data`: Open handles to the source NetCDF files
- `inputs`: Reconstructed scenario input table (one row per scenario/draw)
- `sim_constants`: `SimConstants`
- `model_spec`: Partial model specification derived from `inputs`
- `outcomes`: Cache of computed outcome `YAXArray`s, keyed by metric name
- `coral_size_diameter`: Coral size-class diameters (`ft` × size class); C~scape uses its
  own size classes
"""
struct CScapeResultSet <: ResultSet
    name::String
    RCP::String

    loc_ids
    loc_area::Vector{Float64}
    loc_max_coral_cover::Vector{Float64}
    loc_centroids
    env_layer_md::EnvLayer
    connectivity_data::DataFrame
    loc_data::DataFrame
    raw_data::Vector{NcFile}

    inputs::DataFrame
    sim_constants::SimConstants
    model_spec::DataFrame

    outcomes::Dict{Symbol,YAXArray}
    coral_size_diameter::YAXArray
end

"""
    load_results(::Type{CScapeResultSet}, data_dir::String; result_dir::String=joinpath(data_dir, "results"), result_files::Vector{String}=String[], show_progress::Bool=true)::CScapeResultSet

Interface for loading C~scape model outputs.

See the [Loading C~scape Results](@ref) section for details on expected directory structure.

# Arguments
- `data_dir`: Path to the C~scape data package (`ScenarioID.csv`, `connectivity/`, `site_data/`, `initial_cover/`)

# Keyword Arguments
- `result_dir`: Directory of result NetCDFs. Defaults to the `results` subdirectory of `data_dir`. Ignored when `result_files` is given
- `result_files`: Explicit list of result NetCDF paths to load, instead of scanning `result_dir`
- `show_progress`: Show a progress bar while computing outcomes. Defaults to `true`

# Returns
CScapeResultSet struct compatible with most ADRIA analysis functionality.

# Examples
```julia
## Result NetCDFs read from the data package's own `results/` subdirectory
rs = ADRIA.load_results(CScapeResultSet, "a C~scape data package")

## Result NetCDFs held in a separate directory
rs = ADRIA.load_results(CScapeResultSet, "a C~scape data package"; result_dir="path/to/netcdfs")

## Explicit list of result NetCDFs
rs = ADRIA.load_results(
    CScapeResultSet, "a C~scape data package"; result_files=["NetCDF_Scn_140001.nc"]
)
```
"""
function load_results(
    ::Type{CScapeResultSet}, data_dir::String;
    result_dir::String=joinpath(data_dir, "results"),
    result_files::Vector{String}=String[],
    show_progress::Bool=true
)::CScapeResultSet
    !isdir(data_dir) ? error("Expected a directory but received $(data_dir)") : nothing

    isempty(result_files) && (result_files = _get_result_paths(result_dir))

    scenario_spec_path::String = joinpath(data_dir, "ScenarioID.csv")
    scenario_spec::DataFrame = DataFrame(
        CSV.File(scenario_spec_path; missingstring="NA")
    )

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
    geom_col = _get_geom_col(geodata)
    location_centroids = [centroid(multipoly) for multipoly ∈ geodata[!, geom_col]]

    outcomes = Dict{Symbol,YAXArray}()

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
    _get_scenario_id(nc_handle::NcFile)::Int

Return the scenario id from the NetCDF's `scenario_ID` global attribute, mapped into the
range used by `ScenarioID.csv`. Ids already at or above 100000 are returned unchanged;
smaller ids (older outputs that numbered scenarios from 1) are offset by 100000.

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
    _ft_acronym(ft_name::String)::String

Convert the full functional type names to acronyms for comptability with ADRIA.
"""
function _ft_acronym(ft_name::String)::String
    ft_name_map::Dict{String,String} = Dict{String,String}(
        "acro_corym" => "CA",
        "small_massive" => "SM",
        "corym_non_acro" => "CNA",
        "acro_table" => "TA",
        "large_massive" => "LM"
    )
    return ft_name_map[ft_name]
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
    rcps::Vector{Float64} = parse.(
        Float64, [nc_handle.gatts["ssp"][(end - 1):end] for nc_handle in nc_handles]
    )

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

    fragmented_spec::Vector{DataFrame} = _create_inputs_dataframe.(
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

    dhws::Vector{Float64} = Float64.(fill(input_index, n_draws))
    cyclones::Vector{Float64} = Float64.(fill(input_index, n_draws))
    scenario_id::Vector{Float64} = Int64.(fill(scenario_spec.ID, n_draws))

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
    plasticity::Float64 = scenario_spec.Plasticity

    settle_probability = parse.(Float64, split(scenario_spec.settle_prob, "_"))

    settle_probs_kwargs = Dict(
        Symbol(ft * "_settle_probability") => prob
        for (ft, prob) in zip(functional_types, settle_probability)
    )

    # Intervention factors
    deployment_area = _default_missing(scenario_spec[Symbol("Deployment area")], 0.0)
    total_corals = Float64.(_default_missing(scenario_spec.TotalCorals, 1.0))

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
        Symbol("N_seed_" * _ft_acronym(ft)) => n_corals
        for (ft, n_corals) in zip(functional_types, n_seeded)
    )

    enhancement_mean, enhancement_std = parse.(
        Float64, split(_default_missing(scenario_spec.Enhancement, "0_0"), '_')
    )

    intervention_start = Float64.(
        _default_missing(
            scenario_spec.InterventionYears_start, nc_handle["year"][1]
        ) - nc_handle["year"][1]
    )

    duration = Float64.(_default_missing(scenario_spec.duration, 0.0))
    frequency = Float64.(_default_missing(scenario_spec.frequency, 0.0))
    n_dep_locations = Float64.(
        count(x -> x == '/', _default_missing(scenario_spec.Reef_siteids, ""))
    )

    coral_dict = merge(settle_probs_kwargs, corals_deployed)

    # `scenario_id`, `dhws` and `cyclones` are length-`n_draws` vectors; every other
    # column is a per-scenario scalar that `DataFrame` recycles to match.
    return DataFrame(;
        scenario_id=scenario_id,
        dhw_scenario=dhws,
        cyc_scenario=cyclones,
        RCP=rcp,
        thermal_tol_lb=lb_t,
        thermal_tol_int1=int_t1,
        thermal_tol_int2=int_t2,
        thermal_tol_ub=ub_t,
        init_heat_tol_mean=init_heat_tol_mean,
        init_heat_tol_std=init_heat_tol_std,
        heritability1=heritability1,
        heritability2=heritability2,
        plasticity=plasticity,
        intervention_start=intervention_start,
        intervention_duration=duration,
        intervention_frequency=frequency,
        deployment_area=deployment_area,
        n_deployment_locations=n_dep_locations,
        enhancement_mean=enhancement_mean,
        enhancement_std=enhancement_std,
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
    seeded_names = filter(factor -> contains(factor, "N_seed_"), factor_names)
    seeded_readable = human_readable_name.(seeded_names)
    seeded_names = Symbol.(seeded_names)
    seeded_ptype = repeat(["continuous"], length(settle_names))
    seeded_lb = repeat([0.0], length(seeded_names))
    seeded_ub = repeat([5e7], length(seeded_names))

    # Construct default model spec
    fieldname::Vector{Symbol} = [
        :scenario_id,
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
        "Scenario ID",
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
        "Scenario ID",
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
        150000.0,
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
        repeat(["EnvironmentalLayer"], 3),
        repeat(["Coral"], 9),
        repeat(["Intervention"], 7),
        repeat(["Coral"], length(settle_names)),
        repeat(["Intervention"], length(seeded_names))
    )
    is_constant::Vector{Bool} = fill(false, length(component))
    return DataFrame(;
        component=component,
        fieldname=fieldname,
        description=descriptions,
        name=human_names,
        ptype=ptype,
        dist_params=dist_params,
        lower_bound=lower_bound,
        upper_bound=upper_bound,
        is_constant=is_constant
    )
end

"""
    _get_result_paths(result_dir::String)::Vector{String}

Get the names of all result netcdf result files in the given directory. Assumes filenames
match "NetCDF_Scn.*.nc".
"""
function _get_result_paths(result_dir::String)::Vector{String}
    possible_files = filter(isfile, readdir(result_dir; join=true))
    possible_files = filter(x -> occursin(r"NetCDF_Scn.*.nc", x), possible_files)
    if length(possible_files) == 0
        msg = "Unable to find result netcdf file in subdirectory \"$(result_dir)\""
        throw(ArgumentError(msg))
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
    _get_result_name(ds::NcFile)::String

Get the dataset name from the NetCDF's `title` global attribute, or `"CScape Results"` if
that attribute is missing.
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
    _reformat_cube(::Type{CScapeResultSet}, cscape_cube::YAXArray)::YAXArray

Rename and reorder a C~scape outcome cube's dimensions to ADRIA's convention.

C~scape names its axes `year`, `ft` (functional type), `reef_sites` and `draws`; ADRIA
expects `timesteps`, `groups`, `locations`, `scenarios` in that order. The `ft` axis maps
to ADRIA's `groups`: core metric metadata (`axes_units`, `human_readable_name`) only
recognises `groups`/`sizes`. Axes outside this mapping are kept and appended after the
renamed ones. A `units` property of `"percent"` has its values rescaled to a `[0, 1]`
proportion.
"""
function _reformat_cube(::Type{CScapeResultSet}, cscape_cube::YAXArray)::YAXArray
    dim_names = name.(cscape_cube.axes)
    # C~scape axis name => ADRIA axis name (see docstring for the `ft` -> `groups` mapping)
    cscape_names = [:year, :ft, :reef_sites, :draws]
    adria_names = [:timesteps, :groups, :locations, :scenarios]
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

function Base.show(io::IO, mime::MIME"text/plain", rs::CScapeResultSet)
    return _show_external_resultset(
        io, rs;
        count_label="NetCDFs loaded",
        count=length(rs.raw_data),
        n_locs=length(rs.loc_ids)
    )
end
