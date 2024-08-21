const COMPRESSOR = Zarr.BloscCompressor(; cname="zstd", clevel=2, shuffle=true)

"""
    summarize_env_data(data_cube::AbstractArray)

Summarize environmental data layers (mean and standard deviation).

# Returns
Matrix{Float64, 2}, of mean and standard deviation for each environmental scenario.
"""
function summarize_env_data(data::AbstractArray)::Array{Float64}
    # TODO: Update once
    stats_store::Array{Float64} = zeros(2, size(data, 3), size(data, 2))
    stats_store[1, :, :] .= dropdims(mean(data; dims=1); dims=1)'
    stats_store[2, :, :] .= dropdims(std(data; dims=1); dims=1)'
    return stats_store
end

"""
    store_env_summary(data_cube::AbstractArray, type::String, file_loc::String, rcp::String, compressor::Zarr.Compressor)::ZArray

Retrieve summary statistics matrices from DataFrames of dhws and waves.
Produce summary statistics (mean/std) for given data cube saved to a Zarr data store.

# Arguments
- `data_cube` : Data to summarize
- `type` : Dimension identifier to use
- `file_loc` : Path for Zarr data store
- `rcp` : RCP
- `compressor` : Zarr compressor

# Returns
Zarr data store holding a 2*N*M matrix.

First row is mean over time
Second row is the std over time
N is the number of dhw/wave scenarios.
M is the number of locations.
"""
function store_env_summary(
    data_cube::AbstractArray,
    type::String,
    file_loc::String,
    rcp::String,
    compressor::Zarr.Compressor
)::ZArray
    stats = summarize_env_data(data_cube)

    stats_store = zcreate(
        Float32,
        (2, size(stats, 2), size(stats, 3))...;
        fill_value=nothing, fill_as_missing=false,
        path=joinpath(file_loc, rcp),
        attrs=Dict(
            :structure => ("stat", type, "locations"),
            :stats => ["mean", "std"],
            :scenarios => string.(1:size(stats, 2)),
            :locations => string.(1:size(stats, 3)),
            :rcp => rcp),
        compressor=compressor)

    stats_store[:, :, :] .= stats

    return stats_store
end

"""
    store_conn(conn_data::YAXArray, file_loc::String, rcp::String,
        compressor::Zarr.Compressor)::ZArray

Retrieve connectivity matrices from Domain for storage.
Produce connectivity for a particular RCP saved to a Zarr data store.

# Arguments
- `conn_data` : connectivity data (e.g. `domain.conn`)
- `file_loc` : path for Zarr data store
- `rcp`: RCP associated with connectivity data.

# Returns
Zarr data store holding a M*M matrix.
M is the number of locations.
"""
function store_conn(
    conn_data::YAXArray,
    file_loc::String,
    rcp::String,
    compressor::Zarr.Compressor
)::ZArray
    conn_store = zcreate(
        Float32,
        size(conn_data)...;
        fill_value=nothing, fill_as_missing=false,
        path=joinpath(file_loc, rcp),
        attrs=Dict(
            :structure => ("Source", "Sink"),
            :Source => collect(conn_data.Source),
            :Sink => collect(conn_data.Sink),
            :rcp => rcp),
        compressor=compressor)

    conn_store[:, :] .= Matrix(conn_data)
    return conn_store
end

"""
    scenario_attributes(name, RCP, input_cols, invoke_time, env_layer, sim_constants, unique_loc_ids, area, k, centroids)
    scenario_attributes(domain::Domain, param_df::DataFrame)

Generate dictionary of scenario attributes.
"""
function scenario_attributes(
    name,
    RCP,
    input_cols,
    invoke_time,
    env_layer,
    sim_constants,
    unique_loc_ids,
    area,
    k,
    centroids
)::Dict{Symbol,Any}
    attrs::Dict{Symbol,Any} = Dict(
        :name => name,
        :RCP => RCP,
        :columns => input_cols,
        :invoke_time => invoke_time,
        :ADRIA_VERSION => "v" * string(PkgVersion.Version(@__MODULE__)),
        :site_data_file => env_layer.site_data_fn,
        :site_id_col => env_layer.site_id_col,
        :cluster_id_col => env_layer.cluster_id_col,
        :init_coral_cover_file => env_layer.init_coral_cov_fn,
        :connectivity_file => env_layer.connectivity_fn,
        :DHW_file => env_layer.DHW_fn,
        :wave_file => env_layer.wave_fn,
        :timeframe => env_layer.timeframe,
        :sim_constants => sim_constants,
        :site_ids => unique_loc_ids,
        :site_area => area,
        :site_max_coral_cover => k,
        :site_centroids => centroids
    )

    return attrs
end
function scenario_attributes(domain::Domain, param_df::DataFrame)::Dict{Symbol,Any}
    return scenario_attributes(
        domain.name,
        domain.RCP,
        names(param_df),
        domain.scenario_invoke_time,
        domain.env_layer_md,
        domain.sim_constants,
        unique_loc_ids(domain),
        site_area(domain),
        domain.site_data.k,
        centroids(domain.site_data)
    )
end

"""
    setup_logs(z_store, unique_loc_ids, n_scens, tf, n_locs, n_group_and_size)

Setup logs for ranks, seed_log, fog_log, shade_log and coral_dhw_log.

# Arguments
- `z_store` : ZArray
- `unique_loc_ids` : Unique location ids
- `n_scens` : number of scenarios
- `tf` : timeframe
- `n_locs` : number of location
- `n_group_and_size` : number of function groups ⋅ number of size classes

Note: This setup relies on hardcoded values for number of species represented and seeded.
"""
function setup_logs(z_store, unique_loc_ids, n_scens, tf, n_locs, n_group_and_size)
    # Set up logs for location ranks, seed/fog log
    zgroup(z_store, LOG_GRP)
    log_fn::String = joinpath(z_store.folder, LOG_GRP)

    # Store ranked location
    n_interventions = length(interventions())
    rank_dims::Tuple{Int64,Int64,Int64,Int64} = (tf, n_locs, n_interventions, n_scens)  # locations, location id and rank, no. scenarios
    fog_dims::Tuple{Int64,Int64,Int64} = (tf, n_locs, n_scens)  # timeframe, location, no. scenarios

    # tf, no. species to seed, location id and rank, no. scenarios
    seed_dims::Tuple{Int64,Int64,Int64,Int64} = (tf, 3, n_locs, n_scens)

    attrs = Dict(
        # Here, "intervention" refers to seeding or shading
        :structure => ("timesteps", "locations", "intervention", "scenarios"),
        :unique_site_ids => unique_loc_ids
    )
    ranks = zcreate(
        Float32,
        rank_dims...;
        name="rankings",
        fill_value=nothing,
        fill_as_missing=false,
        path=log_fn,
        chunks=(rank_dims[1:3]..., 1),
        attrs=attrs
    )

    attrs = Dict(
        :structure => ("timesteps", "coral_id", "locations", "scenarios"),
        :unique_site_ids => unique_loc_ids
    )
    seed_log = zcreate(
        Float32,
        seed_dims...;
        name="seed",
        fill_value=nothing,
        fill_as_missing=false,
        path=log_fn,
        chunks=(seed_dims[1:3]..., 1),
        attrs=attrs
    )

    attrs = Dict(
        :structure => ("timesteps", "locations", "scenarios"),
        :unique_site_ids => unique_loc_ids
    )
    fog_log = zcreate(
        Float32,
        fog_dims...;
        name="fog",
        fill_value=nothing,
        fill_as_missing=false,
        path=log_fn,
        chunks=(fog_dims[1:2]..., 1),
        attrs=attrs
    )
    shade_log = zcreate(
        Float32,
        fog_dims...;
        name="shade",
        fill_value=nothing,
        fill_as_missing=false,
        path=log_fn,
        chunks=(fog_dims[1:2]..., 1),
        attrs=attrs
    )

    # TODO: Could log bleaching mortality
    # attrs = Dict(
    #     :structure => ("timesteps", "locations", "scenarios"),
    #     :unique_site_ids => unique_loc_ids,
    # )
    # bleach_log = zcreate(Float32, fog_dims...; name="bleaching_mortality", fill_value=nothing, fill_as_missing=false, path=log_fn, chunks=(fog_dims[1:2]..., 1), attrs=attrs)

    # Log for coral DHW thresholds
    attrs = Dict(
        :structure => ("timesteps", "species", "locations", "scenarios"),
        :unique_site_ids => unique_loc_ids
    )

    local coral_dhw_log
    if parse(Bool, ENV["ADRIA_DEBUG"]) == true
        coral_dhw_log = zcreate(
            Float32,
            tf,
            n_group_and_size,
            n_locs,
            n_scens;
            name="coral_dhw_log",
            fill_value=nothing,
            fill_as_missing=false,
            path=log_fn,
            chunks=(tf, n_group_and_size, n_locs, 1),
            attrs=attrs
        )
    else
        coral_dhw_log = zcreate(
            Float32,
            tf,
            n_group_and_size,
            1,
            n_scens;
            name="coral_dhw_log",
            fill_value=0.0,
            fill_as_missing=false,
            path=log_fn,
            chunks=(tf, n_group_and_size, 1, 1),
            attrs=attrs
        )
    end

    return ranks, seed_log, fog_log, shade_log, coral_dhw_log
end

"""
    setup_result_store!(domain::Domain, scen_spec::DataFrame)

Sets up an on-disk result store.

## Structure of Store

```
├───connectivity
├───env_stats
├───inputs
├───logs
|   ├───coral_dhw_log
│   ├───fog
│   ├───rankings
│   ├───seed
│   └───shade
├───model_spec
├───results
|   ├───absolute_shelter_volume
|   ├───coral_evenness
|   ├───juvenile_indicator
│   ├───relative_cover
|   ├───relative_juveniles
│   ├───relative_shelter_volume
│   └───relative_taxa_cover
└───spatial
```

- `inputs` : includes domain specification metadata including what connectivity/DHW/wave data was used.
- `model_spec` : contains a copy of the ADRIA model specification (as CSV).
- `spatial` : contains a copy of the spatial domain data (as a geopackage).

# Notes
- `domain` is replaced with an identical copy with an updated scenario invoke time.
- -9999.0 is used as an arbitrary fill value.

# Arguments
- `domain` : ADRIA scenario domain
- `scen_spec` : ADRIA scenario specification

# Returns
domain, (relative_cover, relative_shelter_volume, absolute_shelter_volume, relative_juveniles,
juvenile_indicator, relative_taxa_cover, site_ranks, seed_log, fog_log, shade_log)
"""
function setup_result_store!(domain::Domain, scen_spec::DataFrame)::Tuple
    @set! domain.scenario_invoke_time = replace(
        string(now()), "T" => "_", ":" => "_", "." => "_"
    )

    # Collect defined RCPs
    rcps = string.(unique(scen_spec, "RCP")[!, "RCP"])
    log_location::String = _result_location(domain, rcps)

    z_store = DirectoryStore(log_location)

    # Store copy of inputs
    input_loc::String = joinpath(z_store.folder, INPUTS)
    input_dims::Tuple{Int64,Int64} = size(scen_spec)
    attrs::Dict = scenario_attributes(domain, scen_spec)

    # Write a copy of spatial data to the result set
    mkdir(joinpath(log_location, SPATIAL_DATA))
    geo_fn = joinpath(log_location, SPATIAL_DATA, basename(attrs[:name]) * ".gpkg")

    # Writing geopackages out does not currently automatically include CRS
    # so we manually define it as a workaround.
    col = _get_geom_col(domain.site_data)
    ref = AG.getspatialref(domain.site_data[1, col])
    proj4_gft = GFT.ProjString(AG.toPROJ4(ref))
    GDF.write(geo_fn, domain.site_data; crs=proj4_gft, geom_columns=(col,), driver="GPKG")

    # Store copy of model specification as CSV
    mkdir(joinpath(log_location, "model_spec"))
    model_spec(domain, joinpath(log_location, "model_spec", "model_spec.csv"))

    # Create store for scenario spec
    inputs = zcreate(
        Float64,
        input_dims...;
        fill_value=-9999.0,
        fill_as_missing=false,
        path=input_loc,
        chunks=input_dims,
        attrs=attrs
    )

    # Store table of factor values
    inputs[:, :] = Matrix(scen_spec)

    tf, n_locations, _ = size(domain.dhw_scens)
    n_scenarios = nrow(scen_spec)

    # Set up stores for each metric
    function dim_lengths(metric_structure)
        dl = []
        for d in metric_structure
            if d == "timesteps"
                append!(dl, tf)
            elseif d == "species"
                append!(dl, domain.coral_growth.n_group_and_size)
            elseif d == "locations"
                append!(dl, n_locations)
            elseif d == "scenarios"
                append!(dl, n_scenarios)
            end
        end

        return (dl...,)
    end

    met_names = [
        :relative_cover,
        :relative_shelter_volume,
        :absolute_shelter_volume,
        :relative_juveniles,
        :juvenile_indicator,
        :coral_evenness
    ]

    _unique_loc_ids = unique_loc_ids(domain)

    dim_struct = Dict(
        :structure => string.((:timesteps, :locations, :scenarios)),
        :unique_site_ids => _unique_loc_ids
    )
    result_dims::Tuple{Int64,Int64,Int64} = dim_lengths(dim_struct[:structure])

    # Create stores for each metric
    stores = [
        zcreate(
            Float32,
            result_dims...;
            fill_value=nothing,
            fill_as_missing=false,
            path=joinpath(z_store.folder, RESULTS, string(m_name)),
            chunks=(result_dims[1:(end - 1)]..., 1),
            attrs=dim_struct,
            compressor=COMPRESSOR
        ) for m_name in met_names
    ]

    # Handle special case for relative taxa cover
    n_groups::Int64 = domain.coral_growth.n_groups
    push!(
        stores,
        zcreate(
            Float32,
            (result_dims[1], n_groups, result_dims[3])...;
            fill_value=nothing,
            fill_as_missing=false,
            path=joinpath(z_store.folder, RESULTS, "relative_taxa_cover"),
            chunks=((result_dims[1], n_groups)..., 1),
            attrs=Dict(:structure => string.(ADRIA.metrics.relative_taxa_cover.dims)),
            compressor=COMPRESSOR
        )
    )
    push!(met_names, :relative_taxa_cover)

    # dhw and wave zarrays
    dhw_stats = []
    wave_stats = []
    connectivity = []
    dhw_stat_names = []
    wave_stat_names = []
    conn_names = []
    for rcp in rcps
        push!(
            dhw_stats,
            store_env_summary(
                domain.dhw_scens,
                "dhw_scenario",
                joinpath(z_store.folder, ENV_STATS, "dhw"),
                rcp,
                COMPRESSOR
            )
        )
        push!(
            wave_stats,
            store_env_summary(
                domain.wave_scens,
                "wave_scenario",
                joinpath(z_store.folder, ENV_STATS, "wave"),
                rcp,
                COMPRESSOR
            )
        )
        push!(
            connectivity,
            store_conn(
                domain.conn,
                joinpath(z_store.folder, "connectivity"),
                rcp,
                COMPRESSOR
            )
        )

        push!(dhw_stat_names, Symbol("dhw_stat_$rcp"))
        push!(wave_stat_names, Symbol("wave_stat_$rcp"))
        push!(conn_names, Symbol("connectivity_$rcp"))
    end
    stat_store_names = vcat(dhw_stat_names, wave_stat_names)

    n_group_and_size::Int64 = domain.coral_growth.n_group_and_size
    # Group all data stores
    stores = [
        stores...,
        dhw_stats...,
        wave_stats...,
        connectivity...,
        setup_logs(
            z_store, _unique_loc_ids, nrow(scen_spec), tf, n_locations, n_group_and_size
        )...
    ]

    return domain,
    (;
        zip(
            (
                met_names...,
                stat_store_names...,
                conn_names...,
                :site_ranks,
                :seed_log,
                :fog_log,
                :shade_log,
                :coral_dhw_log
            ),
            stores
        )...
    )
end

"""
    _recreate_stats_from_store(zarr_store_path::String)::Dict{String, AbstractArray}

Recreate data structure holding RCP summary statistics from Zarr store.
"""
function _recreate_stats_from_store(zarr_store_path::String)::Dict{String,YAXArray}
    rcp_dirs = filter(d -> isdir(joinpath(zarr_store_path, d)), readdir(zarr_store_path))
    rcp_stat_dirs = joinpath.(zarr_store_path, rcp_dirs)

    stat_d = Dict{String,YAXArray}()
    for (i, sd) in enumerate(rcp_stat_dirs)
        store = zopen(sd; fill_as_missing=false)

        dim_names = Symbol.(store.attrs["structure"])
        stats = string.(store.attrs["stats"])
        scenario_ids = string.(store.attrs["scenarios"])
        loc_ids = string.(store.attrs["locations"])
        stat_set = DataCube(
            store[:, :, :];
            zip(dim_names, [stats, scenario_ids, loc_ids])...
        )

        stat_d[rcp_dirs[i]] = stat_set
    end

    return stat_d
end

"""
    _recreate_conn_from_store(zarr_store_path::String)::Dict{String, AbstractArray}

Recreate data structure holding connectivity for each RCP from Zarr store.
"""
function _recreate_conn_from_store(zarr_store_path::String)::Dict{String,YAXArray}
    rcp_dirs = filter(d -> isdir(joinpath(zarr_store_path, d)), readdir(zarr_store_path))
    rcp_stat_dirs = joinpath.(zarr_store_path, rcp_dirs)

    conn_d = Dict{String,YAXArray}()
    for (i, sd) in enumerate(rcp_stat_dirs)
        store = zopen(sd; fill_as_missing=false)

        dim_names = Symbol.(store.attrs["structure"])
        source_ids = string.(store.attrs["Source"])
        sink_ids = string.(store.attrs["Sink"])
        conn_set = DataCube(
            store[:, :];
            zip(dim_names, [source_ids, sink_ids])...
        )

        conn_d[rcp_dirs[i]] = conn_set
    end

    return conn_d
end

"""
    load_results(result_loc::String)::ResultSet
    load_results(domain::Domain)::ResultSet

Create interface to a given Zarr result set.
"""
function load_results(result_loc::String)::ResultSet
    !isdir(result_loc) ? error("Not a directory: $(result_loc)") : nothing

    # Read in results
    local raw_set
    try
        raw_set = zopen(joinpath(result_loc, RESULTS); fill_as_missing=false)
    catch err
        if !occursin("ArgumentError", sprint(showerror, err))
            rethrow(err)
        end
    end

    # Read in logs
    log_set = zopen(joinpath(result_loc, LOG_GRP); fill_as_missing=false)
    input_set = zopen(joinpath(result_loc, INPUTS); fill_as_missing=false)

    dhw_stat_set = _recreate_stats_from_store(joinpath(result_loc, ENV_STATS, "dhw"))
    wave_stat_set = _recreate_stats_from_store(joinpath(result_loc, ENV_STATS, "wave"))
    conn_set = _recreate_conn_from_store(joinpath(result_loc, "connectivity"))

    result_loc = replace(result_loc, "\\" => "/")
    if endswith(result_loc, "/")
        result_loc = result_loc[1:(end - 1)]
    end

    # Spatial data
    site_data = GDF.read(
        joinpath(result_loc, SPATIAL_DATA, input_set.attrs["name"] * ".gpkg")
    )
    sort!(site_data, [Symbol(input_set.attrs["site_id_col"])])

    # Model specification
    model_spec = CSV.read(
        joinpath(result_loc, MODEL_SPEC, "model_spec.csv"), DataFrame; comment="#"
    )

    # Standardize fieldnames to Symbol
    # TODO: Match all other column data types with original model spec
    model_spec.fieldname .= Symbol.(model_spec.fieldname)

    r_vers_id = input_set.attrs["ADRIA_VERSION"]
    t_vers_id = "v" * string(PkgVersion.Version(@__MODULE__))

    if r_vers_id != t_vers_id
        msg = """Results were produced with a different version of ADRIA ($(r_vers_id)).
        The version of ADRIA in use is $(t_vers_id).\n
        Errors may occur when analyzing data.

        ADRIA v0.8 store results relative to absolute location area, where as v0.9+ now
        stores results relative to available area.
        """

        @warn msg
    end

    # The inputs used
    input_cols::Array{String} = input_set.attrs["columns"]
    inputs_used::DataFrame = DataFrame(input_set[:, :], input_cols)

    # Details of the environmental data layer used for the sims
    env_layer_md::EnvLayer = EnvLayer(
        result_loc,
        input_set.attrs["site_data_file"],
        input_set.attrs["site_id_col"],
        input_set.attrs["cluster_id_col"],
        input_set.attrs["init_coral_cover_file"],
        input_set.attrs["connectivity_file"],
        input_set.attrs["DHW_file"],
        input_set.attrs["wave_file"],
        input_set.attrs["timeframe"]
    )

    outcomes = Dict{Symbol,YAXArray}()
    subdirs = filter(isdir, readdir(joinpath(result_loc, RESULTS); join=true))
    for sd in subdirs
        if !(occursin(LOG_GRP, sd)) && !(occursin(INPUTS, sd))
            res = zopen(sd; fill_as_missing=false)
            sz = size(res)

            # Construct dimension names and metadata
            st = []
            for (i, s) in enumerate(res.attrs["structure"])
                if s == "timesteps"
                    push!(st, input_set.attrs["timeframe"])
                elseif s == "locations"
                    push!(st, res.attrs["unique_site_ids"])
                else
                    push!(st, 1:sz[i])
                end
            end

            try
                outcomes[Symbol(basename(sd))] = DataCube(
                    res; zip(Symbol.(res.attrs["structure"]), st)...
                )
            catch err
                if err isa ArgumentError
                    @warn """Unable to resolve result structure, reverting to numbered keys.
                    For group $(sd)
                    Got: $(size(res))
                    Structure: $(res.attrs["structure"])
                    Generated: $(Array([i[1] for i in size.(st)]))
                    """
                    outcomes[Symbol(basename(sd))] = DataCube(
                        res;
                        zip(Symbol.(res.attrs["structure"]), [1:s for s in size(res)])...
                    )
                else
                    rethrow(err)
                end
            end
        end
    end

    return ResultSet(
        input_set,
        env_layer_md,
        inputs_used,
        outcomes,
        log_set,
        dhw_stat_set,
        wave_stat_set,
        conn_set,
        site_data,
        model_spec
    )
end
function load_results(domain::Domain)::ResultSet
    return load_results(result_location(domain))
end

"""
    result_location(d::Domain)::String

Generate path to the data store of results for the given Domain.
"""
function result_location(d::Domain)::String
    return joinpath(
        ENV["ADRIA_OUTPUT_DIR"], "$(d.name)__RCPs$(d.RCP)__$(d.scenario_invoke_time)"
    )
end

"""
    _result_location(d::Domain, rcps::Vector{String})::String

Generate path to the data store of results for the given Domain and RCPs names.
"""
function _result_location(d::Domain, rcps::Vector{String})::String
    return joinpath(
        ENV["ADRIA_OUTPUT_DIR"],
        "$(d.name)__RCPs_$(join(rcps, "_"))__$(d.scenario_invoke_time)"
    )
end
