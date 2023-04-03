function get_geometry(df::DataFrame)
    if columnindex(df, :geometry) > 0
        return df.geometry
    elseif columnindex(df, :geom) > 0
        return df.geom
    end

    error("No geometry data found")
end

"""
    centroids(df::DataFrame)

Extract and return long/lat from a GeoDataFrame.

# Arguments
- `df` : GeoDataFrame

# Returns
Array of tuples (x, y), where x and y relate to long and lat respectively.
"""
function centroids(df::DataFrame)::Vector{Tuple{Float64,Float64}}
    site_centroids::Vector = AG.centroid.(get_geometry(df))
    return collect(zip(AG.getx.(site_centroids, 0), AG.gety.(site_centroids, 0)))
end


"""
    summarize_env_data(data_cube::AbstractArray)

Summarize environmental data layers (mean and standard deviation).

# Returns
Matrix{Float64, 2}, of mean and standard deviation for each environmental scenario.
"""
function summarize_env_data(data::AbstractArray)::Array{Float64}
    # TODO: Update once
    dc_mean = dropdims(mean(data, dims=(1, 2)), dims=(1, 2))'
    dc_std = dropdims(std(data, dims=(1, 2)), dims=(1, 2))'

    return Array{Float64}(vcat(dc_mean, dc_std))
end

"""
    store_env_summary(data_cube::AbstractArray, type::String, file_loc::String, compressor::Zarr.Compressor)

Retrieve summary statistics matrices from DataFrames of dhws and waves.
Produce summary statistics (mean/std) for given data cube saved to a Zarr data store.

# Arguments
- `data_cube` : data to summarize
- `type` : dimension identifier to use
- `file_loc` : path for Zarr data store

# Returns
Zarr data store holding a 2*N matrix.

First row is mean over time
Second row is the std over time
N is the number of dhw/wave scenarios.
"""
function store_env_summary(data_cube::NamedDimsArray, type::String, file_loc::String, rcp::String, compressor::Zarr.Compressor)::ZArray
    stats = summarize_env_data(data_cube)

    stats_store = zcreate(Float32, (2, size(stats, 2))...;
        fill_value=nothing, fill_as_missing=false,
        path=joinpath(file_loc, rcp),
        attrs=Dict(
            :structure => ("stat", type),
            :rows => ["mean", "std"],
            :cols => string.(1:size(stats, 2)),
            :rcp => rcp),
        compressor=compressor)

    stats_store[:, :] .= stats

    return stats_store
end

"""
    scenario_attributes(name, RCP, input_cols, invoke_time, env_layer, sim_constants, unique_sites, area, k)
    scenario_attributes(domain::Domain, param_df::DataFrame)

Generate dictionary of scenario attributes.
"""
function scenario_attributes(name, RCP, input_cols, invoke_time, env_layer, sim_constants, unique_sites, area, k, centroids)::Dict
    attrs::Dict{Symbol,Any} = Dict(
        :name => name,
        :RCP => RCP,
        :columns => input_cols,
        :invoke_time => invoke_time,
        :ADRIA_VERSION => "v" * string(PkgVersion.Version(@__MODULE__)), :site_data_file => env_layer.site_data_fn,
        :site_id_col => env_layer.site_id_col,
        :unique_site_id_col => env_layer.unique_site_id_col,
        :init_coral_cover_file => env_layer.init_coral_cov_fn,
        :connectivity_file => env_layer.connectivity_fn,
        :DHW_file => env_layer.DHW_fn,
        :wave_file => env_layer.wave_fn,
        :timeframe => env_layer.timeframe,
        :sim_constants => sim_constants,
        :site_ids => unique_sites,
        :site_area => area,
        :site_max_coral_cover => k,
        :site_centroids => centroids
    )

    return attrs
end
function scenario_attributes(domain::Domain, param_df::DataFrame)
    return scenario_attributes(domain.name, domain.RCP, names(param_df), domain.scenario_invoke_time,
        domain.env_layer_md, domain.sim_constants,
        unique_sites(domain), domain.site_data.area, domain.site_data.k, centroids(domain.site_data))
end


function setup_logs(z_store, unique_sites, n_scens, tf, n_sites)
    # Set up logs for site ranks, seed/fog log
    zgroup(z_store, LOG_GRP)
    log_fn::String = joinpath(z_store.folder, LOG_GRP)

    # Store ranked sites
    rank_dims::Tuple{Int64,Int64,Int64,Int64} = (tf, n_sites, 2, n_scens)  # sites, site id and rank, no. scenarios
    fog_dims::Tuple{Int64,Int64,Int64} = (tf, n_sites, n_scens)  # timeframe, sites, no. scenarios

    # tf, no. intervention sites, site id and rank, no. scenarios
    seed_dims::Tuple{Int64,Int64,Int64,Int64} = (tf, 2, n_sites, n_scens)

    attrs = Dict(
        # Here, "intervention" refers to seeding or shading
        :structure => ("timesteps", "sites", "intervention", "scenarios"),
        :unique_site_ids => unique_sites,
    )
    ranks = zcreate(Float32, rank_dims...; name="rankings", fill_value=nothing, fill_as_missing=false, path=log_fn, chunks=(rank_dims[1:3]..., 1), attrs=attrs)

    attrs = Dict(
        :structure => ("timesteps", "coral_id", "sites", "scenarios"),
        :unique_site_ids => unique_sites,
    )
    seed_log = zcreate(Float32, seed_dims...; name="seed", fill_value=nothing, fill_as_missing=false, path=log_fn, chunks=(seed_dims[1:3]..., 1), attrs=attrs)

    attrs = Dict(
        :structure => ("timesteps", "sites", "scenarios"),
        :unique_site_ids => unique_sites,
    )
    fog_log = zcreate(Float32, fog_dims...; name="fog", fill_value=nothing, fill_as_missing=false, path=log_fn, chunks=(fog_dims[1:2]..., 1), attrs=attrs)
    shade_log = zcreate(Float32, fog_dims...; name="shade", fill_value=nothing, fill_as_missing=false, path=log_fn, chunks=(fog_dims[1:2]..., 1), attrs=attrs)

    return ranks, seed_log, fog_log, shade_log
end

"""
    setup_result_store!(domain::Domain, param_df::DataFrame)

Sets up an on-disk result store.

## Structure of Store

```
├───logs
│   ├───fog
│   ├───rankings
│   ├───seed
│   └───shade
├───results
│   ├───relative_cover
│   ├───relative_shelter_volume
│   └───absolute_shelter_volume
├───site_data
├───model_spec
└───inputs
```

- `inputs` : includes domain specification metadata including what connectivity/DHW/wave data was used.
- `site_data` : contains a copy of the spatial domain data (as geopackage).
- `model_spec` : contains a copy of the ADRIA model specification (as CSV).

# Notes
- `domain` is replaced with an identical copy with an updated scenario invoke time.
- -9999.0 is used as an arbitrary fill value.

# Arguments
- `domain` : ADRIA scenario domain
- `param_df` : ADRIA scenario specification

# Returns
domain, (relative_cover, relative_shelter_volume, absolute_shelter_volume, site_ranks, seed_log, fog_log, shade_log)
"""
function setup_result_store!(domain::Domain, param_df::DataFrame)::Tuple
    if "RCP" in names(param_df)
        param_df = param_df[:, Not("RCP")]  # Ignore RCP column if it exists
    end

    # TODO: Support setting up a combined result store.

    # Insert RCP column and populate with this dataset's RCP
    insertcols!(param_df, 1, :RCP => parse(Float64, domain.RCP))

    @set! domain.scenario_invoke_time = replace(string(now()), "T" => "_", ":" => "_", "." => "_")
    log_location::String = joinpath(ENV["ADRIA_OUTPUT_DIR"], "$(domain.name)__RCPs$(domain.RCP)__$(domain.scenario_invoke_time)")

    z_store = DirectoryStore(log_location)

    # Store copy of inputs
    input_loc::String = joinpath(z_store.folder, INPUTS)
    input_dims::Tuple{Int64,Int64} = size(param_df)
    attrs::Dict = scenario_attributes(domain, param_df)

    # Copy site data into result set
    mkdir(joinpath(log_location, "site_data"))
    cp(attrs[:site_data_file], joinpath(log_location, "site_data", basename(attrs[:site_data_file])), force=true)

    inputs = zcreate(Float64, input_dims...; fill_value=-9999.0, fill_as_missing=false, path=input_loc, chunks=input_dims, attrs=attrs)

    # Store copy of model specification as CSV
    mkdir(joinpath(log_location, "model_spec"))
    model_spec(domain, joinpath(log_location, "model_spec", "model_spec.csv"))

    # Store post-processed table of input parameters.
    # +1 skips the RCP column
    integer_params = findall(domain.model[:ptype] .== "integer")
    map_to_discrete!(param_df[:, integer_params.+1], getindex.(domain.model[:bounds], 2)[integer_params])
    inputs[:, :] = Matrix(param_df)

    tf, n_sites, _ = size(domain.dhw_scens)

    # Set up stores for each metric
    function dim_lengths(metric_structure)
        dl = []
        for d in metric_structure
            if d == "timesteps"
                append!(dl, tf)
            elseif d == "species"
                append!(dl, domain.coral_growth.n_species)
            elseif d == "sites"
                append!(dl, n_sites)
            elseif d == "scenarios"
                append!(dl, nrow(param_df))
            end
        end

        return (dl...,)
    end

    compressor = Zarr.BloscCompressor(cname="zstd", clevel=2, shuffle=true)

    met_names = [:total_absolute_cover, :relative_shelter_volume,
        :absolute_shelter_volume, :relative_juveniles, :juvenile_indicator]

    dim_struct = Dict(
        :structure => string.((:timesteps, :sites, :scenarios)),
        :unique_site_ids => unique_sites(domain)
    )
    result_dims::Tuple{Int64,Int64,Int64} = dim_lengths(dim_struct[:structure])

    # Create stores for each metric
    stores = [
        zcreate(Float32, result_dims...;
            fill_value=nothing, fill_as_missing=false,
            path=joinpath(z_store.folder, RESULTS, string(m_name)), chunks=(result_dims[1:end-1]..., 1),
            attrs=dim_struct,
            compressor=compressor)
        for m_name in met_names
    ]

    # Handle special case for relative taxa cover
    push!(
        stores,
        zcreate(Float32, (result_dims[1], 6, result_dims[3])...;
            fill_value=nothing, fill_as_missing=false,
            path=joinpath(z_store.folder, RESULTS, "relative_taxa_cover"), chunks=((result_dims[1], 6)..., 1),
            attrs=Dict(
                :structure => string.(ADRIA.metrics.relative_taxa_cover.dims)
            ),
            compressor=compressor))
    push!(met_names, :relative_taxa_cover)

    dhw_stats_store = store_env_summary(domain.dhw_scens, "dhw_scenario", joinpath(z_store.folder, ENV_STATS, "dhw"), domain.RCP, compressor)
    wave_stats_store = store_env_summary(domain.wave_scens, "wave_scenario", joinpath(z_store.folder, ENV_STATS, "wave"), domain.RCP, compressor)

    # Group all data stores
    stores = [stores..., dhw_stats_store, wave_stats_store, setup_logs(z_store, unique_sites(domain), nrow(param_df), tf, n_sites)...]

    return domain, (; zip((met_names..., :dhw_stats, :wave_stats, :site_ranks, :seed_log, :fog_log, :shade_log,), stores)...)
end

"""
    _recreate_stats_from_store(zarr_store_path::String)::Dict{String, AbstractArray}

Recreate data structure holding RCP summary statistics from Zarr store.
"""
function _recreate_stats_from_store(zarr_store_path::String)::Dict{String,AbstractArray}
    rcp_dirs = filter(d -> isdir(joinpath(zarr_store_path, d)), readdir(zarr_store_path))
    rcp_stat_dirs = joinpath.(zarr_store_path, rcp_dirs)

    stat_d = Dict{String,AbstractArray}()
    for (i, sd) in enumerate(rcp_stat_dirs)
        store = zopen(sd, fill_as_missing=false)

        dims = store.attrs["structure"]
        row_names = string.(store.attrs["rows"])
        col_names = string.(store.attrs["cols"])
        stat_set = NamedDimsArray(store[:, :]; zip(Symbol.(dims), [row_names, col_names])...)

        stat_d[rcp_dirs[i]] = stat_set
    end

    return stat_d
end


"""
    load_results(result_loc::String)::ResultSet
    load_results(domain::Domain)::ResultSet

Create interface to a given Zarr result set.
"""
function load_results(result_loc::String)::ResultSet
    # Read in results
    local raw_set
    try
        raw_set = zopen(joinpath(result_loc, RESULTS), fill_as_missing=false)
    catch err
        if !occursin("ArgumentError", sprint(showerror, err))
            rethrow(err)
        end
    end

    # Read in logs
    log_set = zopen(joinpath(result_loc, LOG_GRP), fill_as_missing=false)
    input_set = zopen(joinpath(result_loc, INPUTS), fill_as_missing=false)

    dhw_stat_set = _recreate_stats_from_store(joinpath(result_loc, ENV_STATS, "dhw"))
    wave_stat_set = _recreate_stats_from_store(joinpath(result_loc, ENV_STATS, "wave"))

    result_loc = replace(result_loc, "\\" => "/")
    if endswith(result_loc, "/")
        result_loc = result_loc[1:end-1]
    end

    # Spatial data
    site_data = GeoDataFrames.read(joinpath(result_loc, SITE_DATA, input_set.attrs["name"] * ".gpkg"))
    sort!(site_data, [Symbol(input_set.attrs["unique_site_id_col"])])

    # Model specification
    model_spec = CSV.read(joinpath(result_loc, MODEL_SPEC, "model_spec.csv"), DataFrame; comment="#")

    r_vers_id = input_set.attrs["ADRIA_VERSION"]
    t_vers_id = "v" * string(PkgVersion.Version(@__MODULE__))

    if r_vers_id != t_vers_id
        msg = """Results were produced with a different version of ADRIA ($(r_vers_id)). The version of ADRIA in use is $(t_vers_id).\n
        Errors may occur when analyzing data."""

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
        input_set.attrs["unique_site_id_col"],
        input_set.attrs["init_coral_cover_file"],
        input_set.attrs["connectivity_file"],
        input_set.attrs["DHW_file"],
        input_set.attrs["wave_file"],
        input_set.attrs["timeframe"]
    )

    outcomes = Dict{Symbol,NamedDimsArray}()
    subdirs = filter(isdir, readdir(joinpath(result_loc, RESULTS), join=true))
    for sd in subdirs
        if !(occursin(LOG_GRP, sd)) && !(occursin(INPUTS, sd))
            res = zopen(sd, fill_as_missing=false)
            sz = size(res)

            # Construct dimension names and metadata
            st = []
            for (i, s) in enumerate(res.attrs["structure"])
                if s == "timesteps"
                    push!(st, input_set.attrs["timeframe"])
                elseif s == "sites"
                    push!(st, res.attrs["unique_site_ids"])
                else
                    push!(st, 1:sz[i])
                end
            end

            try
                outcomes[Symbol(basename(sd))] = NamedDimsArray(res; zip(Symbol.(res.attrs["structure"]), st)...)
            catch err
                if err isa ArgumentError
                    @warn """Unable to resolve result structure, reverting to numbered keys.
                    For group $(sd)
                    Got: $(size(res))
                    Structure: $(res.attrs["structure"])
                    Generated: $(Array([i[1] for i in size.(st)]))
                    """
                    outcomes[Symbol(basename(sd))] = NamedDimsArray(res; zip(Symbol.(res.attrs["structure"]), [1:s for s in size(res)])...)
                else
                    rethrow(err)
                end
            end
        end
    end

    return ResultSet(input_set, env_layer_md, inputs_used, outcomes, log_set, dhw_stat_set, wave_stat_set, site_data, model_spec)
end
function load_results(domain::Domain)::ResultSet
    return load_results(result_location(domain))
end

"""
    result_location(d::Domain)::String

Generate path to the data store of results for the given Domain.
"""
function result_location(d::Domain)::String
    return joinpath(ENV["ADRIA_OUTPUT_DIR"], "$(d.name)__RCPs$(d.RCP)__$(d.scenario_invoke_time)")
end
