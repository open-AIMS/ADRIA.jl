function get_geometry(df::DataFrame)
    if "geometry" in names(df)
        return df.geometry
    elseif "geom" in names(df)
        return df.geom
    end

    error("No geometry data found")
end

"""
    centroids(df::DataFrame)

Extract and return long/lat from a GeoDataFrame.

# Arguments
df : GeoDataFrame

# Returns
Array of tuples (x, y), where x and y relate to long and lat respectively.
"""
function centroids(df::DataFrame)::Array
    site_centroids = AG.centroid.(get_geometry(df::DataFrame))
    return collect(zip(AG.getx.(site_centroids, 0), AG.gety.(site_centroids, 0)))
end

"""
    calculate_wave_dhw_summary(dhw_df::DataFrame,wave_df::DataFrame,param_df::DataFrame)

Retrieve summary statistics matrices from DataFrames of dhws and waves.

# Arguments
dhw_df : DataFrame of dhw data.
waves_df : DataFrame of dhw data.
param_df : DataFrame of scenario parameters.

# Returns
dhw_stats, wave_stats: two 2*N matrices where the first row is mean over time and sites and 2nd row is the std over time and mean of sites.
N is the number of unique dhw/wave scenarios being run.
"""
function store_wave_dhw_summary(data_cube::NamedArray, param_df::DataFrame, type::String, file_loc::String, compressor::Zarr.Compressor)
    unique_scens = unique(param_df[:, [Symbol(type)]])

    scens = data_cube[:, :, convert.(Int64, unique_scens[:, Symbol(type)])]

    scens_stats = vcat(dropdims(mean(scens, dims=[1, 2]), dims=1),
        dropdims(std(mean(scens, dims=2), dims=1), dims=1))

    stats_store = zcreate(Float32, (2, size(unique_scens)[1])...;
        fill_value=nothing, fill_as_missing=false,
        path=file_loc,
        attrs=Dict(:structure => ("stat", type),
            :rows => ["mean", "std"],
            :cols => unique_scens[:, Symbol(type)]),
        compressor=compressor)

    stats_store[:, :] .= scens_stats

    return stats_store
end

"""
    scenario_attributes(name, RCP, input_cols, invoke_time, env_layer, sim_constants, unique_sites, area, k)
    scenario_attributes(domain::Domain, param_df::DataFrame)

Generate dictionary of scenario attributes.
"""
function scenario_attributes(name, RCP, input_cols, invoke_time, env_layer, sim_constants, unique_sites, area, k, centroids)::Dict
    attrs::Dict = Dict(
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
    setup_result_store!(domain::Domain, param_df::DataFrame; metrics::Array=[])

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

# Inputs
- domain : ADRIA scenario domain
- param_df : ADRIA scenario specification

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

    met_names = (:relative_cover, :relative_shelter_volume, :absolute_shelter_volume)
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
    dhw_stats_store = store_wave_dhw_summary(domain.dhw_scens, param_df, "dhw_scenario", joinpath(z_store.folder, RESULTS, DHW_STATS), compressor)
    wave_stats_store = store_wave_dhw_summary(domain.wave_scens, param_df, "wave_scenario", joinpath(z_store.folder, RESULTS, WAVE_STATS), compressor)

    # Set up logs for site ranks, seed/fog log
    stores = [stores..., dhw_stats_store, wave_stats_store, setup_logs(z_store, unique_sites(domain), nrow(param_df), tf, n_sites)...]

    # NamedTuple{(Symbol.(metrics)..., :site_ranks, :seed_log, :fog_log, :shade_log)}(stores)
    return domain, (; zip((met_names..., :dhw_stats, :wave_stats, :site_ranks, :seed_log, :fog_log, :shade_log,), stores)...)
end


"""
    load_results(result_loc::String)::ResultSet
    load_results(domain::Domain)::ResultSet

Create interface to a given Zarr result set.
"""
function load_results(result_loc::String)::ResultSet
    raw_set = nothing
    try
        raw_set = zopen(joinpath(result_loc, RESULTS), fill_as_missing=false)
    catch err
        if !occursin("ArgumentError", sprint(showerror, err))
            rethrow(err)
        end
    end

    outcomes = Dict{Symbol,NamedDimsArray}()
    subdirs = filter(isdir, readdir(joinpath(result_loc, RESULTS), join=true))
    for sd in subdirs
        if !(occursin(LOG_GRP, sd)) && !(occursin(INPUTS, sd))
            res = zopen(sd, fill_as_missing=false)
            outcomes[Symbol(basename(sd))] = NamedDimsArray{Symbol.(Tuple(res.attrs["structure"]))}(res)
        end
    end

    log_set = zopen(joinpath(result_loc, LOG_GRP), fill_as_missing=false)
    input_set = zopen(joinpath(result_loc, INPUTS), fill_as_missing=false)
    dhw_stat_set = zopen(joinpath(result_loc, RESULTS, DHW_STATS), fill_as_missing=false)
    wave_stat_set = zopen(joinpath(result_loc, RESULTS, WAVE_STATS), fill_as_missing=false)

    result_loc = replace(result_loc, "\\" => "/")
    if endswith(result_loc, "/")
        result_loc = result_loc[1:end-1]
    end

    site_data = GeoDataFrames.read(joinpath(result_loc, SITE_DATA, input_set.attrs["name"] * ".gpkg"))
    model_spec = CSV.read(joinpath(result_loc, MODEL_SPEC, "model_spec.csv"), DataFrame; comment="#")

    r_vers_id = input_set.attrs["ADRIA_VERSION"]
    t_vers_id = "v" * string(PkgVersion.Version(@__MODULE__))

    if r_vers_id != t_vers_id
        msg = """Results were produced with a different version of ADRIA ($(r_vers_id)). The version of ADRIA in use is $(t_vers_id).\n
        Errors may occur when analyzing data."""

        @warn msg
    end

    input_cols::Array{String} = input_set.attrs["columns"]
    inputs_used::DataFrame = DataFrame(input_set[:, :], input_cols)

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
