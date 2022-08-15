"""
    centroids(df::DataFrame)

Extract and return long/lat from a GeoDataFrame.

# Arguments
df : GeoDataFrame

# Returns
Array of tuples (x, y), where x and y relate to long and lat respectively.
"""
function centroids(df::DataFrame)::Array
    site_centroids = AG.centroid.(df.geometry)
    return collect(zip(AG.getx.(site_centroids, 0), AG.gety.(site_centroids, 0)))
end


"""
    scenario_attributes(domain::Domain, param_df::DataFrame)

Generate dictionary of scenario attributes.
"""
function scenario_attributes(domain::Domain, param_df::DataFrame)
    attrs::Dict = Dict(
        :name => domain.name,
        :rcp => domain.rcp,
        :columns => names(param_df),
        :invoke_time => domain.scenario_invoke_time,
        :ADRIA_VERSION => "v" * string(PkgVersion.Version(@__MODULE__)),

        :site_data_file => domain.env_layer_md.site_data_fn,
        :site_id_col => domain.env_layer_md.site_id_col,
        :unique_site_id_col => domain.env_layer_md.unique_site_id_col,
        :init_coral_cover_file => domain.env_layer_md.init_coral_cov_fn,
        :connectivity_file => domain.env_layer_md.connectivity_fn,
        :DHW_file => domain.env_layer_md.DHW_fn,
        :wave_file => domain.env_layer_md.wave_fn,
        :sim_constants => Dict(fn=>getfield(domain.sim_constants, fn) for fn ∈ fieldnames(typeof(domain.sim_constants))),

        :site_ids => unique_sites(domain),
        :site_area => domain.site_data.area,
        :site_max_coral_cover => domain.site_data.k
        :site_centroids => centroids(domain.site_data)
    )

    return attrs
end


function setup_logs(z_store, domain, param_df, tf, n_sites, reps)
    # Set up logs for site ranks, seed/fog log
    zgroup(z_store, LOG_GRP)
    log_fn::String = joinpath(z_store.folder, LOG_GRP)

    # Store ranked sites
    n_scens::Int64 = nrow(param_df)
    rank_dims::Tuple{Int64,Int64,Int64,Int64} = (tf, n_sites, 2, n_scens)  # sites, site id and rank, reps, no. scenarios
    fog_dims::Tuple{Int64,Int64,Int64,Int64} = (tf, n_sites, reps, n_scens)  # timeframe, sites, reps, no. scenarios

    # tf, no. intervention sites, site id and rank, no. scenarios
    seed_dims::Tuple{Int64,Int64,Int64,Int64,Int64} = (tf, 2, n_sites, reps, n_scens)

    attrs = Dict(
        :structure=> ("timesteps", "sites", "reps", "scenarios"),
        :unique_site_ids=>unique_sites(domain),
    )
    ranks = zcreate(Float32, rank_dims...; name="rankings", fill_value=nothing, fill_as_missing=false, path=log_fn, chunks=(rank_dims[1:3]..., 1), attrs=attrs)

    attrs = Dict(
        :structure=> ("timesteps", "coral_id", "sites", "reps", "scenarios"),
        :unique_site_ids=>unique_sites(domain),
    )
    seed_log = zcreate(Float32, seed_dims...; name="seed", fill_value=nothing, fill_as_missing=false, path=log_fn, chunks=(seed_dims[1:4]..., 1), attrs=attrs)

    attrs = Dict(
        :structure=> ("timesteps", "sites", "reps", "scenarios"),
        :unique_site_ids=>unique_sites(domain),
    )
    fog_log = zcreate(Float32, fog_dims...; name="fog", fill_value=nothing, fill_as_missing=false, path=log_fn, chunks=(fog_dims[1:3]..., 1), attrs=attrs)
    shade_log = zcreate(Float32, fog_dims...; name="shade", fill_value=nothing, fill_as_missing=false, path=log_fn, chunks=(fog_dims[1:3]..., 1), attrs=attrs)

    return ranks, seed_log, fog_log, shade_log
end


"""
    setup_result_store!(domain::Domain, param_df::DataFrame, reps::Int64=0; metrics::Array=[])

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
└───inputs
```

`inputs` store includes domain specification metadata including what connectivity/DHW/wave data was used.

# Notes
- `domain` is replaced with an identical copy with an updated scenario invoke time.
- -9999.0 is used as an arbitrary fill value.

# Inputs
- domain : ADRIA scenario domain
- param_df : ADRIA scenario specification
- reps : Number of environmental scenario repeats

# Returns
domain, (relative_cover, relative_shelter_volume, absolute_shelter_volume, site_ranks, seed_log, fog_log, shade_log)

"""
function setup_result_store!(domain::Domain, param_df::DataFrame, reps::Int64)::Tuple
    @set! domain.scenario_invoke_time = replace(string(now()), "T"=>"_", ":"=>"_", "."=>"_")
    log_location::String = joinpath(ENV["ADRIA_OUTPUT_DIR"], "$(domain.name)__$(domain.rcp)__$(domain.scenario_invoke_time)")

    z_store = DirectoryStore(log_location)

    # Store copy of inputs
    input_loc::String = joinpath(z_store.folder, INPUTS)
    input_dims::Tuple{Int64,Int64} = size(param_df)
    attrs::Dict = scenario_attributes(domain, param_df)

    inputs = zcreate(Float64, input_dims...; fill_value=-9999.0, fill_as_missing=false, path=input_loc, chunks=input_dims, attrs=attrs)

    # Store post-processed table of input parameters.
    integer_params = findall(domain.model[:ptype] .== "integer")
    map_to_discrete!(param_df[:, integer_params], getindex.(domain.model[:bounds], 2)[integer_params])
    inputs[:, :] = Matrix(param_df)

    # Set up stores for each metric
    tf, n_sites = domain.sim_constants.tf::Int64, domain.coral_growth.n_sites::Int64

    function dim_lengths(metric_structure)
        dl = []
        for d in metric_structure
            if d == "timesteps"
                append!(dl, tf)
            elseif d == "species"
                append!(dl, domain.coral_growth.n_species)
            elseif d == "sites"
                append!(dl, n_sites)
            elseif d == "reps"
                append!(dl, reps)
            elseif d == "scenarios"
                append!(dl, nrow(param_df))
            end
        end

        return (dl..., )
    end

    compressor = Zarr.BloscCompressor(cname="zstd", clevel=2, shuffle=true)

    met_names = (:relative_cover, :relative_shelter_volume, :absolute_shelter_volume)
    dim_struct = Dict(
        :structure => string.((:timesteps, :sites, :reps, :scenarios)),
        :unique_site_ids => unique_sites(domain)
    )
    result_dims::Tuple{Int64,Int64,Int64,Int64} = dim_lengths(dim_struct[:structure])

    # Create stores for each metric
    stores = [
        zcreate(Float32, result_dims...;
                fill_value=nothing, fill_as_missing=false,
                path=joinpath(z_store.folder, "results", string(m_name)), chunks=(result_dims[1:end-1]..., 1),
                attrs=dim_struct,
                compressor=compressor)
        for m_name in met_names
    ]

    # Set up logs for site ranks, seed/fog log
    stores = [stores..., setup_logs(z_store, domain, param_df, tf, n_sites, reps)...]

    # NamedTuple{(Symbol.(metrics)..., :site_ranks, :seed_log, :fog_log, :shade_log)}(stores)
    return domain, (; zip((met_names..., :site_ranks, :seed_log, :fog_log, :shade_log), stores)...)
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

    outcomes = Dict{Symbol, NamedDimsArray}()
    subdirs = filter(isdir, readdir(joinpath(result_loc, RESULTS), join=true))
    for sd in subdirs
        if !(occursin(LOG_GRP, sd)) && !(occursin(INPUTS, sd))
            res = zopen(sd, fill_as_missing=false)
            outcomes[Symbol(basename(sd))] = NamedDimsArray{Symbol.(Tuple(res.attrs["structure"]))}(res)
        end
    end

    log_set = zopen(joinpath(result_loc, LOG_GRP), fill_as_missing=false)
    input_set = zopen(joinpath(result_loc, INPUTS), fill_as_missing=false)

    r_vers_id = input_set.attrs["ADRIA_VERSION"]
    t_vers_id = "v" * string(PkgVersion.Version(@__MODULE__))

    if r_vers_id != t_vers_id
        msg = """Results were produced with an older version of ADRIA ($(r_vers_id)). The installed version of ADRIA is $(t_vers_id).\n
        Errors may occur when analyzing data."""

        @warn msg
    end

    input_cols::Array{String} = input_set.attrs["columns"]
    inputs_used::DataFrame = DataFrame(input_set[:, :], input_cols)

    env_layer_md::EnvLayer = EnvLayer(
        input_set.attrs["site_data_file"],
        input_set.attrs["site_id_col"],
        input_set.attrs["unique_site_id_col"],
        input_set.attrs["init_coral_cover_file"],
        input_set.attrs["connectivity_file"],
        input_set.attrs["DHW_file"],
        input_set.attrs["wave_file"]
    )

    return ResultSet(input_set, env_layer_md, inputs_used, outcomes, log_set)
end
function load_results(domain::Domain)::ResultSet
    log_location = joinpath(ENV["ADRIA_OUTPUT_DIR"], "$(domain.name)__$(domain.rcp)__$(domain.scenario_invoke_time)")
    return load_results(log_location)
end
