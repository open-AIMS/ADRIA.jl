using Dates: now

using PkgVersion, Zarr

using Setfield: @set!
using DataFrames: DataFrame
using ADRIA: Domain, EnvLayer

const CONNECTIVITY = "connectivity"
const ENV_STATS = "env_stats"
const INPUTS = "inputs"
const LOG_GRP = "logs"
const MODEL_SPEC = "model_spec"
const RESULTS = "results"
const SPATIAL_DATA = "spatial"

struct ResultSet{T1,T2,A,B,C,D,G,D1,D2,D3,DF}
    name::String
    RCP::String
    invoke_time::String
    ADRIA_VERSION::String

    site_ids::T1
    site_area::Vector{Float64}
    site_max_coral_cover::Vector{Float64}
    site_centroids::T2
    env_layer_md::EnvLayer
    dhw_stats::D
    wave_stats::D
    connectivity_data::D
    site_data::G

    inputs::G
    sim_constants::D1
    model_spec::DF

    # raw::AbstractArray
    outcomes::D2
    ranks::A
    seed_log::B  # Values stored in m^2
    fog_log::C   # Reduction in bleaching mortality (0.0 - 1.0)
    shade_log::C # Reduction in bleaching mortality (0.0 - 1.0)
    coral_dhw_tol_log::D3
end

function ResultSet(
    input_set::AbstractArray,
    env_layer_md::EnvLayer,
    inputs_used::DataFrame,
    outcomes::Dict,
    log_set::Zarr.ZGroup,
    dhw_stats_set::Dict,
    wave_stats_set::Dict,
    conn_data::Dict,
    site_data::DataFrame,
    model_spec::DataFrame,
)::ResultSet
    rcp = "RCP" in keys(input_set.attrs) ? input_set.attrs["RCP"] : input_set.attrs["rcp"]
    return ResultSet(input_set.attrs["name"],
        string(rcp),
        input_set.attrs["invoke_time"],
        input_set.attrs["ADRIA_VERSION"],
        input_set.attrs["site_ids"],
        convert.(Float64, input_set.attrs["site_area"]),
        convert.(Float64, input_set.attrs["site_max_coral_cover"]),
        input_set.attrs["site_centroids"],
        env_layer_md,
        dhw_stats_set,
        wave_stats_set,
        conn_data,
        site_data,
        inputs_used,
        input_set.attrs["sim_constants"],
        model_spec,
        outcomes,
        _rankings_data(log_set["rankings"]),
        DataCube(log_set["seed"], Symbol.(Tuple(log_set["seed"].attrs["structure"]))),
        DataCube(log_set["fog"], Symbol.(Tuple(log_set["fog"].attrs["structure"]))),
        DataCube(log_set["shade"], Symbol.(Tuple(log_set["shade"].attrs["structure"]))),
        DataCube(log_set["coral_dhw_log"], Symbol.(Tuple(log_set["coral_dhw_log"].attrs["structure"]))),
    )
end

function _rankings_data(rankings_set::ZArray{T})::YAXArray{T} where {T}
    ax_names = Symbol.(Tuple(rankings_set.attrs["structure"]))
    ax_labels::Vector{Union{UnitRange{Int64},Vector{Symbol}}} = range.([1], size(rankings_set))

    # Replace intervention
    intervention_idx = findfirst(x -> x == :intervention, ax_names)
    ax_labels[intervention_idx] = interventions()

    return DataCube(rankings_set; NamedTuple{ax_names}(ax_labels)...)
end

"""
    _copy_env_data(src::String, dst::String, subdir::String)::Nothing

Helper function to copy environmental data layer statistics from data store.
"""
function _copy_env_data(src::String, dst::String, folder_name::String, subdir=""::String)::Nothing
    src_dir = joinpath(src, folder_name, subdir)
    dst_dir = joinpath(dst, folder_name, subdir)
    mkpath(dst_dir)
    src_ds = filter(d -> isdir(joinpath(src_dir, d)), readdir(src_dir))
    for ds in src_ds
        cp(joinpath(src_dir, ds), joinpath(dst_dir, ds); force=true)
    end

    return nothing
end

"""
    combine(result_sets...)::ResultSet
    combine_results(result_set_locs::Array{String})::ResultSet

Combine arbitrary number of ADRIA result sets into a single data store.

Note: Results are stored in Zarr format. Combining data sets can be
  expected to double total disk space requirement.
"""
function combine_results(result_sets...)::ResultSet
    # Make sure results are all for the same domain
    # @assert length(Set([rs.name for rs in result_sets])) == 1

    # Ensure all sim constants are identical
    @assert all([
        result_sets[i].sim_constants == result_sets[i+1].sim_constants
        for i in 1:(length(result_sets)-1)
    ])

    # Ensure all result sets were from the same version of ADRIA
    if length(Set([rs.ADRIA_VERSION for rs in result_sets])) != 1
        @warn "Results were created with different versions of ADRIA so errors may occur!"
        @warn "Results from model runs < 0.9 are no longer compatible!"
    end

    rs1 = result_sets[1]
    canonical_name = rs1.name
    combined_time = replace(string(now()), "T" => "_", ":" => "_", "." => "_")

    rcps = join(unique([rs.RCP for rs in result_sets]), "_")

    # Create new zarr store
    new_loc = joinpath(
        ENV["ADRIA_OUTPUT_DIR"], "$(rs1.name)__RCPs$(rcps)__$(combined_time)"
    )
    z_store = DirectoryStore(new_loc)

    # Get input sizes
    n_scenarios = sum(map((rs) -> size(rs.inputs, 1), result_sets))

    envlayer = rs1.env_layer_md
    env_md = EnvLayer(
        envlayer.dpkg_path,
        envlayer.site_data_fn,
        envlayer.site_id_col,
        envlayer.cluster_id_col,
        envlayer.init_coral_cov_fn,
        envlayer.connectivity_fn,
        dirname(envlayer.DHW_fn),
        dirname(envlayer.wave_fn),
        envlayer.timeframe,
    )

    all_inputs = reduce(vcat, [getfield(rs, :inputs) for rs in result_sets])
    input_dims = size(all_inputs)
    attrs = scenario_attributes(
        canonical_name,
        rcps,
        names(all_inputs),
        combined_time,
        env_md,
        rs1.sim_constants,
        rs1.site_ids,
        rs1.site_area,
        rs1.site_max_coral_cover,
        rs1.site_centroids,
    )

    # Copy site data into result set
    mkdir(joinpath(new_loc, SPATIAL_DATA))
    cp(
        attrs[:site_data_file],
        joinpath(new_loc, SPATIAL_DATA, basename(attrs[:site_data_file]));
        force=true,
    )

    # Store copy of model specification as CSV
    mkdir(joinpath(new_loc, MODEL_SPEC))
    CSV.write(joinpath(new_loc, MODEL_SPEC, "model_spec.csv"), rs1.model_spec)

    input_loc::String = joinpath(z_store.folder, INPUTS)
    input_set = zcreate(
        Float64,
        input_dims...;
        fill_value=-9999.0,
        fill_as_missing=false,
        path=input_loc,
        chunks=(1, input_dims[2]),
        attrs=attrs,
    )

    # Store post-processed table of input parameters.
    input_set[:, :] = Matrix(all_inputs)
    logs = (;
        zip([:ranks, :seed_log, :fog_log, :shade_log, :coral_dhw_tol_log],
            setup_logs(
                z_store,
                rs1.site_ids,
                nrow(all_inputs),
                size(rs1.seed_log, :timesteps),
                size(rs1.seed_log, :sites),
            ))...
    )

    # Copy logs over
    for log in keys(logs)
        scen_id = 1
        for rs in result_sets
            s_log = getfield(rs, log)
            n_log = getfield(logs, log)
            rs_scen_len = size(s_log, :scenarios)

            try
                n_log[:, :, scen_id:(scen_id+(rs_scen_len-1))] .= s_log
            catch
                n_log[:, :, :, scen_id:(scen_id+(rs_scen_len-1))] .= s_log
            end

            scen_id = scen_id + rs_scen_len
        end
    end

    compressor = Zarr.BloscCompressor(; cname="zstd", clevel=4, shuffle=true)
    metrics = keys(rs1.outcomes)
    for m_name in metrics
        m_dim_names = axes_names(rs1.outcomes[m_name])
        dim_struct = Dict{Symbol,Any}(
            :structure => m_dim_names
        )
        if :sites in m_dim_names
            dim_struct[:unique_site_ids] = rs1.site_ids
        end

        result_dims = (size(rs1.outcomes[m_name])[1:(end-1)]..., n_scenarios)
        m_store = zcreate(
            Float32,
            result_dims...;
            fill_value=nothing,
            fill_as_missing=false,
            path=joinpath(
                z_store.folder,
                RESULTS,
                string(m_name)),
            chunks=(result_dims[1:(end-1)]..., 1),
            attrs=dim_struct,
            compressor=compressor,
        )

        # Copy results over
        scen_id = 1
        for rs in result_sets
            rs_scen_len = size(rs.outcomes[m_name], :scenarios)
            try
                m_store[:, :, scen_id:(scen_id+(rs_scen_len-1))] .= rs.outcomes[m_name]
            catch
                m_store[:, :, :, scen_id:(scen_id+(rs_scen_len-1))] .= rs.outcomes[m_name]
            end

            scen_id = scen_id + rs_scen_len
        end
    end

    # Copy env stats
    mkdir(joinpath(new_loc, ENV_STATS))
    for rs in result_sets
        loc::String = rs.env_layer_md.dpkg_path
        _copy_env_data(loc, new_loc, ENV_STATS, "dhw")
        _copy_env_data(loc, new_loc, ENV_STATS, "wave")
        _copy_env_data(loc, new_loc, CONNECTIVITY)
    end

    return load_results(z_store.folder)
end
function combine_results(result_set_locs::Array{String})::ResultSet
    return combine_results(load_results.(result_set_locs)...)
end

"""
    env_stats(rs::ResultSet, s_name::String, rcp::String)
    env_stats(rs::ResultSet, s_name::String, rcp::String, scenario::Int)
    env_stats(rs::ResultSet, s_name::String, stat::String, rcp::String, scenario::Int)

Extract statistics for a given environmental layer ("DHW" or "wave")
"""
function env_stats(rs::ResultSet, s_name::String, rcp::String)
    return getfield(rs, Symbol("$(s_name)_stats"))[rcp]
end
function env_stats(rs::ResultSet, s_name::String, rcp::String, scenario::Int)
    return getfield(rs, Symbol("$(s_name)_stats"))[rcp][:, scenario]
end
function env_stats(rs::ResultSet, s_name::String, stat::String, rcp::String, scenario::Int)
    return getfield(rs, Symbol("$(s_name)_stats"))[rcp][stat, scenario]
end

"""
    store_name(rs::ResultSet)::String

Get name of result set.
"""
function store_name(rs::ResultSet)::String
    return "$(rs.name)__RCPs$(rs.RCP)__$(rs.invoke_time)"
end

"""
    result_location(rs::ResultSet)::String

Get location of result set.
"""
function result_location(rs::ResultSet)::String
    store = ""
    try
        store = joinpath(ENV["ADRIA_OUTPUT_DIR"], store_name(rs))
    catch err
        if isa(err, KeyError)
            @warn "Output directory not yet set. Displaying result directory instead."
            store = store_name(rs)
        else
            rethrow(err)
        end
    end

    return replace(store, "\\" => "/")
end

"""
    select(r::ResultSet, op::String)

Hacky scenario filtering - to be replaced with more robust approach.

Only supports filtering by single attribute.
Should be expanded to support filtering metric results too.

# Examples
```julia
select(result, "guided .> 0.0")

# Above expands to:
# result.inputs.guided .> 0.0
```
"""
function select(r::ResultSet, op::String)
    scens = r.inputs

    col, qry = split(op, " "; limit=2)
    col = Symbol(col)

    df_ss = getproperty(scens, col)

    return eval(Meta.parse("$df_ss $qry"))
end

"""
    timesteps(rs::ResultSet)

Retrieve the time steps represented in the result set.

# Arguments
- `rs` : ResultSet
"""
function timesteps(rs::ResultSet)
    return rs.env_layer_md.timeframe
end

"""
    site_k_area(rs::ResultSet)::Vector{Float64}

Extract vector of a location's coral carrying capacity in terms of absolute area.
"""
function site_k_area(rs::ResultSet)::Vector{Float64}
    return rs.site_max_coral_cover .* rs.site_area
end

"""
    site_k(rs::ResultSet)::Vector{Float64}

Extract vector of a location's coral carrying capacity in as a proportion relative to the
location's total area.
"""
function site_k(rs::ResultSet)::Vector{Float64}
    return rs.site_max_coral_cover
end

"""
    site_area(rs::ResultSet)::Vector{Float64}

Extract vector of a location's total area in its areal unit (m², km², etc).
"""
function site_area(rs::ResultSet)::Vector{Float64}
    return rs.site_area
end

"""
    n_locations(rs::ResultSet)::Int64

Retrieve the number of locations represented in the result set.

# Arguments
- `rs` : ResultSet
"""
function n_locations(rs::ResultSet)::Int64
    return size(rs.site_ids, 1)
end

"""
    n_scenarios(rs::ResultSet)::Int64

Number of scenarios stored in the given result set.

# Arguments
- `rs` : ResultSet
"""
function n_scenarios(rs::ResultSet)::Int64
    return size(rs.inputs, 1)
end

"""
    component_params(spec::DataFrame, component::Type)::DataFrame

Extract parameters for a specific model component from exported model specification.
"""
function component_params(rs::ResultSet, component::T)::DataFrame where {T}
    spec = rs.model_spec
    return spec[spec.component.==string(component), :]
end
function component_params(rs::ResultSet, components::Vector{T})::DataFrame where {T}
    spec = rs.model_spec
    return spec[spec.component.∈[replace.(string.(components), "ADRIA." => "")], :]
end

"""
    model_spec(rs::ResultSet)::DataFrame

Extract model specification from Result Set.
"""
function model_spec(rs::ResultSet)::DataFrame
    return rs.model_spec
end

function Base.show(io::IO, mime::MIME"text/plain", rs::ResultSet)
    vers_id = rs.ADRIA_VERSION

    tf, sites, scens = size(rs.outcomes[:relative_cover])
    # Species/size groups represented: $(species)

    rcps = join(split(rs.RCP, "_"), ", ")

    println("""
    Domain: $(rs.name)

    Run with ADRIA $(vers_id) on $(rs.invoke_time)
    Results stored at: $(result_location(rs))

    RCP(s) represented: $(rcps)
    Scenarios run: $(scens)
    Number of sites: $(sites)
    Timesteps: $(tf)

    Input layers
    ------------""")

    for fn in fieldnames(typeof(rs.env_layer_md))
        if fn == :timeframe
            tf = getfield(rs.env_layer_md, fn)
            println("$(fn) : $(tf[1]) - $(tf[end])")
            continue
        end

        println("$(fn) : $(getfield(rs.env_layer_md, fn))")
    end
end
