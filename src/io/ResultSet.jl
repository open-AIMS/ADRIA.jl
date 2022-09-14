import Dates: now

using PkgVersion, Zarr, NamedDims

import Setfield: @set!
import DataFrames: DataFrame
import ADRIA: Domain, EnvLayer


const RESULTS = "results"
const LOG_GRP = "logs"
const INPUTS = "inputs"
const SITE_DATA = "site_data"
const MODEL_SPEC = "model_spec"


struct ResultSet{S, T1, T2, F, A, B, C, G}
    name::S
    RCP::S
    invoke_time::S
    ADRIA_VERSION::S

    site_ids::T1
    site_area::F
    site_max_coral_cover::F
    site_centroids::T2
    env_layer_md::EnvLayer
    site_data::G

    inputs::G
    sim_constants::Dict
    model_spec::DataFrame

    # raw::AbstractArray
    outcomes::Dict
    ranks::A
    seed_log::B  # Values stored in m^2
    fog_log::C   # Reduction in bleaching mortality (0.0 - 1.0)
    shade_log::C # Reduction in bleaching mortality (0.0 - 1.0)
end


function ResultSet(input_set::Zarr.ZArray, env_layer_md::EnvLayer, 
    inputs_used::DataFrame, outcomes::Dict, log_set::Zarr.ZGroup, 
    site_data::DataFrame, model_spec::DataFrame)::ResultSet

    rcp = "RCP" in keys(input_set.attrs) ? input_set.attrs["RCP"] : input_set.attrs["rcp"]
    ResultSet(input_set.attrs["name"],
              string(rcp),
              input_set.attrs["invoke_time"],
              input_set.attrs["ADRIA_VERSION"],
              input_set.attrs["site_ids"],
              convert.(Float64, input_set.attrs["site_area"]),
              convert.(Float64, input_set.attrs["site_max_coral_cover"]),
              input_set.attrs["site_centroids"],
              env_layer_md,
              site_data,
              inputs_used,
              input_set.attrs["sim_constants"],
              model_spec,
              outcomes,
              NamedDimsArray{Symbol.(Tuple(log_set["rankings"].attrs["structure"]))}(log_set["rankings"]),
              NamedDimsArray{Symbol.(Tuple(log_set["seed"].attrs["structure"]))}(log_set["seed"]),
              NamedDimsArray{Symbol.(Tuple(log_set["fog"].attrs["structure"]))}(log_set["fog"]),
              NamedDimsArray{Symbol.(Tuple(log_set["shade"].attrs["structure"]))}(log_set["shade"]))
end


"""
    combine(result_sets...)

Combine arbitrary number of ADRIA result sets into a single data store.

Note: Results are stored in Zarr format. Combining data sets can be
      expected to double total disk space requirement.
"""
function combine_results(result_sets...)::ResultSet
    # Make sure results are all for the same domain
    # @assert length(Set([rs.name for rs in result_sets])) == 1

    # Ensure all sim constants are identical
    @assert all([result_sets[i].sim_constants == result_sets[i+1].sim_constants for i in 1:length(result_sets)-1])


    # Ensure all result sets were from the same version of ADRIA
    if length(Set([rs.ADRIA_VERSION for rs in result_sets])) != 1
        @warn "Results were created with different versions of ADRIA so errors may occur!"
    end

    rs1 = result_sets[1]
    canonical_name = rs1.name
    combined_time = replace(string(now()), "T"=>"_", ":"=>"_", "."=>"_")

    rcps = join(unique([rs.RCP for rs in result_sets]), "_")

    # Create new zarr store
    new_loc = joinpath(ENV["ADRIA_OUTPUT_DIR"], "$(rs1.name)__RCPs$(rcps)__$(combined_time)")
    z_store = DirectoryStore(new_loc)

    # Get input sizes
    n_scenarios = sum(map((rs) -> size(rs.inputs, 1), result_sets))

    envlayer = rs1.env_layer_md
    env_md = EnvLayer(envlayer.site_data_fn, envlayer.site_id_col, envlayer.unique_site_id_col,
                      envlayer.init_coral_cov_fn, envlayer.connectivity_fn,
                      dirname(envlayer.DHW_fn), dirname(envlayer.wave_fn), envlayer.timeframe)

    all_inputs = reduce(vcat, [getfield(rs, :inputs) for rs in result_sets])
    input_dims = size(all_inputs)
    attrs = scenario_attributes(canonical_name, rcps, names(all_inputs), combined_time, env_md, rs1.sim_constants,
                                rs1.site_ids, rs1.site_area, rs1.site_max_coral_cover, rs1.site_centroids)

    # Copy site data into result set
    mkdir(joinpath(new_loc, SITE_DATA))
    cp(attrs[:site_data_file], joinpath(new_loc, SITE_DATA, basename(attrs[:site_data_file])), force=true)

    # Store copy of model specification as CSV
    mkdir(joinpath(new_loc, MODEL_SPEC))
    CSV.write(joinpath(new_loc, MODEL_SPEC, "model_spec.csv"), rs1.model_spec)
    # model_spec(dom, joinpath(log_location, joinpath(new_loc, MODEL_SPEC, "model_spec.csv")))

    input_loc::String = joinpath(z_store.folder, INPUTS)
    # TODO: Fix issue - ERROR: UndefRefError: access to undefined reference
    input_set = zcreate(Float64, input_dims...; fill_value=-9999.0, fill_as_missing=false, path=input_loc, chunks=(1, input_dims[2]), attrs=attrs)

    # Store post-processed table of input parameters.
    input_set[:, :] = Matrix(all_inputs)

    logs = (; zip([:ranks, :seed_log, :fog_log, :shade_log],
        setup_logs(z_store, rs1.site_ids, nrow(all_inputs), size(rs1.seed_log, :timesteps), size(rs1.seed_log, :sites)))...)

    # Copy logs over
    for log in keys(logs)
        scen_id = 1
        for rs in result_sets
            s_log = getfield(rs, log)
            n_log = getfield(logs, log)
            rs_scen_len = size(s_log, :scenarios)

            try
                n_log[:, :, scen_id:scen_id+(rs_scen_len-1)] .= s_log
            catch
                n_log[:, :, :, scen_id:scen_id+(rs_scen_len-1)] .= s_log
            end

            scen_id = scen_id+rs_scen_len
        end
    end

    compressor = Zarr.BloscCompressor(cname="zstd", clevel=2, shuffle=true)
    metrics = keys(rs1.outcomes)
    result_dims = (size(rs1.outcomes[:relative_cover], :timesteps), size(rs1.outcomes[:relative_cover], :sites), n_scenarios)
    dim_struct = Dict(
        :structure => ("timesteps", "sites", "scenarios"), # string.((:timesteps, :sites, :scenarios)),
        :unique_site_ids => rs1.site_ids
    )
    for m_name in metrics
        m_store = zcreate(Float32, result_dims...;
            fill_value=nothing, fill_as_missing=false,
            path=joinpath(z_store.folder, "results", string(m_name)), chunks=(result_dims[1:end-1]..., 1),
            attrs=dim_struct,
            compressor=compressor)

        # Copy results over
        scen_id = 1
        for rs in result_sets
            rs_scen_len = size(rs.outcomes[m_name], :scenarios)
            try
                m_store[:, :, scen_id:scen_id+(rs_scen_len-1)] .= rs.outcomes[m_name]
            catch
                m_store[:, :, :, scen_id:scen_id+(rs_scen_len-1)] .= rs.outcomes[m_name]
            end

            scen_id = scen_id+rs_scen_len
        end
    end

    return load_results(z_store.folder)
end


"""
    store_name(r::ResultSet)::String

Get name of result set.
"""
function store_name(r::ResultSet)::String
    return "$(r.name)__RCPs$(r.RCP)__$(r.invoke_time)"
end


"""
    store_location(r::ResultSet)::String

Get location of result set.
"""
function store_location(rs::ResultSet)::String
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

    return replace(store, "\\"=>"/")
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

    col, qry = split(op, " ", limit=2)
    col = Symbol(col)

    df_ss = getproperty(scens, col)

    return eval(Meta.parse("$df_ss $qry"))
end


"""
    timesteps(rs::ResultSet)

Retrieve the time steps represented in the result set.

# Arguments
- rs : ResultSet
"""
function timesteps(rs::ResultSet)
    return rs.env_layer_md.timeframe
end

"""
    component_params(spec::DataFrame, component::Type)::DataFrame

Extract parameters for a specific model component from exported model specification.
"""
function component_params(rs::ResultSet, component::Type)::DataFrame
    spec = rs.model_spec
    return spec[spec.component.==string(component), :]
end


function Base.show(io::IO, mime::MIME"text/plain", rs::ResultSet)

    vers_id = rs.ADRIA_VERSION

    tf, sites, scens = size(rs.outcomes[:relative_cover])
    # Species/size groups represented: $(species)

    rcps = join(split(rs.RCP, "_"), ", ")

    println("""
    Domain: $(rs.name)

    Run with ADRIA $(vers_id) on $(rs.invoke_time)
    Results stored at: $(store_location(rs))

    RCP(s) represented: $(rcps)
    Intervention scenarios run: $(scens)
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
