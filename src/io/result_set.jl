import Dates: now

using PkgVersion, Zarr, NamedDims

import Setfield: @set!
import DataFrames: DataFrame
import ADRIA: Domain, EnvLayer


const RESULTS = "results"
const LOG_GRP = "logs"
const INPUTS = "inputs"


struct ResultSet{S, T, F, A, B}
    name::S
    rcp::Int
    invoke_time::S
    ADRIA_VERSION::S
    site_ids::T
    site_area::F
    site_max_coral_cover::F
    env_layer_md::EnvLayer

    inputs::DataFrame
    sim_constants::Dict

    # raw::AbstractArray
    outcomes::Dict
    ranks::A
    seed_log::B  # Values stored in m^2
    fog_log::A   # Reduction in bleaching mortality (0.0 - 1.0)
    shade_log::A # Reduction in bleaching mortality (0.0 - 1.0)    
end


function ResultSet(input_set::Zarr.ZArray, env_layer_md::EnvLayer, inputs_used::DataFrame, outcomes::Dict, log_set::Zarr.ZGroup)::ResultSet
    ResultSet(input_set.attrs["name"],
              input_set.attrs["rcp"],
              input_set.attrs["invoke_time"],
              input_set.attrs["ADRIA_VERSION"],
              input_set.attrs["site_ids"],
              convert.(Float64, input_set.attrs["site_area"]),
              convert.(Float64, input_set.attrs["site_max_coral_cover"]),
              env_layer_md,
              inputs_used,
              input_set.attrs["sim_constants"],
              outcomes,
              NamedDimsArray{Symbol.(Tuple(log_set["rankings"].attrs["structure"]))}(log_set["rankings"]),
              NamedDimsArray{Symbol.(Tuple(log_set["seed"].attrs["structure"]))}(log_set["seed"]),
              NamedDimsArray{Symbol.(Tuple(log_set["fog"].attrs["structure"]))}(log_set["fog"]),
              NamedDimsArray{Symbol.(Tuple(log_set["shade"].attrs["structure"]))}(log_set["shade"]))
end



"""
    store_name(r::ResultSet)::String

Get name of result set.
"""
function store_name(r::ResultSet)::String
    return "$(r.name)__$(r.rcp)__$(r.invoke_time)"
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


function Base.show(io::IO, mime::MIME"text/plain", rs::ResultSet)

    vers_id = rs.ADRIA_VERSION

    tf, sites, reps, scens = size(rs.outcomes[:relative_cover])
    # Species/size groups represented: $(species)

    println("""
    Domain: $(rs.name)

    Run with ADRIA $(vers_id) on $(rs.invoke_time) for RCP $(rs.rcp)
    Results stored at: $(store_location(rs))

    Intervention scenarios run: $(scens)
    Environmental scenarios: $(reps)
    Number of sites: $(sites)
    Timesteps: $(tf)

    Input layers
    ------------""")

    for fn in fieldnames(typeof(rs.env_layer_md))
        println("$(fn) : $(getfield(rs.env_layer_md, fn))")
    end
end
