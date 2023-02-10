using DataFrames
import ModelParameters: Model


"""
    to_coral_spec(m::Model)::DataFrame

Convert Coral Model specification to a coral spec DataFrame
"""
function to_coral_spec(m::Coral)::DataFrame
    _, pnames, spec = coral_spec()
    val_df = DataFrame(Model(m))

    return _update_coral_spec(spec, pnames, val_df)
end


"""
    to_coral_spec(coral_df::DataFrame)::DataFrame

Convert dataframe of model parameters to a coral spec.
"""
function to_coral_spec(coral_df::DataFrame)::DataFrame
    _, pnames, spec = coral_spec()

    return _update_coral_spec(spec, pnames, coral_df)
end

function _update_coral_spec(spec::DataFrame, pnames::Vector{String}, coral_params::DataFrame)::DataFrame
    fnames::Vector{String} = String.(coral_params[!, :fieldname])
    for p in pnames
        target::Bool = occursin.(p, fnames)
        Threads.@threads for tn in fnames[target]
            idx::Vector{Bool} = spec.coral_id .== rsplit(tn, "_$p", keepempty=false)
            spec[idx, [p]] .= coral_params[coral_params.fieldname.==tn, :val]
        end
    end

    return spec
end

function to_coral_spec(inputs::NamedDimsArray)::DataFrame
    _, pnames, spec = coral_spec()

    coral_ids::Vector{String} = spec[:, :coral_id]
    for p in pnames
        # wrapping `p` in an array is necessary so update of DF works
        spec[!, [p]] .= Array(inputs(coral_ids .* "_" .* p))
    end

    return spec
end
function to_coral_spec(inputs::DataFrameRow)::DataFrame
    ins = NamedDimsArray(Vector(inputs), factors=names(inputs))
    return to_coral_spec(ins)
end
