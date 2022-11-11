using DataFrames
import ModelParameters: Model


"""
    to_spec(m::Model)::DataFrame

Convert Coral Model specification to a DataFrame of coral parameters
"""
function to_spec(m::Coral)::DataFrame
    _, pnames, spec = coral_spec()
    val_df = DataFrame(Model(m))

    return _update_coral_spec(spec, pnames, val_df)
end


"""
    to_spec(coral_df::DataFrame)::DataFrame

Convert Coral Model specification to a DataFrame of coral parameters
"""
function to_spec(coral_df::DataFrame)::DataFrame
    _, pnames, spec = coral_spec()

    return _update_coral_spec(spec, pnames, coral_df)
end

function _update_coral_spec(spec::DataFrame, pnames::Vector{String}, coral_params::DataFrame)::DataFrame
    fnames = String.(coral_params[!, :fieldname])
    for p in pnames
        target = [occursin(p, fn) for fn in fnames]
        target_names = map(String, fnames[target])
        for tn in target_names
            c_id = rsplit(tn, "_$p", keepempty=false)
            spec[spec.coral_id.==c_id, p] = coral_params[coral_params.fieldname.==Symbol(tn), :val]
        end
    end

    return spec
end
