using DataFrames
import ModelParameters: Model


"""
    to_spec(m::Model)::DataFrame

Convert Coral Model specification to a coral spec DataFrame
"""
function to_spec(m::Coral)::DataFrame
    _, pnames, spec = coral_spec()
    val_df = DataFrame(Model(m))

    return _update_coral_spec(spec, pnames, val_df)
end


"""
    to_spec(coral_df::DataFrame)::DataFrame

Convert dataframe of model parameters to a coral spec.
"""
function to_spec(coral_df::DataFrame)::DataFrame
    _, pnames, spec = coral_spec()

    return _update_coral_spec(spec, pnames, coral_df)
end

function _update_coral_spec(spec::DataFrame, pnames::Vector{String}, coral_params::DataFrame)::DataFrame
    fnames = String.(coral_params[!, :fieldname])
    Threads.@threads for p in pnames
        target = occursin.(p, fnames)
        for (tn, sym_tn) in zip(fnames[target], Symbol.(fnames[target]))
            idx = spec.coral_id .== rsplit(tn, "_$p", keepempty=false)
            spec[idx, p] = coral_params[coral_params.fieldname.==sym_tn, :val]
        end
    end

    return spec
end
