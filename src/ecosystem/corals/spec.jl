using DataFrames
import ModelParameters: Model


"""
    to_spec(m::Model)::DataFrame

Convert Coral Model specification to a DataFrame of coral parameters
"""
function to_spec(m::Coral)::DataFrame
    _, pnames, spec = coral_spec()

    val_df = DataFrame(Model(m))
    res = copy(spec)
    fnames = String.(val_df[!, :fieldname])
    for p in pnames
        target = [occursin(p, fn) for fn in fnames]
        target_names = map(String, fnames[target])
        for tn in target_names
            c_id = rsplit(tn, "_$p", keepempty=false)
            res[res.coral_id .== c_id, p] = val_df[val_df.fieldname .== Symbol(tn), :val]
        end
    end

    return res
end