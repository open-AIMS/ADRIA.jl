"""
    sample_unguided(d::Domain, n::Int64; sample_method=SobolSample(R=OwenScramble(base=2, pad=32)))::DataFrame

Generate only unguided scenarios.

# Arguments
- `d` : Domain.
- `n` : number of samples to generate.
- `sample_method` : type of sampler to use.

# Returns
Scenario specification
"""
function sample_unguided(
    d::Domain, n::Int64; sample_method=SobolSample(; R=OwenScramble(; base=2, pad=32))
)::DataFrame
    cb_calib_groups::Vector{Int64} = d.loc_data.CB_CALIB_GROUPS
    spec_df = _filtered_model_spec(model_spec(d), cb_calib_groups)

    # Fix guided factor to 0 (i.e., unguided scenarios only)
    guided_col = spec_df.fieldname .== :guided
    cols = [:val, :lower_bound, :upper_bound, :dist_params, :is_constant]
    spec_df[guided_col, cols] .= [0 0 0 (0.0, 0.0) true]

    return sample(spec_df, n, sample_method)
end
