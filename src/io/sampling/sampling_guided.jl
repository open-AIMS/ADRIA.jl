"""
    sample_guided(d::Domain, n::Int64; sample_method=SobolSample(R=OwenScramble(base=2, pad=32)))::DataFrame

Generate only guided scenarios.

# Arguments
- `d` : Domain.
- `n` : number of samples to generate.
- `sample_method` : type of sampler to use.

# Returns
Scenario specification
"""
function sample_guided(
    d::Domain, n::Int64; sample_method=SobolSample(; R=OwenScramble(; base=2, pad=32))
)::DataFrame
    cb_calib_groups::Vector{Int64} = d.loc_data.CB_CALIB_GROUPS
    spec_df = _filtered_model_spec(model_spec(d), cb_calib_groups)

    # Remove unguided scenarios as an option
    # Sample without unguided (i.e., values >= 1), then revert back to original model spec
    if !(spec_df[spec_df.fieldname .== :guided, :is_constant][1])
        _adjust_guided_lower_bound!(spec_df, 1)
        spec_df[!, :is_constant] .= spec_df[!, :lower_bound] .== spec_df[!, :upper_bound]
    end

    return sample(spec_df, n, sample_method)
end
