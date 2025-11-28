"""
    sample_cf(d::Domain, n::Int64, sample_method=SobolSample(R=OwenScramble(base=2, pad=32)))::DataFrame

Generate only counterfactual scenarios.

# Arguments
- `d` : Domain.
- `n` : number of samples to generate.
- `sample_method` : type of sampler to use.

# Returns
Scenario specification
"""
function sample_cf(
    d::Domain, n::Int64, sample_method=SobolSample(; R=OwenScramble(; base=2, pad=32))
)::DataFrame
    cb_calib_groups::Vector{Int64} = d.loc_data.CB_CALIB_GROUPS
    spec_df = _filtered_model_spec(model_spec(d), cb_calib_groups)

    # Unguided scenarios only
    guided_col = spec_df.fieldname .== :guided
    spec_df[guided_col, [:val, :lower_bound, :upper_bound, :dist_params, :is_constant]] .=
        [-1 -1 -1 (-1.0, -1.0) true]

    # Remove intervention scenarios as an option
    deactivate_interventions(spec_df)

    return sample(spec_df, n, sample_method)
end
