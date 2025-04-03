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

"""
    sample_options(d::Domain, pd_frequency::Int64, sample_method=SobolSample(; R=OwenScramble(; base=2, pad=32)))::DataFrame

Generate sample where all default parameters are fixed and only option_ts varies.
"""
function sample_options(
    d::Domain, pd_frequency::Int64,
    sample_method=SobolSample(; R=OwenScramble(; base=2, pad=32))
)::DataFrame
    # Get one guided sample
    scens = sample(d, 2; sample_method=sample_method)[1:1, :]

    # Compute all possible option time series
    options = analysis.option_seed_preference()
    number_changes::Int64 = scens.seed_years[1] ÷ pd_frequency
    max_time::Int64 = size(d.dhw_scens[:, :, 1], 1)
    combinations = options_combinations(options.option_name, number_changes)
    options_ts = options_series(combinations, scens[1, :], pd_frequency, max_time)

    # Add decision frequency to scenario
    scens.pd_frequency = [pd_frequency]

    # Copy initial sample and only change option_ts
    scens = vcat([scens for _ in 1:length(combinations)]...)
    scens.option_ts = options_ts
    return scens
end

