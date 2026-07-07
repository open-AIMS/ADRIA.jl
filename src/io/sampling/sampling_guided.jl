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
    sample_options(d::Domain, pd_frequency::Int64, sample_fraction=1.0; sample_method=...)::DataFrame

Generate sample where all default parameters are fixed and only option_ts varies.

`sample_fraction` (∈ (0, 1]) controls what fraction of all pathway combinations to keep.
At 1.0 (default) all combinations are used; at e.g. 0.3, a random 30% are kept.
"""
function sample_options(
    d::Domain,
    pd_frequency::Int64,
    sample_fraction::Float64=1.0;
    sample_method=SobolSample(; R=OwenScramble(; base=2, pad=32))
)::DataFrame
    # Get one guided sample
    scens = sample(d, 2; sample_method=sample_method)[1:1, :]

    # Compute all possible option time series
    options = analysis.option_seed_preference()
    number_changes::Int64 = scens.seed_years[1] ÷ pd_frequency
    max_time::Int64 = size(d.dhw_scens, :timesteps)
    combinations = options_combinations(options.option_name, number_changes)
    options_ts = build_option_ts(combinations, scens[1, :], pd_frequency, max_time, sample_fraction)

    # Add decision frequency to scenario
    scens.pd_frequency = [pd_frequency]

    # Copy initial sample and only change option_ts
    scens = vcat([scens for _ in 1:length(options_ts)]...)
    scens.option_ts = options_ts
    return scens
end

"""
    options_combinations(options_name::Vector, number_repetitions::Int64)::Vector{Tuple}

Generate all combinations with repetition and considering order of pathway diversity options.
The first decision block uses only the real `options_name`; every subsequent block also allows
the `:nothing` option (do nothing this block — no seeding).
"""
#
function options_combinations(options_name::Vector, number_repetitions::Int64)::Vector{Tuple}
    with_nothing = vcat(options_name, :nothing)
    iterables = [i == 1 ? options_name : with_nothing for i in 1:number_repetitions]
    mat = collect(ADRIA.Iterators.product(iterables...))
    return vec(mat)
end

"""
    build_option_ts(combinations, scen, pd_frequency, max_time[, sample_fraction])::Vector{Int}

Encode each option combination as a base-6 integer. Use `decode_option_ts` to recover the
full time-series at model-run time.

`sample_fraction` (∈ (0, 1], default 1.0) randomly keeps that fraction of all combinations
before encoding. At 1.0 all combinations are kept.
"""
function build_option_ts(
    combinations::Vector{Tuple}, scen::DataFrameRow, pd_frequency::Int64, max_time::Int64,
    sample_fraction::Float64=1.0
)::Vector{Int}
    @assert length(combinations[1]) == scen.seed_years / pd_frequency
    @assert 0.0 < sample_fraction <= 1.0 "sample_fraction must be in (0, 1]"

    if sample_fraction < 1.0
        n_keep = ceil(Int, sample_fraction * length(combinations))
        combinations = StatsBase.sample(combinations, n_keep; replace=false)
    end

    return [ADRIA.analysis.encode_option_ts(combination) for combination in combinations]
end
