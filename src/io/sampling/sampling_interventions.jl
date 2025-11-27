"""
    _adjust_guided_lower_bound!(guided_spec::DataFrame, lower::Int64)::DataFrame

Adjust lower bound of guided parameter spec to alter sampling range.
"""
function _adjust_guided_lower_bound!(spec_df::DataFrame, lower::Int64)::DataFrame
    guided_col = spec_df.fieldname .== :guided
    g_upper = Float64(spec_df[guided_col, :upper_bound][1])

    # Update entries, standardizing values for bounds as floats
    spec_df[guided_col, [:val, :lower_bound, :dist_params]] .=
        [lower Float64(lower) (Float64(lower), g_upper)]

    return spec_df
end

"""
    deactivate_interventions(to_update::DataFrame)::Nothing
    deactivate_interventions!(dom::Domain)::Nothing

Deactivate all intervention factors (excluding `guided`) by settings these to 0.0.

# Arguments
- `to_update` : model specification to modify/update
- `dom` : Domain

# Returns
Nothing
"""
function deactivate_interventions(to_update::DataFrame)::Nothing
    intervs = component_params(to_update, Intervention)
    cols = Symbol[fn for fn in intervs.fieldname if fn != :guided]
    for c in cols
        _row = to_update.fieldname .== c
        _dparams =
            length(to_update[_row, :dist_params][1]) == 2 ? (0.0, 0.0) : (0.0, 0.0, 0.0)

        dval = _is_discrete_factor(to_update[_row, :ptype][1]) ? 0 : 0.0
        to_update[_row, [:val, :lower_bound, :upper_bound, :dist_params, :is_constant]] .=
            [dval 0.0 0.0 _dparams true]
    end

    return nothing
end
function deactivate_interventions!(dom::Domain)::Nothing
    iv_factors = component_params(model_spec(dom), Intervention).fieldname
    criteria_factors =
        component_params(
            model_spec(dom),
            [SeedCriteriaWeights, FogCriteriaWeights, decision.DepthThresholds]
        ).fieldname
    n_iv_factors = length(iv_factors)
    n_criteria_factors = length(criteria_factors)
    ADRIA.fix_factor!(dom; NamedTuple{Tuple(iv_factors)}(zeros(Float64, n_iv_factors))...)
    ADRIA.fix_factor!(
        dom; NamedTuple{Tuple(criteria_factors)}(zeros(Float64, n_criteria_factors))...
    )

    return nothing
end

"""
    sample_guided(d::Domain, n::Int64, sample_method=SobolSample(R=OwenScramble(base=2, pad=32)))::DataFrame

Generate only guided scenarios.

# Arguments
- `d` : Domain.
- `n` : number of samples to generate.
- `sample_method` : type of sampler to use.

# Returns
Scenario specification
"""
function sample_guided(
    d::Domain, n::Int64, sample_method=SobolSample(; R=OwenScramble(; base=2, pad=32))
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