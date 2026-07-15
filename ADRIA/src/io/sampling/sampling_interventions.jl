"""
    sample_set(d::Domain, n::Int64, sample_method=SobolSample(R=OwenScramble(base=2, pad=32)))::DataFrame

Generate scenarios with pre-determined defaults.

# Arguments
- `d` : Domain.
- `n` : number of samples to generate.

# Returns
Scenario specification
"""
function sample_set(d::Domain, n::Int64, rcp::String)::DataFrame
    ADRIA.switch_RCPs!(d, rcp)

    # Get DHW trajectories
    dhws = d.dhw_scens

    # Get mean/std of trajectories
    traj_mean = dropdims((mean(dhws; dims=2)); dims=2)
    # traj_stdev = dropdims((std(dhws; dims=2)); dims=2)

    clusters = ADRIA.analysis.cluster_scenarios(traj_mean, min(50, size(dhws, 3)))

    Random.seed!(ceil(Int64, mean(traj_mean)))

    dhw_scens = [rand(findall(clusters .== c)) for c in unique(clusters)]

    # `guided`'s default bounds already span -1..1 (counterfactual/unguided/guided),
    # so only `mcda_method` needs constraining to a single method here.
    ADRIA.set_factor_bounds!(
        d;
        mcda_method=("COCOSO",),
        dhw_scenario=Tuple(dhw_scens)
    )

    # Deactivate superfluous environmental inputs
    ADRIA.fix_factor!(
        d;
        wave_scenario=0.0,
        cyclone_mortality_scenario=0.0
        # seed_strategy=0.0
    )

    # Assume coral model has been perfectly parameterized
    coral_params = ADRIA.component_params(d.model, ADRIA.Coral).fieldname
    ADRIA.fix_factor!(d, coral_params)

    # Assume the same for growth acceleration parameters
    growth_acc_params = ADRIA.component_params(d.model, ADRIA.GrowthAcceleration).fieldname
    ADRIA.fix_factor!(d, growth_acc_params)

    # Fix coral seeding weights
    seed_criteria_params = ADRIA.component_params(
        d.model,
        ADRIA.SeedCriteriaWeights
    )
    ADRIA.fix_factor!(d, seed_criteria_params.fieldname)

    # Fix moving coral weights
    # Fix coral seeding weights
    mc_criteria_params = ADRIA.component_params(
        d.model,
        ADRIA.MCCriteriaWeights
    )
    ADRIA.fix_factor!(d, mc_criteria_params.fieldname)

    scenarios = ADRIA.sample(d, n)

    return scenarios
end

"""
    _derive_regime(samples::DataFrame, guided_value::Float64)::DataFrame

Produce scenario rows paired to `samples` under a fixed `guided` regime, by setting
`:guided` to `guided_value` and fixing every column that `_PARAM_DEPENDENCIES` renders
inactive under that value.

Non-intervention parameter columns (and any columns still active under the given regime)
are left unchanged, enabling causal attribution by differencing outcomes between paired
rows.

Shared by `derive_cf` (`guided_value=-1.0`) and `derive_unguided` (`guided_value=0.0`) so
both regimes stay consistent with the same dependency table used pre-sampling by
`sample_guided`/`sample_unguided`/`sample_cf`.

# Arguments
- `samples` : intervention scenario rows (output of `sample_guided` or `sample_paired`).
- `guided_value` : the `:guided` sentinel for the target regime (-1.0 CF, 0.0 unguided).

# Returns
Scenario rows for the target regime, paired 1:1 with `samples`.
"""
function _derive_regime(samples::DataFrame, guided_value::Float64)::DataFrame
    regime = copy(samples)
    regime[!, :guided] .= guided_value

    # Mirrors the `haskey(known, rule.child)` guard in _resolve_conditional_spec!:
    # the two :strategy_group rows both evaluate "inactive" under guided=-1.0, and
    # without this guard the second (fix_to=PERIODIC) would overwrite the first
    # (fix_to=-1.0 CF sentinel) applied to the same columns.
    fixed_children = Set{Symbol}()
    for dep in _PARAM_DEPENDENCIES
        dep.child in fixed_children && continue

        active = _CONDITIONS[dep.op](guided_value, dep.value)
        dep.negate && (active = !active)
        active && continue  # only fix rows where this dep is inactive under this regime

        cols =
            dep.child in propertynames(regime) ?
            [dep.child] :
            get(_GROUP_MEMBERS, dep.child, Symbol[])
        for col in cols
            col in propertynames(regime) && (regime[!, col] .= dep.fix_to)
        end
        push!(fixed_children, dep.child)
    end

    return regime
end

"""
    derive_cf(samples::DataFrame, spec::DataFrame)::DataFrame

Produce counterfactual-regime rows paired to `samples` by zeroing all columns
that are inactive under guided = -1.0, using the pre-sampling dependency tables.

Non-intervention parameter columns are left unchanged, enabling causal attribution
by differencing outcomes between paired rows.

# Arguments
- `samples` : intervention scenario rows (output of `sample_guided` or `sample_paired`).
- `spec` : model spec for the same domain (output of `model_spec(d)`).
  Currently unused — the dependency tables are module-level constants — but retained
  for forward-compatibility with domain-specific dependency overrides.

# Returns
Counterfactual-regime scenario rows paired 1:1 with `samples`.
"""
function derive_cf(samples::DataFrame, spec::DataFrame)::DataFrame
    return _derive_regime(samples, -1.0)
end

"""
    derive_unguided(samples::DataFrame, spec::DataFrame)::DataFrame

Produce unguided-regime rows paired to `samples` by fixing all columns that are inactive
under guided = 0.0, using the pre-sampling dependency tables.

Unlike the counterfactual regime, intervention deployment amounts (`:intervention_group`)
and `:depth_thresholds` remain **active** (unchanged) under guided = 0.0 — unguided
scenarios still deploy interventions, just without MCDA-driven site selection. Only
`:criteria_weights`, `:plan_horizon`, and `:projection_confidence` are zeroed, and
`:strategy_group` columns
(`seed_strategy`/`fog_strategy`/`mc_strategy`) are fixed to `DECISION_STRATEGY[:periodic]`
— mirroring what `sample_unguided`'s pre-sampling `_resolve_conditional_spec!(spec,
(guided=0.0,))` does.

# Arguments
- `samples` : intervention scenario rows (output of `sample_guided` or `sample_paired`).
- `spec` : model spec for the same domain (output of `model_spec(d)`).
  Currently unused — the dependency tables are module-level constants — but retained
  for forward-compatibility with domain-specific dependency overrides.

# Returns
Unguided-regime scenario rows paired 1:1 with `samples`.
"""
function derive_unguided(samples::DataFrame, spec::DataFrame)::DataFrame
    return _derive_regime(samples, 0.0)
end

"""
Maps a `regimes` symbol (as accepted by `sample_matched`) to the `:guided` sentinel
value used to derive that regime from a guided draw via `_derive_regime`.
"""
const _REGIME_GUIDED_VALUE = Dict{Symbol,Float64}(
    :counterfactual => -1.0,
    :unguided => 0.0
)

"""
    sample_matched(d::Domain, n::Int; ssp=nothing, sample_method=SobolSample(R=OwenScramble(base=2, pad=32)), regimes=(:counterfactual, :unguided))::DataFrame

Draw `n` guided intervention scenarios and produce paired rows for each regime in
`regimes`. Returns a DataFrame with `(1 + length(regimes)) * n` rows and a `:run_type`
column (`String`) with values in `("intervention", "counterfactual", "unguided")`
(whichever were requested).

All rows share identical non-intervention parameter values (and, for `:unguided`,
identical intervention deployment amounts too — see `derive_unguided`), enabling causal
attribution by differencing outcomes across regimes.

This differs from `sample_balanced`, which draws each regime independently: those rows
are NOT matched row-for-row and cannot be used for this kind of differencing.

`n` should be a power of 2 for Sobol' sequences.

# Arguments
- `d` : Domain.
- `n` : number of guided intervention scenarios to draw.
- `ssp` : RCP/SSP scenario to switch the domain to before sampling, if given.
- `sample_method` : type of sampler to use.
- `regimes` : which regimes to derive alongside the guided draw. Must be a non-empty
  subset of `(:counterfactual, :unguided)`.

# Returns
Scenario specification with paired intervention + regime rows.
"""
function sample_matched(
    d::Domain, n::Int;
    ssp::Union{String,Nothing}=nothing,
    sample_method=SobolSample(; R=OwenScramble(; base=2, pad=32)),
    regimes::Tuple{Vararg{Symbol}}=(:counterfactual, :unguided)
)::DataFrame
    isempty(regimes) && throw(ArgumentError("`regimes` must be non-empty"))
    unknown = filter(r -> !haskey(_REGIME_GUIDED_VALUE, r), regimes)
    isempty(unknown) || throw(
        ArgumentError(
            "Unknown regime(s) $unknown; expected a subset of " *
            "$(Tuple(keys(_REGIME_GUIDED_VALUE)))"
        )
    )

    isnothing(ssp) || switch_RCPs!(d, ssp)

    intervention_samples = sample_guided(d, n; sample_method)
    intervention_samples[!, :run_type] .= "intervention"

    regime_dfs = [
        let df = _derive_regime(intervention_samples, _REGIME_GUIDED_VALUE[r])
            df[!, :run_type] .= string(r)
            df
        end for r in regimes
    ]

    return vcat(intervention_samples, regime_dfs...)
end

"""
    sample_paired(d::Domain, n::Int; ssp=nothing, sample_method=SobolSample(R=OwenScramble(base=2, pad=32)))::DataFrame

Draw `n` guided intervention scenarios and produce paired counterfactual rows.
Returns a DataFrame with `2n` rows and a `:run_type` column (`String`).

The intervention and CF rows share identical non-intervention parameter values,
enabling causal attribution by differencing outcomes.

`n` should be a power of 2 for Sobol' sequences.

A special case of `sample_matched` with `regimes=(:counterfactual,)`; see that function
for a 3-way (guided/CF/unguided) matched comparison.

# Arguments
- `d` : Domain.
- `n` : number of guided intervention scenarios to draw.
- `ssp` : RCP/SSP scenario to switch the domain to before sampling, if given.
- `sample_method` : type of sampler to use.

# Returns
Scenario specification with `2n` paired intervention/counterfactual rows.
"""
function sample_paired(
    d::Domain, n::Int;
    ssp::Union{String,Nothing}=nothing,
    sample_method=SobolSample(; R=OwenScramble(; base=2, pad=32))
)::DataFrame
    return sample_matched(d, n; ssp, sample_method, regimes=(:counterfactual,))
end

"""
    sample_matched_stratified(d::Domain, n_per_stratum::Int, ssp_strata::AbstractVector{String}; regimes=(:counterfactual, :unguided))::DataFrame

For each SSP in `ssp_strata`, draw `n_per_stratum` matched rows (see `sample_matched`)
using an independent Owen scramble. Returns a DataFrame with `:ssp` and `:run_type`
columns.

`n_per_stratum` must be a power of 2.

# Arguments
- `d` : Domain.
- `n_per_stratum` : number of guided intervention scenarios to draw per SSP stratum.
- `ssp_strata` : RCP/SSP scenarios to stratify sampling over.
- `regimes` : which regimes to derive alongside the guided draw; see `sample_matched`.

# Returns
Scenario specification with `(1 + length(regimes)) * n_per_stratum * length(ssp_strata)`
matched rows.
"""
function sample_matched_stratified(
    d::Domain,
    n_per_stratum::Int,
    ssp_strata::AbstractVector{String};
    regimes::Tuple{Vararg{Symbol}}=(:counterfactual, :unguided)
)::DataFrame
    results = map(ssp_strata) do ssp
        m = SobolSample(; R=OwenScramble(; base=2, pad=32))  # fresh scramble per stratum
        df = sample_matched(d, n_per_stratum; ssp, sample_method=m, regimes)
        df[!, :ssp] .= ssp
        df
    end
    return vcat(results...)
end

"""
    sample_paired_stratified(d::Domain, n_per_stratum::Int, ssp_strata::AbstractVector{String})::DataFrame

For each SSP in `ssp_strata`, draw `n_per_stratum` paired rows (see `sample_paired`)
using an independent Owen scramble. Returns a DataFrame with `:ssp` and `:run_type`
columns.

`n_per_stratum` must be a power of 2.

A special case of `sample_matched_stratified` with `regimes=(:counterfactual,)`.

# Arguments
- `d` : Domain.
- `n_per_stratum` : number of guided intervention scenarios to draw per SSP stratum.
- `ssp_strata` : RCP/SSP scenarios to stratify sampling over.

# Returns
Scenario specification with `2 * n_per_stratum * length(ssp_strata)` paired rows.
"""
function sample_paired_stratified(
    d::Domain,
    n_per_stratum::Int,
    ssp_strata::AbstractVector{String}
)::DataFrame
    return sample_matched_stratified(
        d, n_per_stratum, ssp_strata; regimes=(:counterfactual,)
    )
end
