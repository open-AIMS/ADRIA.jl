"""
    sample_stratified(d::Domain, n_env::Int64, n_lever::Int64;
        regimes::Tuple{Vararg{Symbol}}=(:guided, :unguided, :counterfactual),
        env_sample_method=SobolSample(R=OwenScramble(base=2, pad=32)),
        lever_sample_method=SobolSample(R=OwenScramble(base=2, pad=32)))::DataFrame

Draw a nested Environment -> Intervention policy -> Deployment strategy sample.

# Design

Level 1 (blocking, identical across every regime): `EnvironmentalLayer` factors
(`dhw_scenario`, `wave_scenario`, `cyclone_mortality_scenario`) and ecological/biological
factors (`Coral`, `GrowthAcceleration`). `n_env` blocks are drawn once via a space-filling
design and reused, unchanged, across every regime and every lever draw within that block.
Biological parameters are blocked here rather than nested under deployment strategy
because they interact multiplicatively with environment (thermal tolerance changes what
"cooler refugia" means) -- they are a source of exogenous uncertainty, not a lever.

Level 2/3 (independently sampled within each regime, per block): each regime's own active
factors -- MCDA criteria weights/plan_horizon/mcda_method for `:guided`; seeding/fogging/
moving-coral deployment amounts, depth thresholds and strategy timing for all regimes --
are drawn independently via `sample_guided`/`sample_unguided`/`sample_cf`, so RSA/PRIM
retains full independent, space-filling coverage of each regime's own lever space.

This differs from `sample_matched`, which derives unguided/CF rows from a single guided
draw by zeroing columns: that collapses unguided/CF lever coverage down to whatever the
guided draw happened to produce, and is suited to causal differencing (validating a single
candidate region), not open-ended scenario discovery. `sample_stratified` is suited to the
latter: it controls for the environment/ecology confound via blocking while leaving each
regime's own factor space free to vary for RSA/PRIM to search.

Rows carry `:env_block` (integer id of which Level-1 block they belong to, `1:n_env`) so
downstream analysis can facet by block for a 3-way, environment-controlled comparison, or
run RSA per-regime using each block as a repeated observation. Regime is identified by
`:guided` itself (`1` guided, `0` unguided, `-1` counterfactual) -- no separate `:run_type`
column, since it would just duplicate `:guided`.

# Arguments
- `d` : Domain.
- `n_env` : number of Level-1 environment/ecology blocks to draw (power of 2 for Sobol').
- `n_lever` : number of Level-2/3 lever draws per (block, regime) cell (power of 2 for
  Sobol').
- `regimes` : which policy regimes to include; subset of `(:guided, :unguided,
  :counterfactual)`.
- `env_sample_method` : sampler used for the Level-1 block draw.
- `lever_sample_method` : sampler used for each regime's Level-2/3 lever draw.

# Returns
Scenario specification with `n_env * length(regimes) * n_lever` rows.
"""
function sample_stratified(
    d::Domain, n_env::Int64, n_lever::Int64;
    regimes::Tuple{Vararg{Symbol}}=(:guided, :unguided, :counterfactual),
    env_sample_method=SobolSample(; R=OwenScramble(; base=2, pad=32)),
    lever_sample_method=SobolSample(; R=OwenScramble(; base=2, pad=32))
)::DataFrame
    regime_sampler = Dict(
        :guided => sample_guided,
        :unguided => sample_unguided,
        :counterfactual => sample_cf
    )
    unknown = filter(r -> !haskey(regime_sampler, r), regimes)
    isempty(unknown) || throw(
        ArgumentError(
            "Unknown regime(s) $unknown; expected a subset of " *
            "$(Tuple(keys(regime_sampler)))"
        )
    )

    # Level 1: environment + ecological/biological blocking factors, filtered
    # consistently with how sample_guided/sample_unguided/sample_cf filter their spec.
    cb_calib_groups::Vector{Int64} = d.loc_data.CB_CALIB_GROUPS
    full_spec = _filtered_model_spec(model_spec(d), cb_calib_groups)
    env_fieldnames = component_params(
        d.model, [EnvironmentalLayer, Coral, GrowthAcceleration]
    ).fieldname
    env_spec = full_spec[in.(full_spec.fieldname, Ref(env_fieldnames)), :]

    env_blocks = sample(env_spec, n_env, env_sample_method)
    env_cols = names(env_blocks)

    cells = DataFrame[]
    for b = 1:n_env
        for r in regimes
            df = regime_sampler[r](d, n_lever; sample_method=lever_sample_method)
            for col in env_cols
                df[!, col] .= env_blocks[b, col]
            end
            df[!, :env_block] .= b
            push!(cells, df)
        end
    end

    return vcat(cells...)
end
