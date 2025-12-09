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

    clusters = ADRIA.analysis.cluster_scenarios(traj_mean, min(20, size(dhws, 2)))

    Random.seed!(ceil(Int64, mean(traj_mean)))

    dhw_scens = [rand(findall(clusters .== c)) for c in unique(clusters)]

    ADRIA.set_factor_bounds!(
        d;
        guided=("counterfactual", "unguided", "COCOSO"),
        dhw_scenario=Tuple(dhw_scens)
    )

    # Deactivate superfluous environmental inputs
    ADRIA.fix_factor!(
        d;
        wave_scenario=0.0,
        cyclone_mortality_scenario=0.0
    )

    # Assume coral model has been perfectly parameterized
    coral_params = ADRIA.component_params(d.model, ADRIA.Coral).fieldname
    ADRIA.fix_factor!(d, coral_params)

    # Fix coral seeding weights
    seed_criteria_params = ADRIA.component_params(
        d.model,
        ADRIA.SeedCriteriaWeights
    )

    ADRIA.fix_factor!(d, seed_criteria_params.fieldname)

    scenarios = ADRIA.sample(d, n)

    return scenarios
end
