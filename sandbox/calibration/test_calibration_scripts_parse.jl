using Test

@testset "Calibration scripts parse" begin
    for filename in [
        "calibration_run.jl",
        "calibration_parameters.jl",
        "cots_cycle_metrics.jl",
        "calibrate_cots_blackbox.jl",
        "simulate_best_calibration.jl",
        "validate_temporal_holdout.jl",
        "analyze_connectivity_factorial.jl",
        "diagnose_recurrence.jl",
        "run_recurrence_mechanism_factorial.jl",
        "run_separated_production_factorial.jl",
        "compare_lizard_v2_connectivity.jl",
        "compare_lizard_domain_gate.jl",
        "replay_lizard_peak_mechanisms.jl",
        "screen_lizard_stage_survival.jl",
        "reef_observation_mapping.jl",
        "diagnose_observation_smoothing.jl",
        "test_reef_observation_mapping.jl",
        "test_cots_forcing.jl",
        "audit_owen_flux_contract.jl",
        "screen_owen_embedded_mortality.jl",
        "replicate_owen_embedded_mortality.jl",
        "test_owen_boundary_shutdown.jl",
        "screen_owen_crash_parameters.jl",
        "screen_owen_recruitment_blackout.jl",
        "screen_owen_food_survival.jl",
        "diagnose_owen_food_trigger.jl",
        "screen_owen_low_cover_mortality.jl",
        "screen_owen_fecundity_consumption.jl",
        "screen_owen_source_condition_proxy.jl",
        "test_1991_initial_state.jl",
        "run_1991_initialization_test.jl",
        "score_1985_reference_on_1991_window.jl",
    ]
        source = read(joinpath(@__DIR__, filename), String)
        @test Meta.parseall(source; filename=filename) isa Expr
    end
end

@testset "COTS connectivity scripts parse" begin
    for filename in [
        joinpath("..", "domain_building", "build_lizard_cots_connectivity.jl"),
        joinpath("..", "domain_building", "validate_lizard_cots_connectivity.jl"),
        joinpath("..", "domain_building", "validate_lizard_domain_v2.jl"),
        joinpath("..", "domain_building", "smoke_lizard_domain_v2.jl"),
    ]
        source = read(joinpath(@__DIR__, filename), String)
        @test Meta.parseall(source; filename=filename) isa Expr
    end
end
