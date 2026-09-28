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
    ]
        source = read(joinpath(@__DIR__, filename), String)
        @test Meta.parseall(source; filename=filename) isa Expr
    end
end

@testset "COTS connectivity scripts parse" begin
    for filename in [
        joinpath("..", "domain_building", "build_lizard_cots_connectivity.jl"),
        joinpath("..", "domain_building", "validate_lizard_cots_connectivity.jl"),
    ]
        source = read(joinpath(@__DIR__, filename), String)
        @test Meta.parseall(source; filename=filename) isa Expr
    end
end
