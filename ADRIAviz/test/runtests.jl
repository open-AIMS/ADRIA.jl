using Test
using ADRIAviz

@testset "ADRIAviz" begin
    @testset "Package loads without backend" begin
        @test isdefined(ADRIAviz, :COLORS)
        @test isdefined(ADRIAviz, :_get_scenario_groups)
        @test isdefined(ADRIAviz, :_scenario_types)
        @test isdefined(ADRIAviz, :outcome_probability)
        @test isdefined(ADRIAviz, :OPT_TYPE)
        @test isdefined(ADRIAviz, :_time_labels)
    end

    @testset "No Makie in loaded modules" begin
        makie_loaded = any(
            contains(k.name, "Makie") for k in keys(Base.loaded_modules)
        )
        @test !makie_loaded
    end

    if get(ENV, "ADRIA_RUN_VIZ_TESTS", "0") == "1"
        include("annotated_outcomes.jl")
        include("taxa_dynamics.jl")
        include("sensitivity.jl")
        include("clustering.jl")
        # spatial.jl deferred -- requires real domain geometry
    end

    # Plotly backend tests.
    # Gate: ADRIA_RUN_PLOTLY_TESTS=1
    # Constraint: must run in a separate process / before any Makie `using`
    #             statement (Plotly and Makie extensions are mutually exclusive).
    if get(ENV, "ADRIA_RUN_PLOTLY_TESTS", "0") == "1"
        include("plotly.jl")
    end
end
