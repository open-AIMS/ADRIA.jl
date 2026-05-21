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
            contains(string(v.name), "Makie") for v in values(Base.loaded_modules)
        )
        @test !makie_loaded
    end

    if get(ENV, "ADRIA_RUN_VIZ_TESTS", "0") == "1"
        include("annotated_outcomes.jl")
    end
end
