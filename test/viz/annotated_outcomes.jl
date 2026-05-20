using WGLMakie, GeoMakie, GraphMakie

using ADRIA
using ADRIA: DataCube, AnnotatedOutcomes
using OrderedCollections

Makie.inline!(false)

# Synthetic fixtures

function _scenario_ao(; n_timesteps=10, n_scenarios=12)
    data = DataCube(rand(n_timesteps, n_scenarios); timesteps=1:n_timesteps, scenarios=1:n_scenarios)
    groups = OrderedDict{Symbol,BitVector}(
        :counterfactual => vcat(trues(4), falses(8)),
        :unguided       => vcat(falses(4), trues(4), falses(4)),
        :guided         => vcat(falses(8), trues(4)),
    )
    metadata = Dict{Symbol,Any}(
        :scenario_type_groups => groups,
        :scenario_rcp_groups  => OrderedDict{Symbol,BitVector}(
            :rcp45 => vcat(trues(6), falses(6)),
            :rcp60 => vcat(falses(6), trues(6)),
        ),
    )
    return AnnotatedOutcomes(data, metadata)
end

function _taxonomy_ao(; n_timesteps=10, n_groups=6, n_scenarios=12)
    data = DataCube(
        rand(n_timesteps, n_groups, n_scenarios);
        timesteps=1:n_timesteps,
        groups=1:n_groups,
        scenarios=1:n_scenarios,
    )
    sc_groups = OrderedDict{Symbol,BitVector}(
        :counterfactual => vcat(trues(4), falses(8)),
        :unguided       => vcat(falses(4), trues(4), falses(4)),
        :guided         => vcat(falses(8), trues(4)),
    )
    metadata = Dict{Symbol,Any}(
        :scenario_type_groups => sc_groups,
        :scenario_rcp_groups  => nothing,
    )
    return AnnotatedOutcomes(data, metadata)
end

@testset "AvizExt AnnotatedOutcomes dispatch" begin

    @testset "_get_scenario_groups error cases" begin
        @testset "Missing :scenario_type_groups throws ArgumentError" begin
            ao_empty = AnnotatedOutcomes(
                DataCube(rand(5, 4); timesteps=1:5, scenarios=1:4),
                Dict{Symbol,Any}(),
            )
            @test_throws ArgumentError ADRIA.viz.scenarios(ao_empty)
            err = try
                ADRIA.viz.scenarios(ao_empty)
                nothing
            catch e
                e
            end
            @test err isa ArgumentError
            @test occursin("attach_scenario_metadata", err.msg)
        end

        @testset "by_RCP=true with :scenario_rcp_groups => nothing throws ArgumentError" begin
            ao_rme = _taxonomy_ao()  # has :scenario_rcp_groups => nothing
            err = try
                ADRIA.viz.scenarios(ao_rme; by_RCP=true)
                nothing
            catch e
                e
            end
            @test err isa ArgumentError
            @test occursin("RME", err.msg)
        end
    end

    @testset "scenarios smoke tests" begin
        ao = _scenario_ao()

        @testset "scenarios(ao) returns Figure" begin
            @test ADRIA.viz.scenarios(ao) isa Figure
        end

        @testset "scenarios!(g, ao) returns GridLayout or GridPosition" begin
            f = Figure()
            g = f[1, 1] = GridLayout()
            result = ADRIA.viz.scenarios!(g, ao)
            @test result isa Union{GridLayout,GridPosition}
        end

        @testset "scenarios figure has Axis as first content element" begin
            fig = ADRIA.viz.scenarios(ao; legend=false)
            @test fig.content[1] isa Axis
        end

        @testset "scenarios figure has scene plots" begin
            fig = ADRIA.viz.scenarios(ao; legend=false)
            @test !isempty(fig.scene.plots)
        end
    end

    @testset "taxonomy smoke tests" begin
        ao = _taxonomy_ao()

        @testset "taxonomy(ao) returns Figure" begin
            @test ADRIA.viz.taxonomy(ao) isa Figure
        end

        @testset "taxonomy!(g, ao) returns GridLayout or GridPosition" begin
            f = Figure()
            g = f[1, 1] = GridLayout()
            result = ADRIA.viz.taxonomy!(g, ao)
            @test result isa Union{GridLayout,GridPosition}
        end
    end

    # Visual regression tests require committed snapshot baselines generated on the CI
    # platform (Linux + CairoMakie). Skipped here — add test/viz/avizext_visual.jl with
    # ReferenceTests.jl once baselines are committed.
end
