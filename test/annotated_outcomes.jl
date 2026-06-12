using ADRIA: DataCube, AnnotatedOutcomes, attach_scenario_metadata
using OrderedCollections

@testset "AnnotatedOutcomes" begin
    @testset "Direct construction" begin
        dummy_data = DataCube(rand(3, 4); timesteps=1:3, scenarios=1:4)
        meta = Dict{Symbol,Any}(:foo => "bar", :count => 42)
        ao = AnnotatedOutcomes(dummy_data, meta)

        @test ao isa AnnotatedOutcomes
        @test ao.data === dummy_data
        @test ao.metadata === meta
        @test ao.metadata[:foo] == "bar"
        @test ao.metadata[:count] == 42
    end

    @testset "Arbitrary metadata field access" begin
        dummy_data = DataCube(rand(2, 2); timesteps=1:2, scenarios=1:2)
        ao = AnnotatedOutcomes(
            dummy_data,
            Dict{Symbol,Any}(:x => [1, 2, 3], :y => nothing)
        )

        @test ao.metadata[:x] == [1, 2, 3]
        @test isnothing(ao.metadata[:y])
    end

    @testset "attach_scenario_metadata key contract" begin
        # Verify metadata keys without requiring real domain data.
        # attach_scenario_metadata is tested end-to-end in integration tests.
        n_scens = 8
        outcomes = DataCube(rand(5, n_scens); timesteps=1:5, scenarios=1:n_scens)
        groups = OrderedDict{Symbol,BitVector}(
            :counterfactual => vcat(trues(4), falses(4)),
            :guided => vcat(falses(4), trues(4))
        )
        ao = AnnotatedOutcomes(
            outcomes,
            Dict{Symbol,Any}(
                :scenario_type_groups => groups,
                :scenario_rcp_groups => nothing
            )
        )

        @test haskey(ao.metadata, :scenario_type_groups)
        @test haskey(ao.metadata, :scenario_rcp_groups)
        @test ao.metadata[:scenario_type_groups] isa AbstractDict{Symbol,BitVector}
        @test isnothing(ao.metadata[:scenario_rcp_groups])
    end
end
