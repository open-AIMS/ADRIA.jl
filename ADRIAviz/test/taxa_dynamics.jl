using WGLMakie, GeoMakie, GraphMakie
using ADRIA: DataCube, AnnotatedOutcomes
using ADRIAviz
using OrderedCollections

Makie.inline!(false)

function _taxonomy_ao_3d(; n_timesteps=10, n_groups=5, n_scenarios=12)
    data = DataCube(
        rand(n_timesteps, n_groups, n_scenarios);
        timesteps=1:n_timesteps, groups=1:n_groups, scenarios=1:n_scenarios
    )
    metadata = Dict{Symbol,Any}(
        :scenario_type_groups => OrderedDict{Symbol,BitVector}(
            :counterfactual => vcat(trues(4), falses(8)),
            :unguided => vcat(falses(4), trues(4), falses(4)),
            :guided => vcat(falses(8), trues(4))
        ),
        :scenario_rcp_groups => nothing
    )
    return AnnotatedOutcomes(data, metadata)
end

@testset "taxonomy" begin
    ao = _taxonomy_ao_3d()
    @testset "taxonomy(ao) returns Figure" begin
        @test ADRIA.viz.taxonomy(ao) isa Figure
    end
    @testset "taxonomy!(g, ao) returns GridLayout or GridPosition" begin
        f = Figure()
        g = f[1, 1] = GridLayout()
        @test ADRIA.viz.taxonomy!(g, ao) isa Union{GridLayout,GridPosition}
    end
end
