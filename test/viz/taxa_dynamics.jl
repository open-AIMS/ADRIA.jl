using WGLMakie, GeoMakie, GraphMakie

using ADRIA

if !@isdefined(TEST_RS)
    const TEST_DOM, TEST_N_SAMPLES, TEST_SCENS, TEST_RS = test_rs()
end

@testset "taxonomy" begin
    @testset "Default opts" begin
        @test ADRIA.viz.taxonomy(TEST_RS) isa Figure
    end

    # Default opts values
    opts::Dict{Symbol,Any} = Dict(
        :by_functional_groups => true, :by_RCP => false, :show_confints => true
    )

    @testset "Split by taxonomy" begin
        opts[:by_functional_groups] = false
        @test ADRIA.viz.taxonomy(TEST_RS; opts=opts) isa Figure
    end

    @testset "Split by RCP" begin
        opts[:by_functional_groups] = true
        opts[:by_RCP] = true
        @test ADRIA.viz.taxonomy(TEST_RS; opts=opts) isa Figure
    end

    @testset "Don't show confidence intervals" begin
        opts[:by_functional_groups] = true
        opts[:by_RCP] = true
        opts[:show_confints] = false
        @test ADRIA.viz.taxonomy(TEST_RS; opts=opts) isa Figure
    end
end
