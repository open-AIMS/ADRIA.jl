using Test
using TOML, CSV, DataFrames, ADRIA


const TEST_DATA_DIR = joinpath(@__DIR__, "data")


@testset "Domain loading" begin
    dom = ADRIA.load_domain(joinpath(@__DIR__, "..", "examples", "Example_domain"), 45)
end


@testset "Config" begin
    ADRIA.setup()

    # Ensure environment variables are set
    @test haskey(ENV, "ADRIA_OUTPUT_DIR")
    @test haskey(ENV, "ADRIA_NUM_CORES")
    @test haskey(ENV, "ADRIA_THRESHOLD")

    # Check that the correct number of processors have been spun up.
    @eval using Distributed
    @test nprocs() == parse(Int, ENV["ADRIA_NUM_CORES"])
end


@testset "Discrete parameters" begin
    site_path = joinpath(TEST_DATA_DIR, "test_site_data.gpkg")
    conn_path = joinpath(TEST_DATA_DIR, "test_conn_data.csv")
    scen_path = joinpath(TEST_DATA_DIR, "test_scenarios.csv")

    dom = ADRIA.load_domain(joinpath(@__DIR__, "..", "examples", "Example_domain"), 45)

    test_scens = CSV.read(scen_path, DataFrame)
    ADRIA.update_params!(dom, test_scens[5, :])

    @test all(ADRIA.param_table(dom).seed_TA .== 500000)
end


@testset "proportional adjustment" begin
    Y = rand(5, 36, 20)
    tmp = zeros(20)
    max_cover = rand(1, 20)

    for i in axes(Y,1)
        adjusted = ADRIA.proportional_adjustment!(Y[i, :, :], tmp, max_cover)

        @test all(adjusted .<= 1.0)
        @test all(adjusted .>= 0.0)
    end
end

include("site_selection.jl")

include("metrics.jl")
include("modifying_metrics.jl")