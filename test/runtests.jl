using Test
using TOML, CSV, DataFrames, ADRIA


const TEST_DATA_DIR = joinpath(@__DIR__, "data")


@testset "ADRIA.jl" begin
    # Write your tests here.
end


@testset "Config" begin
    ADRIA.setup()

    # Ensure environment variables are set
    @test haskey(ENV, "ADRIA_OUTPUT_DIR")
    @test haskey(ENV, "ADRIA_NUM_CORES")
    @test haskey(ENV, "ADRIA_reps")
    @test haskey(ENV, "ADRIA_THRESHOLD")

    # Check that the correct number of processors have been spun up.
    @eval using Distributed
    @test nprocs() == parse(Int, ENV["ADRIA_NUM_CORES"])
end

@testset "Domain loading" begin
    dom = ADRIA.load_domain(joinpath(@__DIR__, "..", "examples", "Example_domain"), 45)
end


@testset "site selection" begin
    # TODO: Complete tests with @tests

    dom = ADRIA.load_domain(joinpath(@__DIR__, "..", "examples", "Example_domain"), 45)

    p_tbl = ADRIA.param_table(dom)
    p_tbl.depth_offset .= 7.0
    # ranks = ADRIA.site_selection(dom, p_tbl, 1, 10, 1)
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
