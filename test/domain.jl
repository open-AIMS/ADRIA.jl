using Test
using ADRIA

@testset "Domain loading" begin
    @testset "Domain DataFrame" begin
        dom = ADRIA.load_domain(EXAMPLE_DOMAIN_PATH, 45)
        p_df = ADRIA.param_table(dom)
        @test p_df isa DataFrame
    end

    @testset "Config" begin
        ADRIA.setup()

        # Ensure environment variables are set
        @test haskey(ENV, "ADRIA_OUTPUT_DIR")
        @test haskey(ENV, "ADRIA_NUM_CORES")
        @test haskey(ENV, "ADRIA_THRESHOLD")
    end

    #@testset "Discrete parameters" begin
    #    site_path = joinpath(TEST_DATA_DIR, "test_site_data.gpkg")
    #    conn_path = joinpath(TEST_DATA_DIR, "test_conn_data.csv")
    #    scen_path = joinpath(TEST_DATA_DIR, "test_scenarios.csv")
    #    dom = ADRIA.load_domain(EXAMPLE_DOMAIN_PATH)
    #    test_scens = CSV.read(scen_path, DataFrame)
    #    ADRIA.update_params!(dom, test_scens[5, :])
    #
    #    @test all(ADRIA.param_table(dom).N_seed_TA .== 500000)
    #end
end
