using Test
using ADRIA

@testset "Domain loading" begin
    @testset "Domain DataFrame" begin
        dom = ADRIA.load_domain(TEST_DOMAIN_PATH, 45)
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

    @testset "Discrete parameters" begin
        site_path = joinpath(TEST_DATA_DIR, "test_site_data.gpkg")
        conn_path = joinpath(TEST_DATA_DIR, "test_conn_data.csv")
        dom = ADRIA.load_domain(TEST_DOMAIN_PATH)

        # Create scenario spec
        samples = ADRIA.sample(dom, 8)
        samples[!, :N_seed_TA] .= 500_000.0

        # Write out scenario spec
        tmp_dir = mktempdir()
        tmp_fn = joinpath(tmp_dir, "test_scenarios.csv")
        CSV.write(tmp_fn, samples)

        # Read in scenario spec...
        test_scens = CSV.read(tmp_fn, DataFrame)

        # Update model using values from file
        ADRIA.update_params!(dom, test_scens[5, :])

        # Ensure values match
        @test all(ADRIA.param_table(dom).N_seed_TA .== 500000.0)

        # Ensure known discrete values are integer
        @test all(isinteger.(ADRIA.param_table(dom).seed_years))

        # Ensure known categoricals are integer values
        @test all(isinteger.(ADRIA.param_table(dom).guided))
    end
end
