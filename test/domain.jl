using ADRIA.DataFrames

if !@isdefined(TEST_RS)
    const TEST_DOM, TEST_N_SAMPLES, TEST_SCENS, TEST_RS = test_rs()
end

@testset "Domain loading" begin
    @testset "Domain DataFrame" begin
        p_df = ADRIA.param_table(TEST_DOM)
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
        # Create scenario spec
        samples = deepcopy(TEST_SCENS)
        samples[!, :N_seed_TA] .= 500_000.0

        # Write out scenario spec
        tmp_dir = mktempdir()
        tmp_fn = joinpath(tmp_dir, "test_scenarios.csv")
        CSV.write(tmp_fn, samples)

        # Read in scenario spec...
        test_scens = CSV.read(tmp_fn, DataFrame)

        # Update model using values from file
        dom = deepcopy(TEST_DOM)
        ADRIA.update_params!(dom, test_scens[5, :])

        # Ensure values match
        @test all(ADRIA.param_table(dom).N_seed_TA .== 500000.0)

        # Ensure known discrete values are integer
        @test all(isinteger.(ADRIA.param_table(dom).seed_years))

        # Ensure known categoricals are integer values
        @test all(isinteger.(ADRIA.param_table(dom).guided))
    end
end
