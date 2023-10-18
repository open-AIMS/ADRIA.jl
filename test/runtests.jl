using Test
using TOML, CSV, DataFrames, ADRIA
using ADRIA.metrics: total_absolute_cover


const ADRIA_DIR = pkgdir(ADRIA)
const TEST_DATA_DIR = joinpath(ADRIA_DIR, "test", "data")
const EXAMPLE_DOMAIN_PATH = joinpath(ADRIA_DIR, "examples", "Example_domain")


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


    @testset "Discrete parameters" begin
        site_path = joinpath(TEST_DATA_DIR, "test_site_data.gpkg")
        conn_path = joinpath(TEST_DATA_DIR, "test_conn_data.csv")
        dom = ADRIA.load_domain(EXAMPLE_DOMAIN_PATH)

        # Create scenario spec
        samples = ADRIA.sample(dom, 8)
        samples[!, :N_seed_TA] .= 500_000

        # Write out scenario spec
        tmp_dir = mktempdir()
        tmp_fn = joinpath(tmp_dir, "test_scenarios.csv")
        CSV.write(tmp_fn, samples)

        # Read in scenario spec...
        test_scens = CSV.read(tmp_fn, DataFrame)

        # Update model using values from file
        ADRIA.update_params!(dom, test_scens[5, :])

        # Ensure values match
        @test all(ADRIA.param_table(dom).N_seed_TA .== 500000)

        # Ensure known discrete values are integer
        @test all(isinteger.(ADRIA.param_table(dom).seed_years))

        # Ensure known categoricals are integer values
        @test all(isinteger.(ADRIA.param_table(dom).guided))
    end
end


@testset "proportional adjustment" begin
    Y = rand(5, 36, 20)
    tmp = zeros(20)
    max_cover = rand(20)

    for i in axes(Y, 1)
        ADRIA.proportional_adjustment!(Y[i, :, :], tmp, max_cover)

        @test all(0.0 .<= Y[i, :, :] .<= 1.0)
    end
end


include("site_selection.jl")
include("data_loading.jl")
include("seeding.jl")
include("metrics.jl")
include("growth.jl")
include("spec.jl")
include("sampling.jl")
include("clustering.jl")

# Always run this example test case last
# as it sets global environment variables
include("example_run.jl")
