using Test
using TOML, CSV, DataFrames, ADRIA
using ADRIA.metrics: total_absolute_cover

const ADRIA_DIR = pkgdir(ADRIA)
const TEST_DATA_DIR = joinpath(ADRIA_DIR, "test", "data")

const EXAMPLES_PATH = joinpath(ADRIA_DIR, "examples")
const EXAMPLE_DOMAIN_PATH = joinpath(EXAMPLES_PATH, "Example_domain")

function test_rs()
    # Load and apply configuration options
    ADRIA.setup()

    # Set result location to temporary folder within the current path
    ENV["ADRIA_OUTPUT_DIR"] = mktempdir(@__DIR__)

    # Run scenarios with example Domain
    dom = ADRIA.load_domain(EXAMPLE_DOMAIN_PATH)
    dom_path = joinpath(EXAMPLES_PATH, "example_scenarios.csv")
    scens = ADRIA.load_scenarios(dom, dom_path)

    return ADRIA.run_scenarios(dom, scens, "45")
end

const TEST_RS = test_rs()

@testset "proportional adjustment" begin
    Y = rand(5, 36, 20)
    tmp = zeros(20)
    max_cover = rand(20)

    for i in axes(Y, 1)
        ADRIA.proportional_adjustment!(Y[i, :, :], max_cover)

        @test all(0.0 .<= Y[i, :, :] .<= 1.0)
    end
end


include("clustering.jl")
include("data_loading.jl")
include("domain.jl")
include("growth.jl")
include("metrics.jl")
include("sampling.jl")
include("seeding.jl")
include("site_selection.jl")
include("spec.jl")

# Always run this example test case last
# as it sets global environment variables
include("example_run.jl")
