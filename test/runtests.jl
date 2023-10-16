using Test
using TOML, CSV, DataFrames, ADRIA
using ADRIA.metrics: total_absolute_cover

const ADRIA_DIR = pkgdir(ADRIA)
const TEST_DATA_DIR = joinpath(ADRIA_DIR, "test", "data")
const EXAMPLE_DOMAIN_PATH = joinpath(ADRIA_DIR, "examples", "Example_domain")

function test_rs()
    orig_loc = pwd()
    this_loc = @__DIR__
    cd(this_loc)

    # Run full example to make sure nothing errors
    ADRIA.setup()  # Load and apply configuration options

    # Use a temporary directory for result location
    ENV["ADRIA_OUTPUT_DIR"] = mktempdir()

    ex_domain = ADRIA.load_domain("../examples/Example_domain")
    scens = ADRIA.load_scenarios(ex_domain, "../examples/example_scenarios.csv")
    rs = ADRIA.run_scenarios(ex_domain, scens, "45")

    cd(orig_loc)

    return rs
end

const TEST_RS = test_rs()

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
