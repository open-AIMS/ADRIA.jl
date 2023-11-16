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
    ENV["ADRIA_OUTPUT_DIR"] = mktempdir()

    # Run scenarios with example Domain
    dom = ADRIA.load_domain(EXAMPLE_DOMAIN_PATH)
    dom_path = joinpath(EXAMPLES_PATH, "example_scenarios.csv")
    scens = ADRIA.load_scenarios(dom, dom_path)

    return ADRIA.run_scenarios(dom, scens, "45")
end

const TEST_RS = test_rs()

include("clustering.jl")
include("data_loading.jl")
include("domain.jl")
include("growth.jl")
include("io/inputs.jl")
include("metrics.jl")
include("sampling.jl")
include("seeding.jl")
include("site_selection.jl")
include("spec.jl")
include("utils/text_display.jl")

# Always run this example test case last
# as it sets global environment variables
include("example_run.jl")
