using Test
using TOML, CSV, DataFrames, ADRIA
using ADRIA.metrics: total_absolute_cover

const ADRIA_DIR = pkgdir(ADRIA)
const TEST_DATA_DIR = joinpath(ADRIA_DIR, "test", "data")
const EXAMPLE_DOMAIN_PATH = joinpath(ADRIA_DIR, "examples", "Example_domain")

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
