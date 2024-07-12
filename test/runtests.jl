using Test
using Statistics

using ADRIA
using ADRIA.TOML, ADRIA.CSV, ADRIA.DataFrames
using ADRIA.metrics: total_absolute_cover

using WGLMakie, GeoMakie, GraphMakie

include("aqua.jl")
include("test_helpers.jl")

if !@isdefined(TEST_RS)
    # Run full example with figure creation when running full test suite
    const TEST_RS = test_rs_w_fig()
end

include("clustering.jl")
include("data_loading.jl")
include("domain.jl")
include("Ecosystem.jl")
include("growth.jl")
include("io/inputs.jl")
include("metrics.jl")
include("sampling.jl")
include("seeding.jl")
include("spec.jl")
include("mcda.jl")
include("utils/text_display.jl")

# TODO Fix spatial_clustering and site_selection tests
# include("site_selection.jl")
# include("spatial_clustering.jl")

# Always run this example test case last
# as it sets global environment variables
include("example_run.jl")
