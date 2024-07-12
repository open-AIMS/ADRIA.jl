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

include("domain.jl")
include("spec.jl")
include("io/inputs.jl")
include("io/data_loading.jl")
include("scenario_generation/sampling.jl")
include("ecosystem/Ecosystem.jl")
include("ecosystem/growth.jl")
include("ecosystem/seeding.jl")
include("metrics/metrics.jl")
include("analysis/clustering.jl")
include("decision/mcda.jl")
include("decision/rankings.jl")
include("decision/spatial_clustering.jl")
include("utils/text_display.jl")

# Always run this example test case last
# as it sets global environment variables
include("example_run.jl")
