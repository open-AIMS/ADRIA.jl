using Test

using ADRIA

include("mock_data/mock_data.jl")

include("aqua.jl")

include("run_scenarios.jl")

include("analysis.jl")
include("clustering.jl")
include("data_loading.jl")
include("domain.jl")
include("Ecosystem.jl")
include("growth.jl")
include("io/inputs.jl")
include("sampling.jl")
include("seeding.jl")
include("spec.jl")

include("decisions/mcda.jl")
include("decisions/location_spread.jl")

include("metrics/test_metrics_helper.jl")
include("metrics/scenario.jl")
include("metrics/metrics.jl")
include("metrics/spatial.jl")

include("utils/scale.jl")
include("utils/text_display.jl")

include("viz/taxa_dynamics.jl")
include("viz/spatial.jl")

# TODO Fix spatial_clustering and site_selection tests
# include("site_selection.jl")
# include("spatial_clustering.jl")
#
# Always run this example test case last
# as it sets global environment variables
# include("example_run.jl")
