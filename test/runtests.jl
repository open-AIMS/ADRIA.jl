using Test

using ADRIA

include("factories/factories.jl")

include("aqua.jl")

include("result_set.jl")

include("metrics/test_metrics_helper.jl")
# Set result location to temporary folder within the current path
include("clustering.jl")
include("data_loading.jl")
include("domain.jl")
include("Ecosystem.jl")
include("growth.jl")
include("io/inputs.jl")
include("sampling.jl")
include("seeding.jl")
include("spec.jl")
include("mcda.jl")
include("metrics/metrics.jl")
include("utils/scale.jl")
include("utils/text_display.jl")
include("viz.jl")

# TODO Fix spatial_clustering and site_selection tests
# include("site_selection.jl")
# include("spatial_clustering.jl")
#
# Always run this example test case last
# as it sets global environment variables
# include("example_run.jl")
