using Test

using ADRIA

using ADRIA: CSV, DataFrames, YAXArrays, YAXArrays.At

ENV["ADRIA_TEST"] = "true"

include("mock_data/mock_data.jl")

include("aqua.jl")

include("run_scenarios.jl")

# analysis.jl uses ADRIA.sensitivity.* and viz functions that have moved to ADRIAAnalysis.
# Run it from packages/ADRIAAnalysis/test/ instead.
include("clustering.jl")
include("annotated_outcomes.jl")
include("scenario_groups.jl")
include("data_loading.jl")
include("domain.jl")
include("calib_params_loading.jl")
include("growth.jl")
include("io/inputs.jl")
include("io/result_io.jl")
include("sampling.jl")
include("interventions/interventions.jl")
include("interventions/seeding.jl")
include("interventions/moving_corals.jl")
include("spec.jl")

include("example_run.jl")

include("decisions/mcda.jl")
include("decisions/location_spread.jl")

include("metrics/test_metrics_helper.jl")
include("metrics/scenario.jl")
include("metrics/metrics.jl")
include("metrics/spatial.jl")

include("utils/scale.jl")
include("utils/text_display.jl")

if get(ENV, "ADRIA_RUN_VIZ_TESTS", "0") == "1"
    include("viz/taxa_dynamics.jl")
    include("viz/spatial.jl")
    include("viz/annotated_outcomes.jl")
end

# TODO Fix spatial_clustering and site_selection tests
# include("site_selection.jl")
# include("spatial_clustering.jl")
#
# Always run this example test case last
# as it sets global environment variables

ENV["ADRIA_TEST"] = "false"
