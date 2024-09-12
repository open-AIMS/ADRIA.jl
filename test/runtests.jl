using Test

using ADRIA

include("factories/factories.jl")

const ADRIA_DIR = pkgdir(ADRIA)
const TEST_DATA_DIR = joinpath(ADRIA_DIR, "test", "data")

const TEST_DOMAIN_PATH = joinpath(TEST_DATA_DIR, "Test_domain")
const TEST_REEFMOD_ENGINE_DOMAIN_PATH = joinpath(TEST_DATA_DIR, "Reefmod_test_domain")

include("aqua.jl")

"""Test smaller scenario run with example scenario specification"""
function test_small_spec_rs()
    # Load and apply configuration options
    ADRIA.setup()

    # Run scenarios with example Domain
    dom = ADRIA.load_domain(TEST_DOMAIN_PATH)

    # Create scenario spec
    samples = ADRIA.sample(dom, 16)
    samples[!, :N_seed_TA] .= 500_000.0

    # Write out scenario spec
    tmp_dir = mktempdir()
    tmp_fn = joinpath(tmp_dir, "test_scenarios.csv")
    CSV.write(tmp_fn, samples)

    # Test reading in scenarios from a file
    scens = ADRIA.load_scenarios(dom, tmp_fn)

    return ADRIA.run_scenarios(dom, scens, "45")
end

"""Test ReefMod Engine domain loading"""
function test_reefmod_engine_domain()
    # Load and apply configuration options
    ADRIA.setup()

    # Load test RMEDomain
    return ADRIA.load_domain(RMEDomain, TEST_REEFMOD_ENGINE_DOMAIN_PATH, "45")
end

"""Test result set. Return domain, n_samples, scenarios and result set"""
function test_rs()
    # Load and apply configuration options
    ADRIA.setup()

    # Load domain data
    dom = ADRIA.load_domain(TEST_DOMAIN_PATH)

    # Set min_iv_locations upper to 10 as this is the number of locations used in tests
    dom = ADRIA.set_factor_bounds(dom, :min_iv_locations, (5.0, 10.0))

    # Create some scenarios
    # The number of scenarios set here seem to be the rough minimum for SIRUS to produce
    # some results.
    n_samples = 64
    scens = ADRIA.sample(dom, n_samples)

    # Run the model for generated scenarios
    rcp_45 = "45"
    rs = ADRIA.run_scenarios(dom, scens, rcp_45)

    return (dom, n_samples, scens, rs)
end

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
