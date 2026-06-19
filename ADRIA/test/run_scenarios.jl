using ADRIA

const ADRIA_DIR = pkgdir(ADRIA)
const TEST_DATA_DIR = joinpath(ADRIA_DIR, "test", "data")

const TEST_DOMAIN_PATH = joinpath(TEST_DATA_DIR, "Test_domain")

# Load and apply configuration options
ADRIA.setup()

# Set result location to temporary folder within the current path
ENV["ADRIA_OUTPUT_DIR"] = mktempdir()

"""Test smaller scenario run with example scenario specification"""
function test_small_spec_rs()
    # Run scenarios with example Domain
    dom = ADRIA.load_domain(TEST_DOMAIN_PATH, "45")

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
    # Load test RMEDomain
    return ADRIA.load_domain(RMEDomain, TEST_REEFMOD_ENGINE_DOMAIN_PATH, "45")
end

"""Test result set. Return domain, n_samples, scenarios and result set"""
function test_rs()
    # Load domain data
    dom = ADRIA.load_domain(TEST_DOMAIN_PATH, "45")

    # Set min_iv_locations upper to 10 as this is the number of locations used in tests
    dom = ADRIA.set_factor_bounds!(dom, :min_iv_locations, (5.0, 10.0))

    # Create some scenarios
    # Keep this large enough to exercise the shared ResultSet-based tests without paying
    # the full cost of the old 64-scenario run.
    n_samples = 32
    scens = ADRIA.sample(dom, n_samples)

    # Run the model for generated scenarios
    rcp_45 = "45"
    rs = ADRIA.run_scenarios(dom, scens, rcp_45)

    return (dom, n_samples, scens, rs)
end

# ── Globals consumed by all Tier 5 tests ─────────────────────────────────────
if !@isdefined(TEST_RS)
    const TEST_DOM, TEST_N_SAMPLES, TEST_SCENS, TEST_RS = test_rs()
end

# ── Tier 4 gate tests ─────────────────────────────────────────────────────────

@testset "Scenario run produces a valid ResultSet" begin
    @test TEST_RS isa ADRIA.ResultSet
    @test size(TEST_RS.inputs, 1) == TEST_N_SAMPLES
    @test TEST_N_SAMPLES == 32
    # Ensure all expected log arrays are present
    @test !isnothing(TEST_RS.seed_log)
    @test !isnothing(TEST_RS.mc_log)
    @test !isnothing(TEST_RS.shading_log)
    @test !isnothing(TEST_RS.coral_dhw_tol_log)
    @test !isnothing(TEST_RS.ranks)
end

@testset "CSV spec scenario loading" begin
    rs_small = test_small_spec_rs()
    @test rs_small isa ADRIA.ResultSet
    # 16 scenarios were sampled in test_small_spec_rs
    @test size(rs_small.inputs, 1) == 16
end
