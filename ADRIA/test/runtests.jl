using Test

using ADRIA

using ADRIA: CSV, DataFrames, YAXArrays, YAXArrays.At

ENV["ADRIA_TEST"] = "true"

include("test_constants.jl")

# Test sharding, used by `run_tests.jl` to run groups of tiers in separate processes.
# ADRIA_TEST_SHARD: "all" (default) | "fast" (tiers 1-2) | "spec" (tier 3) | "sim" (tiers 4-5)
# Tiers 4 and 5 share globals (TEST_RS etc.), so they must stay in one shard.
const SHARD = get(ENV, "ADRIA_TEST_SHARD", "all")
SHARD in ("all", "fast", "spec", "sim") || error("Unknown ADRIA_TEST_SHARD=$(repr(SHARD))")
_in_shard(names...) = SHARD == "all" || SHARD in names

# Helper: include a test file and report how long it took, on stderr (unbuffered) so the
# progress is still visible in the shard log if the launcher kills a shard at its time limit.
function _timed_include(path::String)
    t0 = time()
    println(stderr, "[test-timing] start  ", path)
    flush(stderr)
    include(path)
    println(
        stderr, "[test-timing] done   ", path, " in ", round(time() - t0; digits=1), " s"
    )
    flush(stderr)
    flush(stdout)
    return nothing
end

# Helper: abort if a tier recorded any failures or errors.
# Uses Test.get_test_counts() available in Julia 1.9+.
function _abort_if_failed(ts, tier_name::String)
    tc = Test.get_test_counts(ts)
    if tc.fails + tc.errors > 0
        @error "$tier_name failed ($(tc.fails) failures, $(tc.errors) errors) — aborting to skip expensive tests"
        exit(1)
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# Tier 1: Pre-Smoke (~30–60s)
# Pure-function and synthetic-data tests. No file I/O, no domain loads,
# no simulations. If these fail, there is no point running anything else.
# ─────────────────────────────────────────────────────────────────────────────
if _in_shard("fast")
    tier1 = @testset "Tier 1: Pre-Smoke" begin
        _timed_include("mock_data/mock_data.jl")
        _timed_include("aqua.jl")
        _timed_include("clustering.jl")
        _timed_include("Ecosystem.jl")
        _timed_include("growth.jl")
        _timed_include("ecosystem/cots.jl")
        _timed_include("ecosystem/low_dhw_bleaching.jl")
        _timed_include("scenario_groups.jl")
        _timed_include("annotated_outcomes.jl")
        _timed_include("utils/scale.jl")
        _timed_include("utils/text_display.jl")
    end
    _abort_if_failed(tier1, "Tier 1 (Pre-Smoke)")
end

# ─────────────────────────────────────────────────────────────────────────────
# Tier 2: I/O Smokescreen (~2–3 min)
# File I/O against real test data files. No simulations. Tests the
# data-loading layer before committing to expensive simulation runs.
# ─────────────────────────────────────────────────────────────────────────────
if _in_shard("fast")
    tier2 = @testset "Tier 2: I/O Smokescreen" begin
        _timed_include("data_loading.jl")
        _timed_include("connectivity.jl")
        _timed_include("io/inputs.jl")
        _timed_include("io/cots_inputs.jl")
        _timed_include("calib_params_loading.jl")
        # Cold-JIT budget for the work done in ADRIA's @compile_workload (~1 min; spawns
        # subprocesses). Guards against 20-30 min fresh-install precompile regressions.
        _timed_include("precompile_budget/precompile_budget.jl")
    end
    _abort_if_failed(tier2, "Tier 2 (I/O Smokescreen)")
end

# ─────────────────────────────────────────────────────────────────────────────
# Tier 3: Spec & Domain Loads (~2–3 min)
# Domain loading without full simulation. spec.jl MUST precede decisions/mcda.jl
# so that the ADRIA_DOM_45 global is shared via the @isdefined guard in mcda.jl.
# ─────────────────────────────────────────────────────────────────────────────
if _in_shard("spec")
    tier3 = @testset "Tier 3: Spec & Domain Loads" begin
        @info "Starting Tier 3"
        _timed_include("spec.jl")
        _timed_include("sampling.jl")
        _timed_include("growth_limits.jl")
        _timed_include("decisions/mcda.jl")
        _timed_include("decisions/location_spread.jl")
        _timed_include("decisions/revisit_cadence.jl")
    end
    _abort_if_failed(tier3, "Tier 3 (Spec & Domain Loads)")
end

# ─────────────────────────────────────────────────────────────────────────────
# Tier 4: Simulation Gate (~3–8 min)
# Runs 32 full coral-reef scenarios. Defines the TEST_RS, TEST_DOM, TEST_SCENS,
# and TEST_N_SAMPLES globals consumed by all Tier 5 tests.
# ─────────────────────────────────────────────────────────────────────────────
if _in_shard("sim")
    # Tier 5 metrics tests need MockData, which Tier 1 loads in the "fast" shard.
    @isdefined(MockData) || _timed_include("mock_data/mock_data.jl")

    tier4 = @testset "Tier 4: Simulation Gate" begin
        _timed_include("run_scenarios.jl")
    end
    _abort_if_failed(tier4, "Tier 4 (Simulation Gate)")

    # ─────────────────────────────────────────────────────────────────────────────
    # Tier 5: Simulation-Dependent Tests (~5–10 min)
    # All depend on TEST_RS produced by Tier 4.
    # analysis.jl has moved to ADRIAanalysis — run via ADRIAanalysis/test/ instead.
    # viz tests have moved to ADRIAviz/test/ — run via ADRIAviz test suite instead.
    # ─────────────────────────────────────────────────────────────────────────────
    @testset "Tier 5: Simulation-Dependent" begin
        _timed_include("domain.jl")
        _timed_include("io/result_io.jl")
        _timed_include("io/result_set_metadata.jl")
        _timed_include("example_run.jl")
        _timed_include("interventions/interventions.jl")
        _timed_include("interventions/CAq.jl")
        _timed_include("interventions/LvM.jl")
        _timed_include("determinism.jl")
        _timed_include("ecosystem/cots_scenarios.jl")
        _timed_include("ecosystem/low_dhw_scenarios.jl")
        _timed_include("metrics/test_metrics_helper.jl")
        _timed_include("metrics/scenario.jl")
        _timed_include("metrics/metrics.jl")
        _timed_include("metrics/spatial.jl")
        _timed_include("metrics/reef_indices.jl")
    end
end

# TODO Fix spatial_clustering and site_selection tests
# _timed_include("site_selection.jl")
# _timed_include("spatial_clustering.jl")

ENV["ADRIA_TEST"] = "false"
