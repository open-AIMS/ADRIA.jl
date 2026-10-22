# build/sim/precompile_sim.jl
#
# PackageCompiler execution file for adria_sim.so.
# Every code path exercised here is baked into the sysimage.
# Keep this representative of the actual EC2/Batch hot path.

using ADRIA

# ── Domain + simulation (mirrors existing build/precompile_script.jl) ──────
const FIXTURE_DIR = joinpath(pkgdir(ADRIA), "test", "data", "Test_domain")

if isdir(FIXTURE_DIR)
    ENV["ADRIA_DEBUG"] = "false"

    dom = ADRIA.load_domain(FIXTURE_DIR, "45")
    p_df = ADRIA.sample(dom, 32)
    rs = ADRIA.run_scenarios(dom, p_df, "45")

    # Core metrics – bake the most-called specialisations
    ADRIA.metrics.scenario_total_cover(rs)
    ADRIA.metrics.scenario_relative_cover(rs)
    ADRIA.metrics.scenario_coral_evenness(rs)
    ADRIA.metrics.scenario_absolute_shelter_volume(rs)

    delete!(ENV, "ADRIA_DEBUG")
    @info "Simulation precompile fixture executed successfully"
else
    @warn "Test fixture not found – sysimage will have reduced specialisation coverage" FIXTURE_DIR
end
