# build/analysis/precompile_analysis.jl
#
# PackageCompiler execution file for adria_analysis.so.
# Exercises the hot paths for load-results → cluster/stats → plot workflows.
# MLJ and SIRUS are intentionally NOT imported here.

using ADRIA
using ADRIAviz
using PlotlyBase           # triggers ADRIAvizPlotlyExt
using PlotlyKaleido        # triggers ADRIAvizPlotlyKaleidoExt (ADRIA.viz.savefig)

using ADRIAanalysis
using Statistics

# ── Simulation fixture (needed to produce a ResultSet for analysis) ─────────
const FIXTURE_DIR = joinpath(pkgdir(ADRIA), "test", "data", "Test_domain")

if isdir(FIXTURE_DIR)
    ENV["ADRIA_DEBUG"] = "false"

    dom = ADRIA.load_domain(FIXTURE_DIR, "45")
    p_df = ADRIA.sample(dom, 32)
    rs = ADRIA.run_scenarios(dom, p_df, "45")

    # ── Metrics ─────────────────────────────────────────────────────────────
    ADRIA.metrics.scenario_relative_cover(rs)

    metric_fns = [
        ADRIA.metrics.scenario_total_cover, ADRIA.metrics.scenario_absolute_shelter_volume
    ]
    outcomes = ADRIA.metrics.scenario_outcomes(rs, metric_fns)
    tc = outcomes[:, :, 1]
    s_tac = dropdims(mean(tc; dims=:timesteps); dims=:timesteps)
    s_sv = dropdims(mean(outcomes[:, :, 2]; dims=:timesteps); dims=:timesteps)

    # ── Clustering (ADRIA.analysis wraps Clustering.jl + Distances.jl) ───────
    ADRIA.analysis.cluster_scenarios(outcomes, 3)

    # ── RSA (ADRIAanalysis wraps HypothesisTests.jl) ─────────────────────────
    fs = ADRIAanalysis.feature_set(rs)
    ADRIAanalysis.rsa(fs, s_tac)

    # ── Counterfactual delta + bootstrap CI (ADRIAanalysis wraps Bootstrap.jl) ──
    p_cf = ADRIA.sample_cf(dom, 8)
    rs_cf = ADRIA.run_scenarios(dom, p_cf, "45")
    ADRIAanalysis.counterfactual_delta(
        rs, rs_cf,
        r -> dropdims(
            mean(ADRIA.metrics.scenario_total_cover(r); dims=:timesteps);
            dims=:timesteps
        );
        bootstrap_n=200
    )

    # ── DEA (ADRIAanalysis wraps DataEnvelopmentAnalysis.jl) ─────────────────
    cost = rand(Float64, length(s_tac))
    ADRIAanalysis.data_envelopment_analysis(cost, hcat(s_tac, s_sv))

    # ── Plotting ─────────────────────────────────────────────────────────────
    PlotlyKaleido.start()
    p = ADRIA.viz.scenarios(rs.inputs, tc)
    ADRIA.viz.savefig(p, joinpath(tempdir(), "_adria_precompile_analysis.png"))

    delete!(ENV, "ADRIA_DEBUG")
    @info "Analysis precompile fixture executed successfully"
else
    @warn "Test fixture not found – sysimage will have reduced specialisation coverage" FIXTURE_DIR
end
