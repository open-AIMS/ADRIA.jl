# build/analysis/precompile_analysis.jl
#
# PackageCompiler execution file for adria_analysis.so.
# Exercises the hot paths for load-results → cluster/stats → plot workflows.
# MLJ and SIRUS are intentionally NOT imported here.

using ADRIA
using ADRIAviz
using CairoMakie           # triggers ADRIAvizMakieExt

using ADRIAanalysis
using Bootstrap
using Clustering
using HypothesisTests
using DataEnvelopmentAnalysis

# ── Simulation fixture (needed to produce a ResultSet for analysis) ─────────
const FIXTURE_DIR = joinpath(dirname(dirname(@__DIR__)), "test", "data", "Test_domain")

if isdir(FIXTURE_DIR)
    ENV["ADRIA_DEBUG"] = "false"

    dom = ADRIA.load_domain(FIXTURE_DIR, "45")
    p_df = ADRIA.sample(dom, 32)
    rs = ADRIA.run_scenarios(dom, p_df, "45")

    # ── Metrics ─────────────────────────────────────────────────────────────
    tc = ADRIA.metrics.scenario_total_cover(rs)
    ADRIA.metrics.scenario_relative_cover(rs)

    # ── Clustering (generic Clustering.jl call; ADRIAanalysis wraps these) ──
    data_mat = rand(Float64, 10, 32)
    km = kmeans(data_mat, 3; maxiter=20)
    _ = km.assignments

    hc = hclust(pairwise(Distances.Euclidean(), data_mat); linkage=:ward)
    cutree(hc; k=3)

    # ── Bootstrap ───────────────────────────────────────────────────────────
    bs = bootstrap(mean, rand(Float64, 100), BasicSampling(200))
    confint(bs, BasicConfInt(0.95))

    # ── HypothesisTests ─────────────────────────────────────────────────────
    x, y = rand(Float64, 50), rand(Float64, 50)
    pvalue(MannWhitneyUTest(x, y))
    pvalue(KruskalWallisTest(x, y))

    # ── DEA ─────────────────────────────────────────────────────────────────
    X = rand(Float64, 5, 10)
    Y = rand(Float64, 2, 10)
    dea(X, Y)

    # ── Plotting ─────────────────────────────────────────────────────────────
    fig = ADRIAviz.ADRIA_outcomes(rs)
    save(joinpath(tempdir(), "_adria_precompile_analysis.png"), fig)

    delete!(ENV, "ADRIA_DEBUG")
    @info "Analysis precompile fixture executed successfully"
else
    @warn "Test fixture not found – sysimage will have reduced specialisation coverage" FIXTURE_DIR
end
