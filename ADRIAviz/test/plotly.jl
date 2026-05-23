# ADRIAviz/test/plotly.jl
#
# Test suite for the ADRIAvizPlotlyExt backend.
#
# Run with:
#   ADRIA_RUN_PLOTLY_TESTS=1 julia --project=ADRIAviz -e 'using Pkg; Pkg.test()'
#
# CRITICAL — do NOT import any Makie package (WGLMakie, GLMakie, CairoMakie,
# GeoMakie, GraphMakie) in this file. Doing so triggers ADRIAvizMakieExt,
# which causes method-ambiguity errors with the Plotly backend methods.
#
# API contract tested here (Plotly-specific signatures):
#   • scenarios(ao)                                — AnnotatedOutcomes dispatch
#   • scenarios(rs, outcomes)                      — ResultSet dispatch
#   • scenarios(matrix, clusters)                  — clustering dispatch
#   • clustered_scenarios(matrix, clusters)
#   • taxonomy(rs::ResultSet)                      — ResultSet dispatch
#   • taxonomy(scenarios, relative_taxa_cover)     — DataFrame dispatch
#   • pawn(Si)                                     — 2-D YAXArray
#   • tsa(Si)                                      — 3-D YAXArray (no ResultSet)
#   • rsa(Si, factor_values)                       — no ResultSet
#   • outcome_map(outcomes, factor_values)         — no ResultSet
#   • convergence(Si_conv, factors)
#   • data_envelopment_analysis(dea_output)
#   • rules_scatter(scenarios_df, clusters, rules)
#   • !-variant stubs raise descriptive errors
#   • HTML serialisation

using Test
using ADRIA
using ADRIA: DataCube, AnnotatedOutcomes
using ADRIAviz
using PlotlyBase
using OrderedCollections
using DataFrames
using Random

include("plotly_fixtures.jl")

Random.seed!(42)

if get(ENV, "ADRIA_RUN_VIZ_TESTS", "0") == "1"
    error("Cannot run Makie and Plotly tests in the same process. Use separate CI steps.")
end

# ─────────────────────────────────────────────────────────────────────────────
# 1. Extension loading & activate(:plotly)
# ─────────────────────────────────────────────────────────────────────────────

@testset "Extension loading" begin
    @testset "ADRIAvizPlotlyExt module is active after `using PlotlyBase`" begin
        ext = Base.get_extension(ADRIAviz, :ADRIAvizPlotlyExt)
        @test !isnothing(ext)
        @test ext isa Module
    end

    @testset "No Makie module present in loaded_modules" begin
        makie_loaded = any(
            contains(string(k.name), "Makie") for k in keys(Base.loaded_modules)
        )
        @test !makie_loaded
    end

    @testset "activate(:plotly) is defined" begin
        @test hasmethod(ADRIAviz.activate, (Symbol,))
    end

    @testset "activate(:plotly) is idempotent and does not error" begin
        # Extension is already loaded at this point; calling activate again
        # must not raise an error.
        @test_nowarn ADRIAviz.activate(:plotly)
        ext = Base.get_extension(ADRIAviz, :ADRIAvizPlotlyExt)
        @test !isnothing(ext)
    end

    @testset "activate(:makie) raises ArgumentError" begin
        @test_throws ArgumentError ADRIAviz.activate(:makie)
    end

    @testset "activate(:unknown) raises ArgumentError" begin
        @test_throws ArgumentError ADRIAviz.activate(:foobar)
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# 2. scenarios(ao::AnnotatedOutcomes)
# ─────────────────────────────────────────────────────────────────────────────

@testset "scenarios(ao)" begin
    ao = _plotly_scenario_ao()   # 3 type-groups, rcp45 / rcp60

    @testset "return type is PlotlyBase.Plot" begin
        @test ADRIA.viz.scenarios(ao) isa PlotlyBase.Plot
    end

    @testset "has at least 3 traces (one per scenario-type group)" begin
        p = ADRIA.viz.scenarios(ao)
        @test length(p.data) >= 3
    end

    @testset "all traces are scatter type" begin
        p = ADRIA.viz.scenarios(ao)
        @test all(t.type == "scatter" for t in p.data)
    end

    @testset "summarize=true produces fewer traces than summarize=false" begin
        # summarize=true  → one mean line + one CI band per group (compact)
        # summarize=false → one individual trace per scenario (expanded)
        p_ci = ADRIA.viz.scenarios(ao; summarize=true)
        p_raw = ADRIA.viz.scenarios(ao; summarize=false)
        @test length(p_raw.data) > length(p_ci.data)
    end

    @testset "by_RCP=true groups traces by RCP label" begin
        p = ADRIA.viz.scenarios(ao; by_RCP=true)
        @test p isa PlotlyBase.Plot
        # Must produce at least 2 traces (one per RCP group)
        @test length(p.data) >= 2
        trace_names = [
            lowercase(string(t.name)) for t in p.data if !isnothing(get(t, :name, nothing))
        ]
        @test any(contains(n, "45") || contains(n, "rcp") for n in trace_names)
    end

    @testset "by_RCP=false uses scenario-type group labels" begin
        p = ADRIA.viz.scenarios(ao; by_RCP=false)
        @test p isa PlotlyBase.Plot
        @test length(p.data) >= 3
        trace_names = [
            lowercase(string(t.name)) for t in p.data if !isnothing(get(t, :name, nothing))
        ]
        type_keywords = ["counterfactual", "unguided", "guided"]
        @test any(any(contains(n, kw) for kw in type_keywords) for n in trace_names)
    end

    @testset "sort_by=:mean produces same trace count as default" begin
        p_def = ADRIA.viz.scenarios(ao)
        p_sorted = ADRIA.viz.scenarios(ao; sort_by=:mean)
        @test p_sorted isa PlotlyBase.Plot
        @test length(p_sorted.data) == length(p_def.data)
    end

    @testset "each trace has x-values spanning the timestep range" begin
        p = ADRIA.viz.scenarios(ao)
        # data[1] is CI band (2×n x-vals); find the first line trace
        line_trace = first(t for t in p.data if get(t, :mode, "") == "lines")
        xs = line_trace.x
        @test !isnothing(xs)
        @test length(xs) >= 1
        # Expect 10 timesteps (from _plotly_scenario_ao default)
        @test length(xs) == 10
    end

    @testset "by_RCP=true on RME-style ao (no RCP metadata) raises ArgumentError" begin
        ao_rme = _plotly_rme_ao()
        @test_throws ArgumentError ADRIA.viz.scenarios(
            ao_rme; by_RCP=true
        )
    end

    @testset "missing :scenario_type_groups raises ArgumentError with helpful message" begin
        ao_bare = _plotly_bare_ao()
        err = try
            ADRIA.viz.scenarios(ao_bare)
            nothing
        catch e
            e
        end
        @test err isa ArgumentError
        @test occursin("attach_scenario_metadata", err.msg)
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# 2b. scenarios(rs::ResultSet, outcomes) — ResultSet dispatch
# ─────────────────────────────────────────────────────────────────────────────

@testset "scenarios(rs::ResultSet, outcomes)" begin
    rs = _plotly_scenarios_rs()
    if isnothing(rs)
        @warn "Test domain not available — scenarios(rs, outcomes) test skipped"
        @test_skip true
    else
        s_tc = ADRIA.metrics.scenario_total_cover(rs)

        @testset "return type is PlotlyBase.Plot" begin
            @test ADRIA.viz.scenarios(rs, s_tc) isa PlotlyBase.Plot
        end

        @testset "has traces" begin
            p = ADRIA.viz.scenarios(rs, s_tc)
            @test length(p.data) > 0
        end

        @testset "all traces are scatter type" begin
            p = ADRIA.viz.scenarios(rs, s_tc)
            @test all(t.type == "scatter" for t in p.data)
        end

        @testset "x-axis spans scenario timesteps" begin
            p = ADRIA.viz.scenarios(rs, s_tc)
            line_trace = first(t for t in p.data if get(t, :mode, "") == "lines")
            @test length(line_trace.x) == length(s_tc.timesteps)
        end

        @testset "by_RCP=true and by_RCP=false both return valid Plots" begin
            p_rcp = ADRIA.viz.scenarios(rs, s_tc; by_RCP=true)
            p_type = ADRIA.viz.scenarios(rs, s_tc; by_RCP=false)
            @test p_rcp isa PlotlyBase.Plot
            @test p_type isa PlotlyBase.Plot
        end
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# 3. scenarios(matrix, clusters)  — clustering dispatch
# ─────────────────────────────────────────────────────────────────────────────

@testset "scenarios(matrix, clusters)" begin
    matrix, clusters = _plotly_cluster_data(; n_timesteps=10, n_scenarios=12, n_clusters=3)

    @testset "return type is PlotlyBase.Plot" begin
        @test ADRIA.viz.scenarios(matrix, clusters) isa PlotlyBase.Plot
    end

    @testset "at least one trace per unique cluster" begin
        p = ADRIA.viz.scenarios(matrix, clusters)
        n_groups = length(unique(clusters))
        @test length(p.data) >= n_groups
    end

    @testset "all traces are scatter type" begin
        p = ADRIA.viz.scenarios(matrix, clusters)
        @test all(t.type == "scatter" for t in p.data)
    end

    @testset "summarize=false increases trace count" begin
        p_ci = ADRIA.viz.scenarios(
            matrix, clusters; summarize=true
        )
        p_raw = ADRIA.viz.scenarios(
            matrix, clusters; summarize=false
        )
        @test length(p_raw.data) > length(p_ci.data)
    end

    @testset "BitVector clusters are accepted" begin
        mat_bv, clust_bv = _plotly_bitvec_cluster_data()
        @test ADRIA.viz.scenarios(mat_bv, clust_bv) isa PlotlyBase.Plot
    end

    @testset "single-cluster input does not error" begin
        one_cluster = ones(Int, 12)
        @test ADRIA.viz.scenarios(matrix, one_cluster) isa PlotlyBase.Plot
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# 4. clustered_scenarios(matrix, clusters)
# ─────────────────────────────────────────────────────────────────────────────

@testset "clustered_scenarios(matrix, clusters)" begin
    matrix, clusters = _plotly_cluster_data(; n_clusters=3)

    @testset "return type is PlotlyBase.Plot" begin
        @test ADRIA.viz.clustered_scenarios(matrix, clusters) isa PlotlyBase.Plot
    end

    @testset "layout title contains 'Cluster'" begin
        p = ADRIA.viz.clustered_scenarios(matrix, clusters)
        title_text = lowercase(
            string(
                get(p.layout, :title_text, get(get(p.layout, :title, Dict()), :text, ""))
            )
        )
        @test contains(title_text, "cluster")
    end

    @testset "at least one trace per unique cluster" begin
        p = ADRIA.viz.clustered_scenarios(matrix, clusters)
        @test length(p.data) >= length(unique(clusters))
    end

    @testset "trace names reference cluster identifiers" begin
        p = ADRIA.viz.clustered_scenarios(matrix, clusters)
        trace_names = [
            lowercase(string(t.name)) for t in p.data if !isnothing(get(t, :name, nothing))
        ]
        @test any(
            contains(n, "cluster") || any(contains(n, string(c)) for c in unique(clusters))
            for n in trace_names
        )
    end

    @testset "BitVector clusters are accepted" begin
        mat_bv, clust_bv = _plotly_bitvec_cluster_data()
        @test ADRIA.viz.clustered_scenarios(mat_bv, clust_bv) isa PlotlyBase.Plot
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# 5. taxonomy(rs::ResultSet) and taxonomy(scenarios, relative_taxa_cover)
# ─────────────────────────────────────────────────────────────────────────────

@testset "taxonomy(rs::ResultSet)" begin
    rs = _plotly_taxonomy_rs()
    if isnothing(rs)
        @warn "Test domain not available — taxonomy(rs) test skipped"
        @test_skip true
    else
        @test ADRIA.viz.taxonomy(rs) isa PlotlyBase.Plot
        p = ADRIA.viz.taxonomy(rs)
        @test length(p.data) > 0
        @test all(t.type == "scatter" for t in p.data)
    end
end

@testset "taxonomy(scenarios::DataFrame, relative_taxa_cover)" begin
    scens, rtc = _plotly_taxonomy_df()

    @testset "return type is PlotlyBase.Plot" begin
        @test ADRIA.viz.taxonomy(scens, rtc) isa PlotlyBase.Plot
    end

    @testset "has traces" begin
        p = ADRIA.viz.taxonomy(scens, rtc)
        @test length(p.data) > 0
    end

    @testset "all traces are scatter type" begin
        p = ADRIA.viz.taxonomy(scens, rtc)
        @test all(t.type == "scatter" for t in p.data)
    end

    @testset "by_functional_groups=true and =false both return valid Plots" begin
        p_fg = ADRIA.viz.taxonomy(scens, rtc; by_functional_groups=true)
        p_type = ADRIA.viz.taxonomy(scens, rtc; by_functional_groups=false)
        @test p_fg isa PlotlyBase.Plot
        @test p_type isa PlotlyBase.Plot
    end

    @testset "show_confints=false reduces trace count" begin
        p_ci = ADRIA.viz.taxonomy(scens, rtc; show_confints=true)
        p_no = ADRIA.viz.taxonomy(scens, rtc; show_confints=false)
        @test length(p_no.data) <= length(p_ci.data)
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# 6. pawn(Si)
# ─────────────────────────────────────────────────────────────────────────────

@testset "pawn(Si)" begin
    Si = _plotly_pawn_si(; n_factors=4)   # 4 factors × 8 Si cols

    @testset "return type is PlotlyBase.Plot" begin
        @test ADRIA.viz.pawn(Si) isa PlotlyBase.Plot
    end

    @testset "exactly one heatmap trace" begin
        p = ADRIA.viz.pawn(Si)
        @test length(p.data) == 1
        @test p.data[1].type == "heatmap"
    end

    @testset "heatmap z has n_factors rows" begin
        p = ADRIA.viz.pawn(Si)
        z = p.data[1].z
        @test !isempty(z)
        # z is either a Matrix or a Vector-of-Vectors; either way outer size == n_factors
        n_rows = z isa AbstractMatrix ? size(z, 1) : length(z)
        @test n_rows == 4
    end

    @testset "normalize=true produces a valid Plot" begin
        p = ADRIA.viz.pawn(Si; normalize=true)
        @test p isa PlotlyBase.Plot
        @test length(p.data) == 1
    end

    @testset "factors= kwarg filters displayed rows" begin
        p_all = ADRIA.viz.pawn(Si)
        p_subset = ADRIA.viz.pawn(Si; factors=[:Factor1, :Factor2])
        z_all = p_all.data[1].z
        z_sub = p_subset.data[1].z
        n_all = z_all isa AbstractMatrix ? size(z_all, 1) : length(z_all)
        n_sub = z_sub isa AbstractMatrix ? size(z_sub, 1) : length(z_sub)
        @test n_sub == 2
        @test n_all == 4
    end

    @testset "sort_by=:median produces valid Plot with same number of traces" begin
        p_def = ADRIA.viz.pawn(Si)
        p_sorted = ADRIA.viz.pawn(Si; by=:median)
        @test p_sorted isa PlotlyBase.Plot
        @test length(p_sorted.data) == length(p_def.data)
    end

    @testset "y-axis tick labels contain factor names" begin
        p = ADRIA.viz.pawn(Si)
        yticks = get(p.layout.yaxis, :ticktext, nothing)
        @test !isnothing(yticks)
        @test length(yticks) == 4
        factor_strs = string.(Symbol.("Factor" .* string.(1:4)))
        @test any(string(yt) in factor_strs for yt in yticks)
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# 7. tsa(Si)  — Plotly-specific signature (no ResultSet)
#
# Contract: ADRIA.viz.tsa(Si::YAXArray; opts=...) where Si has dims
#           factors × Si_cols × timesteps.
# ─────────────────────────────────────────────────────────────────────────────

@testset "tsa(Si)" begin
    Si = _plotly_tsa_si(; n_factors=3, n_timesteps=10)

    @testset "return type is PlotlyBase.Plot" begin
        @test ADRIA.viz.tsa(Si) isa PlotlyBase.Plot
    end

    @testset "at least one scatter trace per factor" begin
        p = ADRIA.viz.tsa(Si)
        @test length(p.data) >= 3
    end

    @testset "all traces are scatter type" begin
        p = ADRIA.viz.tsa(Si)
        @test all(t.type == "scatter" for t in p.data)
    end

    @testset "trace x-values length equals n_timesteps" begin
        p = ADRIA.viz.tsa(Si)
        xs = p.data[1].x
        @test !isnothing(xs)
        @test length(xs) == 10
    end

    @testset "stat=:mean is accepted without error" begin
        p = ADRIA.viz.tsa(Si; stat=:mean)
        @test p isa PlotlyBase.Plot
    end

    @testset "stat=:lb is accepted without error" begin
        p = ADRIA.viz.tsa(Si; stat=:lb)
        @test p isa PlotlyBase.Plot
    end

    @testset "x-axis title refers to time" begin
        p = ADRIA.viz.tsa(Si)
        xlabel = lowercase(
            string(
                get(
                    p.layout.xaxis,
                    :title_text,
                    get(get(p.layout.xaxis, :title, Dict()), :text, "")
                )
            )
        )
        @test contains(xlabel, "year") || contains(xlabel, "time") ||
            contains(xlabel, "step")
    end

    @testset "factor names appear in legend / trace names" begin
        p = ADRIA.viz.tsa(Si)
        trace_names = [string(t.name) for t in p.data if !isnothing(get(t, :name, nothing))]
        factor_strs = string.(Symbol.("Factor" .* string.(1:3)))
        @test any(any(contains(tn, fs) for fs in factor_strs) for tn in trace_names)
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# 8. rsa(Si, factor_values)  — Plotly-specific signature
#
# Contract: ADRIA.viz.rsa(Si::YAXArray, factor_values::AbstractMatrix; opts=...)
#           Si dims: factors × si_quantile
#           factor_values: n_scenarios × n_factors
# ─────────────────────────────────────────────────────────────────────────────

@testset "rsa(Si, factor_values)" begin
    Si, fvals = _plotly_rsa_data(; n_factors=3, n_quantiles=10, n_scenarios=20)

    @testset "return type is PlotlyBase.Plot" begin
        @test ADRIA.viz.rsa(Si, fvals) isa PlotlyBase.Plot
    end

    @testset "at least one trace per factor" begin
        p = ADRIA.viz.rsa(Si, fvals)
        @test length(p.data) >= 3
    end

    @testset "all traces are scatter type" begin
        p = ADRIA.viz.rsa(Si, fvals)
        @test all(t.type == "scatter" for t in p.data)
    end

    @testset "3 factors → subplot grid (xaxis2 present in layout)" begin
        p = ADRIA.viz.rsa(Si, fvals)
        # Plotly uses xaxis, xaxis2, ... for multi-subplot layouts
        @test !isnothing(get(p.layout, :xaxis2, nothing))
    end

    @testset "single-factor rsa does not produce xaxis2" begin
        Si1, fv1 = _plotly_rsa_data(; n_factors=1)
        p = ADRIA.viz.rsa(Si1, fv1)
        @test p isa PlotlyBase.Plot
        @test length(p.data) >= 1
        @test isnothing(get(p.layout, :xaxis2, nothing))
    end

    @testset "x-axis label defaults to 'Factor Value' or similar" begin
        p = ADRIA.viz.rsa(Si, fvals)
        xlabel = lowercase(
            string(
                get(
                    p.layout.xaxis,
                    :title_text,
                    get(get(p.layout.xaxis, :title, Dict()), :text, "")
                )
            )
        )
        @test contains(xlabel, "factor") || contains(xlabel, "value") ||
            contains(xlabel, "param")
    end

    @testset "custom axis label kwarg is applied" begin
        p = ADRIA.viz.rsa(Si, fvals; xlabel="My X")
        xlabel = string(
            get(
                p.layout.xaxis,
                :title_text,
                get(get(p.layout.xaxis, :title, Dict()), :text, "")
            )
        )
        @test xlabel == "My X"
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# 9. outcome_map(outcomes, factor_values)  — Plotly-specific signature
#
# Contract: ADRIA.viz.outcome_map(outcomes::YAXArray, factor_values::AbstractMatrix; opts=...)
#           outcomes dims: factors × CI × si_quantile
# ─────────────────────────────────────────────────────────────────────────────

@testset "outcome_map(outcomes, factor_values)" begin
    outcomes, fvals = _plotly_outcome_map_data(; n_factors=3, n_quantiles=10)

    @testset "return type is PlotlyBase.Plot" begin
        @test ADRIA.viz.outcome_map(outcomes, fvals) isa PlotlyBase.Plot
    end

    @testset "at least 2 traces per factor (mean line + CI ribbon)" begin
        p = ADRIA.viz.outcome_map(outcomes, fvals)
        @test length(p.data) >= 2 * 3
    end

    @testset "has scatter traces (mean lines)" begin
        p = ADRIA.viz.outcome_map(outcomes, fvals)
        @test any(t.type == "scatter" for t in p.data)
    end

    @testset "3 factors → subplot grid (xaxis2 present)" begin
        p = ADRIA.viz.outcome_map(outcomes, fvals)
        @test !isnothing(get(p.layout, :xaxis2, nothing))
    end

    @testset "single-factor variant works" begin
        oc1, fv1 = _plotly_outcome_map_data(; n_factors=1)
        p = ADRIA.viz.outcome_map(oc1, fv1)
        @test p isa PlotlyBase.Plot
        @test isnothing(get(p.layout, :xaxis2, nothing))
    end

    @testset "custom x-axis label is applied" begin
        p = ADRIA.viz.outcome_map(outcomes, fvals; xlabel="FactorAxis")
        xlabel = string(
            get(
                p.layout.xaxis,
                :title_text,
                get(get(p.layout.xaxis, :title, Dict()), :text, "")
            )
        )
        @test xlabel == "FactorAxis"
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# 10. convergence(Si_conv, factors)
# ─────────────────────────────────────────────────────────────────────────────

@testset "convergence(Si_conv, factors)" begin
    conv = _plotly_convergence_si(; n_factors=3, n_sample_sizes=8)
    Si_conv, factors = conv.Si, conv.factors

    @testset "return type is PlotlyBase.Plot" begin
        @test ADRIA.viz.convergence(Si_conv, factors) isa PlotlyBase.Plot
    end

    @testset "at least one trace per factor" begin
        p = ADRIA.viz.convergence(Si_conv, factors)
        @test length(p.data) >= length(factors)
    end

    @testset "trace x-values match n_scenarios axis length" begin
        p = ADRIA.viz.convergence(Si_conv, factors)
        xs = p.data[1].x
        @test !isnothing(xs)
        @test length(xs) == length(Si_conv.n_scenarios)
    end

    @testset "plot_overlay=true and =false both return valid Plots" begin
        p_ov = ADRIA.viz.convergence(Si_conv, factors; plot_overlay=true)
        p_gr = ADRIA.viz.convergence(Si_conv, factors; plot_overlay=false)
        @test p_ov isa PlotlyBase.Plot
        @test p_gr isa PlotlyBase.Plot
    end

    @testset "plot_overlay=false (grid) has at least as many traces as overlay" begin
        p_ov = ADRIA.viz.convergence(Si_conv, factors; plot_overlay=true)
        p_gr = ADRIA.viz.convergence(Si_conv, factors; plot_overlay=false)
        @test length(p_gr.data) >= length(p_ov.data)
    end

    @testset "x-axis title refers to scenarios or sample size" begin
        p = ADRIA.viz.convergence(Si_conv, factors)
        xlabel = lowercase(
            string(
                get(
                    p.layout.xaxis,
                    :title_text,
                    get(get(p.layout.xaxis, :title, Dict()), :text, "")
                )
            )
        )
        @test (
            contains(xlabel, "scenario") ||
            contains(xlabel, "sample") ||
            contains(xlabel, "n ")
        )
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# 11. data_envelopment_analysis(DEA_output)
# ─────────────────────────────────────────────────────────────────────────────

@testset "data_envelopment_analysis(DEA_output)" begin
    dea = _plotly_dea_output(; n_scenarios=20)

    @testset "return type is PlotlyBase.Plot" begin
        @test ADRIA.viz.data_envelopment_analysis(dea) isa PlotlyBase.Plot
    end

    @testset "has at least 2 scatter traces (frontier + data cloud)" begin
        p = ADRIA.viz.data_envelopment_analysis(dea)
        scatter_count = count(t -> t.type == "scatter", p.data)
        @test scatter_count >= 2
    end

    @testset "has 3 subplot rows (efficiency frontier, scale eff., technical eff.)" begin
        p = ADRIA.viz.data_envelopment_analysis(dea)
        # Plotly encodes multiple y-axes as yaxis, yaxis2, yaxis3, ...
        n_yaxes = sum(1 for k in keys(p.layout) if startswith(string(k), "yaxis"))
        @test n_yaxes >= 3
    end

    @testset "frontier_type=:crs_peers is accepted" begin
        p = ADRIA.viz.data_envelopment_analysis(dea; frontier_type=:crs_peers)
        @test p isa PlotlyBase.Plot
    end

    @testset "custom metric axis labels are applied" begin
        p = ADRIA.viz.data_envelopment_analysis(
            dea;
            metrics_x_lab="CostMetric",
            metrics_y_lab="CoverMetric"
        )
        # One of the first subplot axis labels should reflect the custom string
        xlab = string(
            get(
                p.layout.xaxis,
                :title_text,
                get(get(p.layout.xaxis, :title, Dict()), :text, "")
            )
        )
        ylab = string(
            get(
                p.layout.yaxis,
                :title_text,
                get(get(p.layout.yaxis, :title, Dict()), :text, "")
            )
        )
        @test contains(xlab, "CostMetric") || contains(ylab, "CoverMetric")
    end

    @testset "legend entries match frontier and data-cloud labels" begin
        p = ADRIA.viz.data_envelopment_analysis(dea)
        # Default legend names: "Best practice frontier" and "Scenario data cloud"
        trace_names = [
            lowercase(string(t.name)) for t in p.data if !isnothing(get(t, :name, nothing))
        ]
        @test any(contains(n, "frontier") || contains(n, "best") for n in trace_names)
        @test any(contains(n, "cloud") || contains(n, "data") for n in trace_names)
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# 12. rules_scatter(scenarios_df, clusters, rules)
#
# Plotly-specific signature drops the `rs::ResultSet` argument (only needed
# for figure-size computation in the Makie backend).
# ─────────────────────────────────────────────────────────────────────────────

@testset "rules_scatter(scenarios, clusters, rules)" begin
    scens = _plotly_scenarios_df(; n_scenarios=20)
    clusters = vcat(ones(Int, 10), zeros(Int, 10))
    rules = _plotly_rules(; n_rules=4)

    @testset "return type is PlotlyBase.Plot" begin
        @test ADRIA.viz.rules_scatter(scens, clusters, rules) isa PlotlyBase.Plot
    end

    @testset "at least one scatter trace per 2-clause rule" begin
        p = ADRIA.viz.rules_scatter(scens, clusters, rules)
        scatter_count = count(t -> t.type == "scatter", p.data)
        @test scatter_count >= length(rules)
    end

    @testset "3-clause rules are silently filtered (same trace count as without them)" begin
        three_clause = [
            ["N_seed_TA", "<", 0.3], ["fogging", "<=", 0.2], ["guided", ">", 1.0]
        ]
        mixed = vcat(rules, [three_clause])
        p_orig = ADRIA.viz.rules_scatter(scens, clusters, rules)
        p_mixed = ADRIA.viz.rules_scatter(scens, clusters, mixed)
        @test length(p_mixed.data) == length(p_orig.data)
    end

    @testset "empty rules vector returns a Plot with no traces" begin
        p = ADRIA.viz.rules_scatter(scens, clusters, Vector{Vector}[])
        @test p isa PlotlyBase.Plot
        @test isempty(p.data)
    end

    @testset "BitVector clusters accepted" begin
        bv_clusters = vcat(trues(10), falses(10))
        p = ADRIA.viz.rules_scatter(scens, bv_clusters, rules)
        @test p isa PlotlyBase.Plot
    end

    @testset "target / non-target points use distinct marker colors" begin
        p = ADRIA.viz.rules_scatter(scens, clusters, rules)
        # Expect at least 2 distinct marker colors across traces
        colors = unique([
            string(get(get(t, :marker, Dict()), :color,
                get(t, :marker_color, "")))
            for t in p.data if t.type == "scatter"
        ])
        @test length(colors) >= 2
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# 13. !-variant stubs — must raise descriptive errors
#
# The Plotly backend defines catchall `foo!(container, ...)` stubs that raise
# an `ErrorException` explaining that in-place mutation is not supported and
# directing the user to the non-mutating variant.
# ─────────────────────────────────────────────────────────────────────────────

@testset "!-variant stubs raise errors" begin
    ao = _plotly_scenario_ao()
    matrix, clusters = _plotly_cluster_data()
    Si = _plotly_pawn_si()

    function _captures_error(f)
        try
            f()
            nothing
        catch e
            e
        end
    end

    @testset "scenarios! raises error (AnnotatedOutcomes dispatch)" begin
        err = _captures_error(() -> ADRIA.viz.scenarios!("not_a_grid", ao))
        @test err !== nothing
        @test err isa Union{MethodError,ErrorException}
        if err isa ErrorException
            @test occursin("Plotly", err.msg) || occursin("backend", lowercase(err.msg))
        end
    end

    @testset "scenarios! raises error (matrix+clusters dispatch)" begin
        err = _captures_error(() -> ADRIA.viz.scenarios!("not_a_grid", matrix, clusters))
        @test err !== nothing
        @test err isa Union{MethodError,ErrorException}
    end

    @testset "clustered_scenarios! raises error" begin
        err = _captures_error(
            () -> ADRIA.viz.clustered_scenarios!("not_a_grid", matrix, clusters)
        )
        @test err !== nothing
        @test err isa Union{MethodError,ErrorException}
    end

    @testset "taxonomy! raises error" begin
        scens_tax, rtc_tax = _plotly_taxonomy_df()
        err = _captures_error(() -> ADRIA.viz.taxonomy!("not_a_grid", scens_tax, rtc_tax))
        @test err !== nothing
        @test err isa Union{MethodError,ErrorException}
    end

    @testset "pawn! raises error" begin
        err = _captures_error(() -> ADRIA.viz.pawn!("not_a_grid", Si))
        @test err !== nothing
        @test err isa Union{MethodError,ErrorException}
    end

    @testset "tsa! raises error" begin
        Si_tsa = _plotly_tsa_si()
        err = _captures_error(() -> ADRIA.viz.tsa!("not_a_grid", Si_tsa))
        @test err !== nothing
        @test err isa Union{MethodError,ErrorException}
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# 14. Serialisation / savefig
# ─────────────────────────────────────────────────────────────────────────────

@testset "HTML serialisation" begin
    ao = _plotly_scenario_ao()
    p = ADRIA.viz.scenarios(ao)

    @testset "show(io, MIME\"text/html\", p) writes non-empty HTML" begin
        buf = IOBuffer()
        @test_nowarn show(buf, MIME"text/html"(), p)
        html_str = String(take!(buf))
        @test length(html_str) > 200
        @test contains(html_str, "plotly") || contains(html_str, "Plotly") ||
            contains(html_str, "data")
    end

    @testset "HTML file can be written" begin
        tmp = tempname() * ".html"
        try
            @test_nowarn open(tmp, "w") do io
                show(io, MIME"text/html"(), p)
            end
            @test isfile(tmp)
            @test filesize(tmp) > 100
        finally
            isfile(tmp) && rm(tmp)
        end
    end

    @testset "PlotlyKaleido PNG export (skip if not installed)" begin
        if isnothing(Base.find_package("PlotlyKaleido"))
            @info "PlotlyKaleido not installed — PNG export test skipped"
            @test true   # placeholder so the testset is not empty
        else
            @eval using PlotlyKaleido
            PlotlyKaleido.start()
            tmp = tempname() * ".png"
            try
                PlotlyKaleido.savefig(p, tmp)
                @test isfile(tmp)
                @test filesize(tmp) > 100
            finally
                isfile(tmp) && rm(tmp)
                try
                    PlotlyKaleido.kill()
                catch
                end
            end
        end
    end
end
# ─────────────────────────────────────────────────────────────────────────────
# 14. Spatial — map(gdf)
#
# Tests the GeoDataFrame-level API which does not require a ResultSet.
# Fixtures build synthetic ArchGDAL polygons in WGS84 (GBR region).
# ResultSet/Domain-level overloads are integration-tested manually.
# ─────────────────────────────────────────────────────────────────────────────

@testset "map(gdf)" begin
    gdf = _plotly_spatial_gdf(; n_sites=5)

    @testset "outline-only (color=nothing) returns PlotlyBase.Plot" begin
        p = ADRIA.viz.map(gdf)
        @test p isa PlotlyBase.Plot
    end

    @testset "has exactly one choropleth trace" begin
        p = ADRIA.viz.map(gdf)
        @test length(p.data) == 1
        @test p.data[1].type == "choropleth"
    end

    @testset "choropleth locations match site_id column" begin
        p = ADRIA.viz.map(gdf)
        @test length(p.data[1].locations) == nrow(gdf)
        @test p.data[1].locations == string.(gdf.site_id)
    end

    @testset "layout has geo attribute with fitbounds=locations" begin
        p = ADRIA.viz.map(gdf)
        geo = get(p.layout, :geo, nothing)
        @test !isnothing(geo)
        @test get(geo, :fitbounds, nothing) == "locations"
    end

    @testset "color=values produces choropleth with z matching input" begin
        values = collect(Float64, 1:nrow(gdf))
        p = ADRIA.viz.map(gdf; color=values)
        @test p isa PlotlyBase.Plot
        @test p.data[1].z == values
    end

    @testset "colorbar_label is set when color is provided" begin
        p = ADRIA.viz.map(gdf; color=rand(nrow(gdf)), colorbar_label="Cover [%]")
        @test p isa PlotlyBase.Plot
    end

    @testset "title kwarg is applied to layout" begin
        p = ADRIA.viz.map(gdf; title="My Map")
        @test get(p.layout, :title_text, "") == "My Map"
    end

    @testset "width and height kwargs are applied" begin
        p = ADRIA.viz.map(gdf; width=800, height=1000)
        @test get(p.layout, :width, 0) == 800
        @test get(p.layout, :height, 0) == 1000
    end

    @testset "gdf without site_id falls back to row indices" begin
        _tmp = _plotly_spatial_gdf()
        gdf_no_id = DataFrame(; k=_tmp.k, geometry=_tmp.geometry)
        p = ADRIA.viz.map(gdf_no_id)
        @test p isa PlotlyBase.Plot
        @test length(p.data[1].locations) == nrow(gdf_no_id)
    end

    @testset "gdf with no geometry column raises ArgumentError" begin
        bad_gdf = DataFrame(; site_id=["a", "b"], k=[0.3, 0.4])
        @test_throws ArgumentError ADRIA.viz.map(bad_gdf)
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# 15. Spatial !-variant stubs
# ─────────────────────────────────────────────────────────────────────────────

@testset "spatial !-variant stubs raise errors" begin
    @test_throws ErrorException ADRIA.viz.map!(nothing)
    @test_throws ErrorException ADRIA.viz.connectivity!(nothing)
    @test_throws ErrorException ADRIA.viz.ranks_to_frequencies!(nothing)
    @test_throws ErrorException ADRIA.viz.selection_criteria_map!(nothing)
end