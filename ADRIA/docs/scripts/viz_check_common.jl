# ADRIA/docs/scripts/viz_check_common.jl
#
# Shared setup and visualization checks for all ADRIA viz backends.
#
# Include this file from a backend-specific wrapper that first defines:
#
#   _activate_backend()            -> calls ADRIAviz.activate(...)
#   _save_fig(fig, name)           -> saves figure, returns path string
#   _viz_tsa(rs, si)               -> ADRIA.viz.tsa with backend signature
#   _viz_rsa(X, y, foi)            -> ADRIA.viz.rsa with backend signature
#                                     foi may be AbstractVector{Symbol} (multi-panel) or
#                                     NTuple{2,Symbol} (factor-vs-factor scatter)
#   _viz_outcome_map(X, y, foi)     -> ADRIA.viz.outcome_map with backend signature
#   _viz_rules_scatter(rs, X, tgt, rules) -> ADRIA.viz.rules_scatter with backend signature
#   _viz_backend_extras()          -> optional extra checks (no-op default)
#   RESULTS::Vector{NamedTuple} -> (name::String, ok::Bool, msg::String, elapsed::Float64)

# ----------------------------------------------------------------------------
# Shared packages
# ----------------------------------------------------------------------------

using Statistics
using DataFrames
using Dates
using Printf

using MLJ, SIRUS

@time using ADRIA
@time using ADRIAanalysis
@time using ADRIAviz

using ADRIA.analysis: cluster_scenarios, cluster_series
using ADRIAanalysis:
    find_scenarios, target_clusters, cluster_rules,
    data_envelopment_analysis, feature_set
using ADRIAanalysis.sensitivity: pawn, tsa, convergence, rsa

# Activate the backend (hook defined by the including file)
_activate_backend()

# ----------------------------------------------------------------------------
# check() — runs a viz closure, saves it, records pass/fail
# ----------------------------------------------------------------------------

function check(build_fn, name::String)
    println("\n>>> $name")
    start_time = time()
    try
        fig = build_fn()
        path = _save_fig(fig, name)
        elapsed = time() - start_time
        push!(RESULTS, (name=name, ok=true, msg=string(path), elapsed=elapsed))
        println("    [OK]   -> $path ($(Printf.@sprintf("%.2f", elapsed))s)")
        return fig
    catch err
        elapsed = time() - start_time
        msg = sprint(showerror, err)
        push!(RESULTS, (name=name, ok=false, msg=msg, elapsed=elapsed))
        println("    [FAIL] $msg ($(Printf.@sprintf("%.2f", elapsed))s)")
        println(sprint(showerror, err, catch_backtrace()))
        return nothing
    end
end

# ----------------------------------------------------------------------------
# Load domain, sample / reuse result set
# ----------------------------------------------------------------------------

dom_path = get(ENV, "ADRIA_TEST_DOMAIN", "")
if isempty(dom_path) || !isdir(dom_path)
    ArgumentError(
        """
        ADRIA_TEST_DOMAIN must be set to a valid domain directory path.
        Example: ADRIA_TEST_DOMAIN=/path/to/domain
        julia --project=. <script>"
        """
    )
end

@time example_dom = ADRIA.load_domain(dom_path, "45")

ADRIA.setup()
const OUTPUTS_DIR = get(ENV, "ADRIA_OUTPUT_DIR", "./Outputs")
const RS_PREFIX = "$(example_dom.name)__RCPs_45__"

existing_rs = if isdir(OUTPUTS_DIR)
    sort!(
        filter(readdir(OUTPUTS_DIR)) do d
            startswith(d, RS_PREFIX) && isdir(joinpath(OUTPUTS_DIR, d))
        end
    )
else
    String[]
end

rs = if !isempty(existing_rs)
    rs_path = joinpath(OUTPUTS_DIR, last(existing_rs))
    @info "Reusing existing result set (skipping run)" rs_path
    ADRIA.load_results(rs_path)
else
    n_scens = 2048
    @info "No existing result set found; sampling and running $(n_scens) scenarios."
    @info "Using $(Threads.nthreads()) threads."

    ADRIA.set_factor_bounds!(
        example_dom;
        N_seed_TA=(1000000.0, 15000000.0, 100000.0),
        N_seed_CA=(1000000.0, 15000000.0, 100000.0),
        N_seed_CNA=(1000000.0, 15000000.0, 100000.0),
        N_seed_SM=(1000000.0, 15000000.0, 100000.0),
        N_seed_LM=(1000000.0, 15000000.0, 100000.0)
    )

    run_scens = ADRIA.sample_set(example_dom, n_scens, "45")
    ADRIA.run_scenarios(example_dom, run_scens, "45")
end

scens = DataFrame(rs.inputs)

# ----------------------------------------------------------------------------
# Shared metrics / derived quantities
# ----------------------------------------------------------------------------

s_tac = ADRIA.metrics.scenario_total_cover(rs)
mean_s_tac = vec(mean(s_tac; dims=1))

# Feature set: rs.inputs enriched with DHW statistics and deployment summary stats,
# with N_seed_* replaced by actual deployment volume/location summaries.
# Used for RSA, outcome mapping, and rule induction.
fs = feature_set(rs)

# Outcome: number of years where mean relative cover exceeded 20%.
s_rc = ADRIA.metrics.scenario_relative_cover(rs)
y_outcome = collect(ADRIA.metrics.years_above_threshold(s_rc; threshold=0.20))

# Deployment-location-only outcome: same metric restricted to seeded locations.
deployed_locs = ADRIA.metrics.deployed_locations(rs; intervention=:seed)
s_rc_dep = ADRIA.metrics.scenario_relative_cover(rs; locations=deployed_locs)
y_dep = collect(ADRIA.metrics.years_above_threshold(s_rc_dep; threshold=0.20))

# Rank all features by Mann-Whitney RSA against the deployment-location outcome,
# then take the top 6 as the factors of interest. This replaces a hardcoded
# candidate list and ensures the most outcome-discriminating features drive all
# downstream visualizations (RSA, outcome map, rule induction).
rsa_si = rsa(fs, y_dep)
foi = rsa_si.feature[1:min(6, nrow(rsa_si))]

# ----------------------------------------------------------------------------
# 1. Scenario outcomes
# ----------------------------------------------------------------------------

check("scenarios_tac") do
    ADRIA.viz.scenarios(rs, s_tac)
end

# ----------------------------------------------------------------------------
# 2. Time series clustering
# ----------------------------------------------------------------------------

clusters4 = cluster_scenarios(s_tac, 4)

check("tsc") do
    ADRIA.viz.clustered_scenarios(s_tac, clusters4)
end

check("scenarios_by_cluster") do
    ADRIA.viz.scenarios(s_tac, clusters4)
end

# ----------------------------------------------------------------------------
# 3. Taxonomy (relative cover by taxa group)
# ----------------------------------------------------------------------------

check("taxonomy") do
    ADRIA.viz.taxonomy(rs)
end

# ----------------------------------------------------------------------------
# 4. PAWN sensitivity (heatmap)
# ----------------------------------------------------------------------------

check("pawn_si") do
    tac_Si = pawn(rs, mean_s_tac)
    ADRIA.viz.pawn(tac_Si)
end

# ----------------------------------------------------------------------------
# 5. Temporal sensitivity analysis  — backend-specific signature
# ----------------------------------------------------------------------------

check("tsa") do
    tsa_s = tsa(rs, s_tac)
    _viz_tsa(rs, tsa_s)
end

# ----------------------------------------------------------------------------
# 6. Convergence analysis
# ----------------------------------------------------------------------------

check("convergence_factors_series") do
    outcome = dropdims(mean(s_tac; dims=:timesteps); dims=:timesteps)
    # dhw_scenario is removed by feature_set; prefer dhw_mean as its replacement.
    # guided survives into fs unchanged.
    conv_foi = Symbol[f for f in (:dhw_mean, :guided) if string(f) in names(fs)]
    isempty(conv_foi) && (conv_foi = foi[1:min(2, length(foi))])
    Si_conv = convergence(fs, outcome, conv_foi)
    ADRIA.viz.convergence(Si_conv, conv_foi)
end

# ----------------------------------------------------------------------------
# 7. Regional Sensitivity Analysis  — backend-specific signature
# ----------------------------------------------------------------------------

check("rsa") do
    length(foi) < 2 && return nothing
    _viz_rsa(fs, y_dep, (foi[1], foi[2]))
end

check("rsa_multi_panel") do
    length(foi) < 4 && return nothing
    n = length(foi)
    # Span the importance ranking: pair top factor with 2nd and 3rd, then
    # pair mid-ranked factors together. Duplicates collapse if foi is short.
    panels = unique(
        NTuple{2,Symbol}[
            (foi[1], foi[2]),
            (foi[1], foi[min(3, n)]),
            (foi[2], foi[min(4, n)]),
            (foi[min(3, n)], foi[min(5, n)])
        ]
    )
    _viz_rsa(fs, y_dep, panels)
end

# ----------------------------------------------------------------------------
# 8. Outcome mapping  — backend-specific signature
# ----------------------------------------------------------------------------

check("outcome_map") do
    isempty(foi) && return nothing
    _viz_outcome_map(fs, y_outcome, foi)
end

check("outcome_map_deployed_locs") do
    isempty(deployed_locs) && return nothing
    isempty(foi) && return nothing
    _viz_outcome_map(fs, y_dep, foi)
end

# ----------------------------------------------------------------------------
# 9. Data Envelopment Analysis
# ----------------------------------------------------------------------------

check("example_dea_fig") do
    seed_cols = String[c for c in ("N_seed_TA", "N_seed_CA") if c in names(scens)]
    cost = if isempty(seed_cols)
        ones(Float64, nrow(scens))
    else
        Float64.(vec(sum(Matrix(scens[:, seed_cols]); dims=2))) .+ 1.0
    end

    s_tac_mean = dropdims(mean(s_tac; dims=:timesteps); dims=:timesteps)
    asv = ADRIA.metrics.absolute_shelter_volume(rs)
    s_sv = dropdims(
        mean(mean(asv; dims=:timesteps); dims=:locations);
        dims=(:timesteps, :locations)
    )

    function _norm01(v::AbstractVector{Float64})
        lo, hi = extrema(v)
        return hi - lo < eps() ? ones(length(v)) : (v .- lo) ./ (hi - lo)
    end

    X = _norm01(cost)
    Y = hcat(
        _norm01(Array{Float64}(s_tac_mean)),
        _norm01(Array{Float64}(s_sv))
    )
    DEA_out = data_envelopment_analysis(X, Y)
    ADRIA.viz.data_envelopment_analysis(rs, DEA_out)
end

# ----------------------------------------------------------------------------
# 10. Rule induction scatter  — backend-specific signature
# ----------------------------------------------------------------------------

check("rules_scatter") do
    clusters6 = cluster_scenarios(s_tac, 6)
    tgt = target_clusters(clusters6, s_tac)

    rule_foi = try
        # component_params returns raw input parameter names; filter to those
        # that survive feature_set post-processing (e.g. N_seed_* are removed).
        fs_cols = Set(names(fs))
        raw = ADRIA.component_params(rs, [Intervention, SeedCriteriaWeights]).fieldname
        Symbol[f for f in raw if string(f) in fs_cols]
    catch
        foi
    end
    isempty(rule_foi) && (rule_foi = foi)

    rules_iv = cluster_rules(rs, tgt, fs, rule_foi, 10; remove_duplicates=true)
    _viz_rules_scatter(rs, fs, tgt, rules_iv)
end

# ----------------------------------------------------------------------------
# 11. Spatial maps
# ----------------------------------------------------------------------------

check("tsc_map") do
    tac = ADRIA.metrics.total_absolute_cover(rs)
    tac_site_series = ADRIA.metrics.loc_trajectory(median, tac)
    tac_clusters = cluster_scenarios(tac_site_series, 6)
    tac_sites = ADRIA.metrics.per_loc(median, tac)
    _viz_tsc_map(rs, tac_sites, tac_clusters)
end

check("tsc_asv") do
    asv = ADRIA.metrics.absolute_shelter_volume(rs)
    asv_site_series = ADRIA.metrics.loc_trajectory(median, asv)
    asv_clusters = cluster_series(asv_site_series, 6)
    lowest = x -> x .∈ [sort(x; rev=true)[1:2]]
    asv_target = find_scenarios(asv_site_series, asv_clusters, lowest)
    ADRIA.viz.clustered_scenarios(asv_site_series, asv_target)
end

# ----------------------------------------------------------------------------
# 12. Location selection frequencies
# ----------------------------------------------------------------------------

const INTERVENTION_TYPES = (:seed, :fog, :shade, :mc)
_intervention_name(iv) = get(
    Dict(:seed => "Seed", :fog => "Fog", :shade => "Shade", :mc => "Moving Corals"),
    iv, titlecase(string(iv))
)

check("single_rank_plot") do
    seed_freq = ADRIA.decision.selection_frequency(rs.ranks, :seed)
    ADRIA.viz.map(rs, seed_freq)
end

check("ranks_by_intervention") do
    panels = Pair{String,Vector{Float64}}[]
    for iv in INTERVENTION_TYPES
        freq = try
            collect(Float64, ADRIA.decision.selection_frequency(rs.ranks, iv))
        catch
            continue
        end
        any(x -> isfinite(x) && x > 0, freq) || continue
        push!(panels, _intervention_name(iv) => freq)
    end
    isempty(panels) && error("No interventions with recorded selections")

    freq_matrix = reduce(hcat, last.(panels))
    ADRIA.viz.map(rs, freq_matrix, first.(panels))
end

# ----------------------------------------------------------------------------
# 13. Location selection criteria maps
# ----------------------------------------------------------------------------

check("criteria_spatial_plots") do
    mcda_funcs = ADRIA.decision.mcda_methods()

    guided_scens = ADRIA.sample_guided(example_dom, 2^2)
    scen = guided_scens[1, :]

    seed_pref = ADRIA.decision.SeedPreferences(example_dom, scen)

    sum_cover = vec(sum(example_dom.init_coral_cover; dims=1).data)
    dhw_scens = example_dom.dhw_scens[:, :, Int64(scen["dhw_scenario"])]
    plan_horizon = Int64(scen["plan_horizon"])
    decay = 0.99 .^ (1:(plan_horizon + 1)) .^ 2
    dhw_projection = ADRIA.decision.weighted_projection(
        dhw_scens, 1, plan_horizon, decay, 75
    )
    area_weighted_conn = example_dom.conn.data .* ADRIA.loc_k_area(example_dom)
    conn_cache = similar(area_weighted_conn)
    in_conn, out_conn, network = ADRIA.connectivity_strength(
        area_weighted_conn, sum_cover, conn_cache
    )

    seed_decision_mat = ADRIA.decision.decision_matrix(
        example_dom.loc_ids,
        seed_pref.names;
        seed_in_connectivity=in_conn,
        seed_out_connectivity=out_conn,
        seed_heat_stress=dhw_projection,
        seed_coral_cover=sum_cover
    )

    crit_agg = ADRIA.decision.criteria_aggregated_scores(
        seed_pref, seed_decision_mat, mcda_funcs[1]
    )

    is_const = Bool[length(x) == 1 for x in unique.(eachcol(seed_decision_mat.data))]

    ADRIA.viz.selection_criteria_map(
        example_dom,
        seed_decision_mat[criteria = .!is_const],
        crit_agg.scores ./ maximum(crit_agg.scores)
    )
end

# ----------------------------------------------------------------------------
# 14. Environmental drivers
# ----------------------------------------------------------------------------

check("dhw_scenario") do
    ADRIA.viz.dhw_scenario(example_dom, 1)
end

check("dhw_scenarios") do
    ADRIA.viz.dhw_scenarios(example_dom)
end

check("cyclone_scenario") do
    ADRIA.viz.cyclone_scenario(example_dom, 1)
end

# ----------------------------------------------------------------------------
# Backend-specific extras (no-op unless overridden)
# ----------------------------------------------------------------------------

_viz_backend_extras()

# ----------------------------------------------------------------------------
# 15. Connectivity graph
# ----------------------------------------------------------------------------

check("connectivity") do
    ADRIA.viz.connectivity(example_dom)
end

# ----------------------------------------------------------------------------
# Summary and report generation
# ----------------------------------------------------------------------------

function _format_time(elapsed::Float64)
    if elapsed < 1.0
        return "$(Printf.@sprintf("%.0f", elapsed * 1000))ms"
    else
        return "$(Printf.@sprintf("%.2f", elapsed))s"
    end
end

function _format_summary_line(result)
    status = result.ok ? "PASS" : "FAIL"
    time_str = _format_time(result.elapsed)
    detail = result.ok ? "" : "\n      Error: $(result.msg)"
    return "  [$status] $(result.name) [$time_str]$detail"
end

function _save_report_markdown(results, backend::String; filename="viz_check_report.md")
    report_name = replace(filename, r"\.md$" => "")
    filename = "$(backend)_$(report_name).md"

    n_ok = count(r -> r.ok, results)
    total_time = sum(r -> r.elapsed, results)

    lines = [
        "# Visualization Check Report",
        "",
        "**Date:** $(now())",
        "",
        "## Summary",
        "- **Total Checks:** $(length(results))",
        "- **Passed:** $n_ok",
        "- **Failed:** $(length(results) - n_ok)",
        "- **Total Time:** $(_format_time(total_time))",
        "",
        "## Results",
        ""
    ]

    for result in results
        status = result.ok ? "✓ PASS" : "✗ FAIL"
        time_str = _format_time(result.elapsed)
        push!(lines, "### $status — $(result.name) [$time_str]")
        if result.ok
            push!(lines, "")
            push!(lines, "Location: `$(result.msg)`")
        else
            push!(lines, "")
            push!(lines, "```")
            push!(lines, result.msg)
            push!(lines, "```")
        end
        push!(lines, "")
    end

    push!(lines, "## Timings by Check")
    push!(lines, "")
    sorted_by_time = sort(results; by=r -> r.elapsed, rev=true)
    for result in sorted_by_time
        time_str = _format_time(result.elapsed)
        status = result.ok ? "✓" : "✗"
        push!(lines, "- $status $(result.name): $time_str")
    end

    report_content = join(lines, "\n")
    write(filename, report_content)
    return filename
end

println("\n" * "="^70)
println("Visualization run summary")
println("="^70)
n_ok = count(r -> r.ok, RESULTS)
for result in RESULTS
    println(_format_summary_line(result))
end
println("-"^70)
total_time = sum(r -> r.elapsed, RESULTS)
println(
    "$(n_ok) / $(length(RESULTS)) visualizations succeeded in $(_format_time(total_time))"
)

# Save markdown report
report_path = _save_report_markdown(RESULTS, BACKEND_NAME)
println("\nReport saved to: $report_path")
