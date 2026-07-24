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
using StatsBase: corspearman
using DataFrames
using Dates
using Printf

using MLJ, SIRUS

@time using ADRIA
@time using ADRIAanalysis
@time using ADRIAviz

const At = ADRIA.YAXArrays.DD.At

using ADRIA.analysis: cluster_scenarios, cluster_series
using ADRIAanalysis:
    find_scenarios, target_clusters, cluster_rules,
    data_envelopment_analysis, feature_set
using ADRIAanalysis.sensitivity: pawn, tsa, convergence, rsa, stratified_rsa

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
const OUTPUTS_DIR = get(ENV, "ADRIA_OUTPUT_DIR", joinpath(@__DIR__, "Outputs"))
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

# Counterfactual scenarios are embedded in rs and marked by guided == -1.
rs_cf = rs.inputs.guided .== -1
ug_mask = rs.inputs.guided .== 0
guided_mask = rs.inputs.guided .> 0   # MCDA-guided interventions only
iv_all_mask = .!rs_cf                  # guided + unguided (all non-CF)

scens = DataFrame(rs.inputs)

# ----------------------------------------------------------------------------
# Shared metrics / derived quantities
# ----------------------------------------------------------------------------

s_tac = ADRIA.metrics.scenario_total_cover(rs)
mean_s_tac = vec(mean(s_tac; dims=1))

# Deployed seeding locations.
# CF scenarios have rank 0, so the union below covers only guided+unguided sites.
deployed_locs = ADRIA.metrics.deployed_locations(rs; intervention=:seed)

# Guided-only deployed locations: exclude sites that only unguided scenarios ever
# used, so the metric for guided analysis isn't diluted by random-selection sites.
let _iv_ranks = Array(rs.ranks[intervention = At(:seed)])
    # dims: (timesteps × locations × scenarios)
    _guided_ever = vec(any(_iv_ranks[:, :, guided_mask] .> 0.0; dims=(1, 3)))
    global deployed_locs_guided = findall(_guided_ever)
end

# Counterfactual delta: per-scenario intervention effect (IV minus CF baseline).
# Outcome = years where relative cover at seeded locations exceeded 20%.
# Run once for ALL intervention scenarios (guided + unguided) so we can split
# afterwards; DHW stat columns are dropped from the returned feature matrix.
delta_all = ADRIAanalysis.sensitivity.counterfactual_delta(
    rs, rs_cf,
    function (_rs)
        src = ADRIA.metrics.scenario_relative_cover(_rs; locations=deployed_locs_guided)
        collect(ADRIA.metrics.years_above_threshold(src; threshold=0.20))
    end
)

# Post-filter to guided-only and unguided-only rows.
# iv_all_mask defines which rows of rs map to delta_all entries (rows ordered as
# they appear in rs, with CF rows removed).
guided_in_iv = guided_mask[iv_all_mask]   # Bool, length = sum(iv_all_mask)
ug_in_iv = ug_mask[iv_all_mask]

fs_guided = delta_all.fs_intervention[guided_in_iv, :]
y_guided = delta_all.y_delta[guided_in_iv]

fs_ug = delta_all.fs_intervention[ug_in_iv, :]
y_ug = delta_all.y_delta[ug_in_iv]

# Guided vs unguided lift: residual benefit of MCDA selection above random selection.
y_guided_vs_ug = y_guided .- mean(y_ug)

# Primary RSA analysis: guided scenarios vs CF at guided deployment locations.
fs = fs_guided
y_dep = y_guided

# Feature set retaining DHW stat columns for stratified_rsa (which stratifies on
# dhw_mean and drops DHW cols internally before each within-stratum RSA call).
# Must have the same row count as y_dep, so filter to guided_mask rows only.
fs_strat = feature_set(rs)[guided_mask, :]

# Rank features by the guided-vs-CF delta outcome; take top 6 as factors of interest.
rsa_si = rsa(fs, y_dep)
foi = rsa_si.feature[1:min(6, nrow(rsa_si))]

# ----------------------------------------------------------------------------
# 1. Scenario outcomes
# ----------------------------------------------------------------------------

check("scenarios_tac") do
    ADRIA.viz.scenarios(rs, s_tac)
end

# Delta time series: guided scenarios minus the mean CF trajectory.
check("iv_cf_delta_timeseries") do
    cf_mean_ts = vec(mean(Array(s_tac[scenarios = rs_cf]); dims=2))
    iv_arr = Array(s_tac[scenarios = guided_mask])   # (n_timesteps × n_guided)
    delta_arr = iv_arr .- cf_mean_ts                 # broadcast: subtract mean CF per timestep
    n_iv = sum(guided_mask)
    ts_labels = collect(s_tac.timesteps)
    delta_yax = ADRIA.DataCube(delta_arr; timesteps=ts_labels, scenarios=1:n_iv)
    guided_scens = scens[guided_mask, :]
    ADRIA.viz.scenarios(
        guided_scens, delta_yax;
        axis_opts=Dict{Symbol,Any}(:ylabel => "\u0394Total Cover (Guided \u2212 CF mean)")
    )
end

# Delta time series restricted to guided deployment locations only.
check("iv_cf_delta_timeseries_deployed") do
    s_tac_dep = ADRIA.metrics.scenario_total_cover(rs; locations=deployed_locs_guided)
    cf_mean_ts_dep = vec(mean(Array(s_tac_dep[scenarios = rs_cf]); dims=2))
    iv_arr_dep = Array(s_tac_dep[scenarios = guided_mask])
    delta_arr_dep = iv_arr_dep .- cf_mean_ts_dep
    n_iv = sum(guided_mask)
    ts_labels = collect(s_tac_dep.timesteps)
    delta_yax_dep = ADRIA.DataCube(delta_arr_dep; timesteps=ts_labels, scenarios=1:n_iv)
    guided_scens = scens[guided_mask, :]
    ADRIA.viz.scenarios(
        guided_scens, delta_yax_dep;
        axis_opts=Dict{Symbol,Any}(
            :ylabel => "\u0394Total Cover at Guided Deployment Locations (Guided \u2212 CF mean)"
        )
    )
end

# ----------------------------------------------------------------------------
# 2. Time series clustering
# ----------------------------------------------------------------------------

clusters4 = cluster_scenarios(s_tac, 4)

check("tsc") do
    ADRIA.viz.clustered_scenarios(s_tac, clusters4)
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
    tac_Si = pawn(feature_set(rs), mean_s_tac)
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
    # Filter outcome to intervention scenarios only (fs is already guided-only).
    outcome = dropdims(mean(s_tac; dims=:timesteps); dims=:timesteps)[guided_mask]
    # dhw_scenario is removed by feature_set; prefer dhw_mean as its replacement.
    # guided survives into fs unchanged.
    conv_foi = Symbol[
        f for f in (:dhw_mean, :guided, :mcda_method) if string(f) in names(fs)
    ]
    isempty(conv_foi) && (conv_foi = foi[1:min(2, length(foi))])
    Si_conv = convergence(fs, outcome, conv_foi)
    ADRIA.viz.convergence(Si_conv, conv_foi)
end

# ----------------------------------------------------------------------------
# 7. Regional Sensitivity Analysis  — backend-specific signature
# ----------------------------------------------------------------------------

check("rsa") do
    length(foi) < 2 && return nothing
    _viz_rsa(fs, y_dep, (foi[1], foi[2]); binary_mode=true, outcome_threshold=0.8)
end

check("rsa_multi_panel") do
    length(foi) < 2 && return nothing
    # Greedily select n_factors decorrelated factors from the RSA-ranked feature list,
    # then plot up to 9 pairwise combinations as factor-space panels
    # (x = factor A, y = factor B, colour = outcome).
    #
    # n_factors = 5  →  C(5,2) = 10 pairs; capped at 9 panels
    # cor_threshold  : skip a candidate if |ρ_Spearman| >= this with any selected factor
    # n_candidates   : pool size drawn from rsa_si; wider than foi so the greedy pass
    #                  has room to find 5 sufficiently decorrelated factors even when
    #                  the top-ranked features are mutually correlated.
    n_factors = 5
    cor_threshold = 0.9
    n_candidates = min(15, nrow(rsa_si))

    candidates = rsa_si.feature[1:n_candidates]
    cand_mat = Matrix{Float64}(fs[:, collect(candidates)])
    ρ = corspearman(cand_mat)   # (n_candidates × n_candidates)

    selected = Int[1]
    for i = 2:n_candidates
        if all(abs(ρ[i, j]) < cor_threshold for j in selected)
            push!(selected, i)
        end
        length(selected) >= n_factors && break
    end

    top = candidates[selected]
    panels = NTuple{2,Symbol}[
        (top[i], top[j]) for i = 1:length(top) for j = (i + 1):length(top)
    ]
    panels = first(panels, 9)
    isempty(panels) && return nothing
    _viz_rsa(fs, y_dep, panels; binary_mode=true, outcome_threshold=0.8)
end

# Classic RSA: Hornberger-Spear empirical CDF plots (one panel per factor).
# Uses the same top-ranked factors and delta outcome as the scatter RSA.
check("rsa_cdf") do
    isempty(foi) && return nothing
    ADRIA.viz.rsa_cdf(fs, y_dep, foi; outcome_threshold=0.8)
end

# ----------------------------------------------------------------------------
# 8. Outcome mapping  — backend-specific signature
# ----------------------------------------------------------------------------

check("outcome_map") do
    isempty(deployed_locs_guided) && return nothing
    isempty(foi) && return nothing
    _viz_outcome_map(fs, y_dep, foi)
end

# ----------------------------------------------------------------------------
# Task A -- DHW-stratified RSA (guided vs CF)
# ----------------------------------------------------------------------------

# fs_strat retains dhw_mean (required by stratified_rsa for binning); it matches
# y_dep row-for-row (both restricted to guided_mask scenarios).
check("rsa_stratified") do
    strat_result = stratified_rsa(fs_strat, y_dep)
    isempty(strat_result) && return nothing
    ADRIA.viz.stratified_rsa(strat_result)
end

# ----------------------------------------------------------------------------
# Task B -- Guided vs unguided RSA
# ----------------------------------------------------------------------------

# y_guided_vs_ug = guided delta minus mean(unguided delta): residual MCDA benefit.
# This comparison isolates MCDA selection quality from general intervention effects.

# Rank features by guided-vs-unguided lift; drop guided column (it is constant
# within this guided-only subset post-split, so dropping it is a no-op) —
# mcda_method is intentionally RETAINED here as the signal being ranked
# (which MCDA method drives the residual guided-vs-unguided benefit).
rsa_si_ug = let _fs = DataFrames.select(fs, Not(:guided))
    rsa(_fs, y_guided_vs_ug)
end
foi_ug = rsa_si_ug.feature[1:min(6, nrow(rsa_si_ug))]

check("rsa_guided_vs_unguided") do
    length(foi_ug) < 2 && return nothing
    _viz_rsa(
        DataFrames.select(fs, Not(:guided)), y_guided_vs_ug,
        (foi_ug[1], foi_ug[2]); binary_mode=true, outcome_threshold=0.8
    )
end

check("rsa_cdf_guided_vs_unguided") do
    isempty(foi_ug) && return nothing
    ADRIA.viz.rsa_cdf(
        DataFrames.select(fs, Not(:guided)), y_guided_vs_ug,
        foi_ug; outcome_threshold=0.8
    )
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
    # Delta TAC time series: guided minus CF (matched if equal counts, mean-CF otherwise).
    s_tac_iv = Array(s_tac[scenarios = guided_mask])
    s_tac_cf = Array(s_tac[scenarios = rs_cf])
    s_tac_delta =
        sum(guided_mask) == sum(rs_cf) ?
        s_tac_iv .- s_tac_cf :
        s_tac_iv .- mean(s_tac_cf; dims=2)

    clusters6 = cluster_scenarios(s_tac_delta, 6)
    tgt = target_clusters(clusters6, s_tac_delta)

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
    projection_confidence = scen["projection_confidence"]
    decay = ADRIA.decision.build_decay(plan_horizon, projection_confidence)
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
    filename = "$(report_name)_$(backend).md"

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
