# ADRIA/docs/scripts/viz_check_common.jl
#
# Shared setup and visualization checks for all ADRIA viz backends.
#
# Include this file from a backend-specific wrapper that first defines:
#
#   _activate_backend()            -> calls ADRIAviz.activate(...)
#   _save_fig(fig, name)           -> saves figure, returns path string
#   _viz_tsa(rs, si)               -> ADRIA.viz.tsa with backend signature
#   _viz_rsa(rs, rsa_ds, foi)      -> ADRIA.viz.rsa with backend signature
#   _viz_outcome_map(rs, om_ds, foi) -> ADRIA.viz.outcome_map with backend signature
#   _viz_rules_scatter(rs, scens, tgt, rules) -> ADRIA.viz.rules_scatter with backend signature
#   _viz_backend_extras()          -> optional extra checks (no-op default)
#   RESULTS::Vector{Tuple{String,Bool,String}}

# ----------------------------------------------------------------------------
# Shared packages
# ----------------------------------------------------------------------------

using Statistics
using DataFrames

using MLJ, SIRUS

@time using ADRIA
@time using ADRIAanalysis
@time using ADRIAviz

using ADRIAanalysis:
    find_scenarios, target_clusters, cluster_rules,
    data_envelopment_analysis
using ADRIAanalysis.sensitivity: pawn, tsa, rsa, outcome_map, convergence
using ADRIA.analysis: cluster_scenarios, cluster_series

# Activate the backend (hook defined by the including file)
_activate_backend()

# ----------------------------------------------------------------------------
# check() — runs a viz closure, saves it, records pass/fail
# ----------------------------------------------------------------------------

function check(build_fn, name::String)
    println("\n>>> $name")
    try
        fig = build_fn()
        path = _save_fig(fig, name)
        push!(RESULTS, (name, true, string(path)))
        println("    [OK]   -> $path")
        return fig
    catch err
        msg = sprint(showerror, err)
        push!(RESULTS, (name, false, msg))
        println("    [FAIL] $msg")
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
    @info "No existing result set found; sampling and running 128 scenarios."
    run_scens = ADRIA.sample_set(example_dom, 128, "45")
    ADRIA.run_scenarios(example_dom, run_scens, "45")
end

scens = DataFrame(rs.inputs)

# ----------------------------------------------------------------------------
# Shared metrics / derived quantities
# ----------------------------------------------------------------------------

s_tac = ADRIA.metrics.scenario_total_cover(rs)
mean_s_tac = vec(mean(s_tac; dims=1))

const CANDIDATE_FOI = [:N_seed_TA, :N_seed_CA, :fogging, :SRM, :a_adapt]
foi = Symbol[f for f in CANDIDATE_FOI if string(f) in names(rs.inputs)]

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
    conv_foi = Symbol[f for f in (:dhw_scenario, :guided) if string(f) in names(rs.inputs)]
    isempty(conv_foi) && (conv_foi = foi[1:min(2, length(foi))])
    Si_conv = convergence(scens, outcome, conv_foi)
    ADRIA.viz.convergence(Si_conv, conv_foi)
end

# ----------------------------------------------------------------------------
# 7. Regional Sensitivity Analysis  — backend-specific signature
# ----------------------------------------------------------------------------

check("rsa") do
    rsa_ds = rsa(rs, mean_s_tac, foi; S=10)
    _viz_rsa(rs, rsa_ds, foi)
end

# ----------------------------------------------------------------------------
# 8. Outcome mapping  — backend-specific signature
# ----------------------------------------------------------------------------

check("outcome_map") do
    om_ds = outcome_map(rs, mean_s_tac, x -> any(x .>= 0.5), foi; S=20)
    _viz_outcome_map(rs, om_ds, foi)
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
        ADRIA.component_params(rs, [Intervention, SeedCriteriaWeights]).fieldname
    catch
        foi
    end

    rules_iv = cluster_rules(rs, tgt, scens, rule_foi, 10; remove_duplicates=true)
    _viz_rules_scatter(rs, scens, tgt, rules_iv)
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
# 12. Connectivity graph
# ----------------------------------------------------------------------------

check("connectivity") do
    ADRIA.viz.connectivity(example_dom)
end

# ----------------------------------------------------------------------------
# 13. Location selection frequencies
# ----------------------------------------------------------------------------

const INTERVENTION_TYPES = (:seed, :fog, :shade, :mc)
_intervention_name(iv) = get(
    Dict(:seed => "Seed", :fog => "Fog", :shade => "Shade", :mc => "Move Corals"),
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
# 14. Location selection criteria maps
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
# 15. Environmental drivers
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
# Summary
# ----------------------------------------------------------------------------

println("\n" * "="^70)
println("Visualization run summary")
println("="^70)
n_ok = count(r -> r[2], RESULTS)
for (name, ok, msg) in RESULTS
    status = ok ? "OK  " : "FAIL"
    detail = ok ? "" : "  ($msg)"
    println("  [$status] $name$detail")
end
println("-"^70)
println("$(n_ok) / $(length(RESULTS)) visualizations succeeded.")
