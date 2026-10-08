using Pkg
const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
Pkg.activate(joinpath(REPO_ROOT, "sandbox"))
cd(REPO_ROOT)

using ADRIA
using CSV, DataFrames, Statistics, Random, Dates

println("=== Unified ADRIA COTS Calibration Simulation Driver ===")
ENV["COTS_ALLEE_THRESHOLD"] = get(ENV, "COTS_ALLEE_THRESHOLD", "3.0")

# 1. Load historical domain & reef mapping
println("Loading Lizard Island historical domain...")
domain_version = get(ENV, "LIZARD_DOMAIN_VERSION", "Lizard_Historical_v0.1")
occursin(r"^[A-Za-z0-9_.-]+$", domain_version) || error("Unsupported LIZARD_DOMAIN_VERSION")
domain_path = joinpath(REPO_ROOT, "sandbox", "data", domain_version)
dom = ADRIA.load_domain(ADRIA.LizardDomain, domain_path, "historical")
site_to_reef = CSV.read(joinpath(domain_path, "site_to_reef.csv"), DataFrame)
site_to_reef.reef_name_clean = [split(r, " (")[1] for r in site_to_reef.reef_name]
unique_reefs_clean = unique(site_to_reef.reef_name_clean)

# 2. Load best calibrated candidate from BBO summary
requested_run_dir = get(ENV, "BBO_RUN_DIR", "")
if isempty(requested_run_dir) && haskey(ENV, "BBO_RUN_ID")
    requested_run_dir = joinpath(REPO_ROOT, "sandbox", "calibration", "runs", ENV["BBO_RUN_ID"])
end
artifact_dir = isempty(requested_run_dir) ? joinpath(REPO_ROOT, "sandbox", "data") : normpath(requested_run_dir)
summary_file = isempty(requested_run_dir) ?
    joinpath(artifact_dir, "bbo_best_summary.csv") : joinpath(artifact_dir, "best_summary.csv")
if !isfile(summary_file)
    summary_file = isempty(requested_run_dir) ?
        joinpath(artifact_dir, "bbo_evaluated_candidates.csv") : joinpath(artifact_dir, "evaluated_candidates.csv")
end

if !isfile(summary_file)
    error("No optimization results found at $summary_file. Run calibrate_cots_blackbox.jl first.")
end

best_summary_df = CSV.read(summary_file, DataFrame)
sort!(best_summary_df, :loss)
best_row = best_summary_df[1, :]

println("Loaded Best Candidate ID: ", best_row.eval_id)
println("Best Candidate Loss Score: ", round(best_row.loss; digits=4))

# 3. Determine Execution Mode (Deterministic vs Stochastic Ensemble)
N_scens = parse(Int, get(ENV, "COTS_N_STOCHASTIC_SCENS", "1"))
stochastic_mode = lowercase(get(ENV, "COTS_STOCHASTIC_MODE", "demographic"))
stochastic_seed = parse(Int, get(ENV, "COTS_STOCHASTIC_SEED", "20260902"))
stochastic_mode in ["environmental", "demographic", "combined"] || error("Unknown COTS_STOCHASTIC_MODE: $stochastic_mode")
Random.seed!(stochastic_seed)

println("Simulation Mode: ", N_scens > 1 ? "Stochastic Ensemble ($N_scens runs, mode=$stochastic_mode)" : "Deterministic Single-Run")

# Build scenario parameter template
scen_template = ADRIA.sample(dom, max(2, N_scens))[1:N_scens, :]
p_df = ADRIA.param_table(dom)

for col in names(scen_template)
    preserve_environment = N_scens > 1 && stochastic_mode in ["environmental", "combined"] &&
        (startswith(col, "surv_") || col in ["dhw_scenario", "wave_scenario", "cyclone_mortality_scenario", "fecundity", "coral_recruitment"])
    if !preserve_environment
        scen_template[!, col] .= p_df[1, col]
    end
end

# Exclude metadata columns
meta_cols = [
    "eval_id", "timestamp", "loss", "cycle_loss", "legacy_loss", "mean_rmse",
    "mean_abs_percent_bias", "mean_pearson", "mean_spearman", "mean_peak_count_penalty",
    "mean_peak_timing_penalty", "mean_period_penalty", "mean_amplitude_penalty",
    "mean_flatline_penalty", "mean_lag_correlation_penalty", "mean_best_lag_years",
    "mean_best_lag_pearson", "mean_best_lag_spearman"
]

# Apply optimal parameter values across scenario table
for col in names(best_summary_df)
    if col in meta_cols
        continue
    end
    if hasproperty(scen_template, Symbol(col))
        scen_template[!, Symbol(col)] .= best_row[col]
    end
end
scen_template[!, :allee_threshold] .= parse(Float64, ENV["COTS_ALLEE_THRESHOLD"])

# Handle seed multiplier & pulse environment variables
if hasproperty(best_row, :seed_mult)
    ENV["COTS_INITIAL_MULTIPLIER"] = string(best_row.seed_mult)
end

if hasproperty(best_row, :pulse_relative_magnitude) && best_row.pulse_relative_magnitude > 0.05
    ENV["COTS_EXTERNAL_PULSE"] = "true"
    ENV["COTS_PULSE_START"] = string(round(Int, best_row.pulse_start))
    ENV["COTS_PULSE_DURATION"] = string(round(Int, best_row.pulse_duration))
    ENV["COTS_PULSE_RELATIVE_MAGNITUDE"] = string(best_row.pulse_relative_magnitude)
else
    ENV["COTS_EXTERNAL_PULSE"] = "false"
end

# Stochastic jitter setup if N_scens > 1
rng = MersenneTwister(stochastic_seed)
lognormal_multiplier(rng::AbstractRNG, cv::Float64) = exp(randn(rng) * cv)
clamp_param(x::Float64, lo::Float64, hi::Float64) = min(max(x, lo), hi)
cots_demographic_cv = parse(Float64, get(ENV, "COTS_DEMOGRAPHIC_CV", "0.08"))
cots_seed_cv = parse(Float64, get(ENV, "COTS_SEED_MULT_CV", "0.05"))
cots_pulse_cv = parse(Float64, get(ENV, "COTS_PULSE_MAGNITUDE_CV", "0.05"))
base_seed_multiplier = hasproperty(best_row, :seed_mult) ? Float64(best_row.seed_mult) : 1.0
base_pulse_magnitude = hasproperty(best_row, :pulse_relative_magnitude) ? Float64(best_row.pulse_relative_magnitude) : 0.0
seed_multipliers = fill(base_seed_multiplier, N_scens)
pulse_magnitudes = fill(base_pulse_magnitude, N_scens)

jitter_specs = [
    (:a_ricker, 2.0, 10.0), (:b_ricker, 0.01, 0.5),
    (:m1, 0.1, 0.9), (:m2, 0.05, 0.5), (:m3, 0.05, 0.3),
    (:p_tilde, 0.8, 1.0), (:C_max, 0.4, 1.0),
    (:tau_condition, 1.0, 10.0),
    (:imm_threshold, 0.1, 0.8), (:eta_imm, 1.0, 5.0),
]
if N_scens > 1 && stochastic_mode in ["demographic", "combined"]
    for s in 1:N_scens
        for (name, lower, upper) in jitter_specs
            hasproperty(scen_template, name) || continue
            base_value = hasproperty(best_row, name) ? Float64(best_row[name]) : Float64(scen_template[1, name])
            scen_template[s, name] = clamp_param(base_value * lognormal_multiplier(rng, cots_demographic_cv), lower, upper)
        end
        seed_multipliers[s] = base_seed_multiplier * lognormal_multiplier(rng, cots_seed_cv)
        pulse_magnitudes[s] = base_pulse_magnitude * lognormal_multiplier(rng, cots_pulse_cv)
    end
end

simulation_metadata = copy(scen_template)
insertcols!(simulation_metadata, 1, :sim_id => collect(1:N_scens))
simulation_metadata[!, :cots_initial_multiplier] = seed_multipliers
simulation_metadata[!, :cots_pulse_relative_magnitude] = pulse_magnitudes
simulation_metadata[!, :stochastic_mode] = fill(stochastic_mode, N_scens)
simulation_metadata[!, :stochastic_seed] = fill(stochastic_seed, N_scens)
simulation_metadata[!, :cots_connectivity_mode] = fill(
    lowercase(get(ENV, "ADRIA_COTS_CONNECTIVITY_MODE", "auto")), N_scens
)
for key in [
    "COTS_ALLEE_THRESHOLD", "COTS_CONNECTIVITY_TEMPORAL_MODE",
    "COTS_CONNECTIVITY_SEED", "COTS_APPLY_LARVAL_SURVIVAL",
    "COTS_EXTERNAL_SOURCE_DENSITY", "COTS_JUVENILE_STORAGE",
    "COTS_HABITAT_MEDIATION", "COTS_CALENDAR_START_YEAR",
    "COTS_SEPARATED_RECRUITMENT", "COTS_LARVAL_FECUNDITY",
    "COTS_SETTLEMENT_PROBABILITY", "COTS_EXTERNAL_OUTBREAK_PRODUCTION"
]
    simulation_metadata[!, Symbol(lowercase(key))] = fill(get(ENV, key, ""), N_scens)
end

# 4. Execute Simulation Runs & Collect Trajectories
sim_df = DataFrame(
    sim_id = Int[],
    year = Int[],
    reef_name = String[],
    sim_cots_recruits = Float64[],
    sim_cots_juveniles = Float64[],
    sim_cots_adult = Float64[],
    sim_cots_condition = Float64[],
    sim_cots_norm = Float64[],
    sim_coral_cover = Float64[]
)

site_df = DataFrame(
    sim_id = Int[],
    year = Int[],
    reef_name = String[],
    site_index = Int[],
    sim_cots_recruits = Float64[],
    sim_cots_juveniles = Float64[],
    sim_cots_adult = Float64[],
    sim_cots_condition = Float64[],
    sim_coral_cover = Float64[]
)

flow_df = DataFrame(
    sim_id=Int[], year=Int[], forcing_year=Int[], reef_name=String[], site_index=Int[],
    local_fecundity=Float64[], background_immigration=Float64[],
    internal_immigration=Float64[], external_immigration=Float64[],
    maturation=Float64[], retained_juveniles=Float64[], settlement_gate=Float64[],
    local_retention=Float64[], pelagic_survivors=Float64[],
    settled_recruits=Float64[], external_potential=Float64[],
    external_pelagic=Float64[]
)

for s in 1:N_scens
    if N_scens > 1
        println("Simulating run $s / $N_scens ...")
    end

    ENV["COTS_INITIAL_MULTIPLIER"] = string(seed_multipliers[s])
    if pulse_magnitudes[s] > 0.05
        ENV["COTS_EXTERNAL_PULSE"] = "true"
        ENV["COTS_PULSE_RELATIVE_MAGNITUDE"] = string(pulse_magnitudes[s])
    else
        ENV["COTS_EXTERNAL_PULSE"] = "false"
    end
    p = scen_template[s, :]
    rs = ADRIA.run_scenario(dom, p)

    recruit_cots_site = rs.cots_log[:, 1, :]
    juvenile_cots_site = rs.cots_log[:, 2, :]
    adult_cots_site = rs.cots_log[:, 3, :]
    condition_cots_site = rs.cots_condition_log
    cots_flow_site = rs.cots_flow_log
    total_cover_site = dropdims(sum(rs.raw, dims=(2, 3)), dims=(2, 3))
    n_timesteps = size(adult_cots_site, 1)
    years = 1985:(1984 + n_timesteps)

    for reef in unique_reefs_clean
        site_indices = findall(site_to_reef.reef_name_clean .== reef)
        isempty(site_indices) && continue

        reef_sim_recruits = [mean(recruit_cots_site[t, site_indices]) for t in 1:n_timesteps]
        reef_sim_juveniles = [mean(juvenile_cots_site[t, site_indices]) for t in 1:n_timesteps]
        reef_sim_cots = [mean(adult_cots_site[t, site_indices]) for t in 1:n_timesteps]
        reef_sim_condition = [mean(condition_cots_site[t, site_indices]) for t in 1:n_timesteps]
        reef_sim_coral = [mean(total_cover_site[t, site_indices]) for t in 1:n_timesteps]

        for t in 1:n_timesteps
            push!(sim_df, (
                s, years[t], reef,
                reef_sim_recruits[t], reef_sim_juveniles[t], reef_sim_cots[t],
                reef_sim_condition[t], 0.0, reef_sim_coral[t]
            ))
            for site_idx in site_indices
                push!(
                    site_df,
                    (
                        s, years[t], reef, site_idx,
                        recruit_cots_site[t, site_idx], juvenile_cots_site[t, site_idx],
                        adult_cots_site[t, site_idx], condition_cots_site[t, site_idx],
                        total_cover_site[t, site_idx]
                    )
                )
                push!(
                    flow_df,
                    (
                        s, years[t], rs.cots_forcing_year[t], reef, site_idx,
                        cots_flow_site[t, 1, site_idx], cots_flow_site[t, 2, site_idx],
                        cots_flow_site[t, 3, site_idx], cots_flow_site[t, 4, site_idx],
                        cots_flow_site[t, 5, site_idx], cots_flow_site[t, 6, site_idx],
                        cots_flow_site[t, 7, site_idx], cots_flow_site[t, 8, site_idx],
                        cots_flow_site[t, 9, site_idx], cots_flow_site[t, 10, site_idx],
                        cots_flow_site[t, 11, site_idx], cots_flow_site[t, 12, site_idx]
                    )
                )
            end
        end
    end
end

# Normalize adult COTS per reef across simulation runs
for reef in unique(sim_df.reef_name)
    reef_rows = sim_df.reef_name .== reef
    max_sim = maximum(sim_df[reef_rows, :sim_cots_adult])
    if max_sim > 0
        sim_df[reef_rows, :sim_cots_norm] .= sim_df[reef_rows, :sim_cots_adult] ./ max_sim
    end
end

# Apply the fitted calibration observation model for direct COTS-per-tow plots.
sim_df[!, :sim_cots_cpue] = fill(NaN, nrow(sim_df))
by_reef_file = joinpath(artifact_dir, "evaluated_by_reef.csv")
if isfile(by_reef_file)
    by_reef_scores = CSV.read(by_reef_file, DataFrame)
    best_scores = by_reef_scores[by_reef_scores.eval_id .== best_row.eval_id, :]
    for row in eachrow(best_scores)
        reef_rows = sim_df.reef_name .== row.reef_name
        sim_df[reef_rows, :sim_cots_cpue] .= sim_df[reef_rows, :sim_cots_adult] .* row.observation_scale
    end
end

# 5. Export Standardized Datasets
out_traj_path = joinpath(artifact_dir, "best_calibrated_trajectories.csv")
out_site_path = joinpath(artifact_dir, "best_calibrated_site_trajectories.csv")
out_flow_path = joinpath(artifact_dir, "best_calibrated_cots_flows.csv")
out_metadata_path = joinpath(artifact_dir, "simulation_metadata.csv")

CSV.write(out_traj_path, sim_df)
CSV.write(out_site_path, site_df)
CSV.write(out_flow_path, flow_df)
CSV.write(out_metadata_path, simulation_metadata)

println("Saved calibrated trajectories to: $out_traj_path")
println("Saved site-level trajectories to: $out_site_path")
println("Saved decomposed COTS flows to: $out_flow_path")
println("Saved simulation metadata to: $out_metadata_path")
println("=== Simulation completed successfully ===")
