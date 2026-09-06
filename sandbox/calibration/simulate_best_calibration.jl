using Pkg
const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
Pkg.activate(joinpath(REPO_ROOT, "sandbox"))
cd(REPO_ROOT)

using ADRIA
using CSV, DataFrames, Statistics, Random, Dates

println("=== Unified ADRIA COTS Calibration Simulation Driver ===")

# 1. Load historical domain & reef mapping
println("Loading Lizard Island historical domain...")
dom = ADRIA.load_domain(ADRIA.LizardDomain, "sandbox/data/Lizard_Historical_v0.1", "historical")
site_to_reef = CSV.read("sandbox/data/Lizard_Historical_v0.1/site_to_reef.csv", DataFrame)
site_to_reef.reef_name_clean = [split(r, " (")[1] for r in site_to_reef.reef_name]
unique_reefs_clean = unique(site_to_reef.reef_name_clean)

# 2. Load best calibrated candidate from BBO summary
summary_file = joinpath(REPO_ROOT, "sandbox", "data", "bbo_best_summary.csv")
if !isfile(summary_file)
    summary_file = joinpath(REPO_ROOT, "sandbox", "data", "bbo_evaluated_candidates.csv")
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

println("Simulation Mode: ", N_scens > 1 ? "Stochastic Ensemble ($N_scens runs, mode=$stochastic_mode)" : "Deterministic Single-Run")

# Build scenario parameter template
scen_template = ADRIA.sample(dom, max(2, N_scens))[1:N_scens, :]
p_df = ADRIA.param_table(dom)

for col in names(scen_template)
    if N_scens == 1 || (!startswith(col, "surv_") && !in(col, ["dhw_scenario", "wave_scenario", "cyclone_mortality_scenario", "fecundity", "coral_recruitment"]))
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

if N_scens > 1 && stochastic_mode == "demographic"
    for s in 1:N_scens
        if hasproperty(scen_template, :a_ricker)
            scen_template.a_ricker[s] = clamp_param(best_row.a_ricker * lognormal_multiplier(rng, cots_demographic_cv), 2.0, 10.0)
        end
        if hasproperty(scen_template, :b_ricker)
            scen_template.b_ricker[s] = clamp_param(best_row.b_ricker * lognormal_multiplier(rng, cots_demographic_cv), 0.01, 0.5)
        end
    end
end

# 4. Execute Simulation Runs & Collect Trajectories
sim_df = DataFrame(
    sim_id = Int[],
    year = Int[],
    reef_name = String[],
    sim_cots_adult = Float64[],
    sim_cots_norm = Float64[],
    sim_coral_cover = Float64[]
)

site_df = DataFrame(
    sim_id = Int[],
    year = Int[],
    reef_name = String[],
    site_index = Int[],
    sim_cots_adult = Float64[],
    sim_coral_cover = Float64[]
)

for s in 1:N_scens
    if N_scens > 1
        println("Simulating run $s / $N_scens ...")
    end

    p = scen_template[s, :]
    rs = ADRIA.run_scenario(dom, p)

    adult_cots_site = rs.cots_log[:, 3, :]
    total_cover_site = dropdims(sum(rs.raw, dims=(2, 3)), dims=(2, 3))
    n_timesteps = size(adult_cots_site, 1)
    years = 1985:(1984 + n_timesteps)

    for reef in unique_reefs_clean
        site_indices = findall(site_to_reef.reef_name_clean .== reef)
        isempty(site_indices) && continue

        reef_sim_cots = [mean(adult_cots_site[t, site_indices]) for t in 1:n_timesteps]
        reef_sim_coral = [mean(total_cover_site[t, site_indices]) for t in 1:n_timesteps]

        for t in 1:n_timesteps
            push!(sim_df, (s, years[t], reef, reef_sim_cots[t], 0.0, reef_sim_coral[t]))
            for site_idx in site_indices
                push!(
                    site_df,
                    (s, years[t], reef, site_idx, adult_cots_site[t, site_idx], total_cover_site[t, site_idx])
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

# 5. Export Standardized Datasets
out_traj_path = joinpath(REPO_ROOT, "sandbox", "data", "best_calibrated_trajectories.csv")
out_site_path = joinpath(REPO_ROOT, "sandbox", "data", "best_calibrated_site_trajectories.csv")

CSV.write(out_traj_path, sim_df)
CSV.write(out_site_path, site_df)

println("Saved calibrated trajectories to: $out_traj_path")
println("Saved site-level trajectories to: $out_site_path")
println("=== Simulation completed successfully ===")
