using Pkg
const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
Pkg.activate(joinpath(REPO_ROOT, "sandbox"))

using CSV
using DataFrames
using Statistics
using TOML

include(joinpath(@__DIR__, "cots_cycle_metrics.jl"))

run_id = get(ENV, "BBO_RUN_ID", "expanded_cotsconn_pilot_seed20260930")
run_dir = joinpath(REPO_ROOT, "sandbox", "calibration", "runs", run_id)
trajectory_path = joinpath(run_dir, "best_calibrated_trajectories.csv")
observation_path = joinpath(REPO_ROOT, "sandbox", "data", "reef_cots.csv")
isfile(trajectory_path) || error("Missing trajectories: $trajectory_path")

trajectories = CSV.read(trajectory_path, DataFrame)
observations = CSV.read(observation_path, DataFrame)
split_year = parse(Int, get(ENV, "COTS_TEMPORAL_SPLIT_YEAR", "2005"))

sim_to_obs = Dict(
    "Lizard Island Reef" => "Lizard Isles",
    "MacGillivray Reef" => "Macgillivray Reef",
    "North Direction Reef" => "North Direction Island",
    "Eyrie Reef" => "Eyrie Reef"
)

results = DataFrame(
    reef_name=String[],
    training_end_year=Int[],
    validation_start_year=Int[],
    observation_scale=Float64[],
    validation_loss=Float64[],
    validation_rmse=Float64[],
    validation_abs_percent_bias=Float64[],
    validation_peak_count_penalty=Float64[],
    validation_peak_timing_penalty=Float64[],
    validation_amplitude_penalty=Float64[],
    validation_peak_height_penalty=Float64[],
    validation_peak_prominence_penalty=Float64[],
    validation_matched_peaks=Int[],
    validation_sim_peak_years=String[],
    validation_obs_peak_years=String[],
    passed_peak_gate=Bool[]
)

for (sim_reef, obs_reef) in sim_to_obs
    sim = trajectories[trajectories.reef_name .== sim_reef, :]
    obs = observations[observations.reef_name .== obs_reef, :]
    sort!(sim, :year)
    sort!(obs, :year)

    train_sim = sim[sim.year .<= split_year, :]
    train_obs = obs[obs.year .<= split_year, :]
    training_score = cots_cycle_score(
        train_sim.year, train_sim.sim_cots_adult,
        train_obs.year, train_obs.cotsptow
    )

    validation_sim = sim[sim.year .> split_year, :]
    validation_obs = obs[obs.year .> split_year, :]
    validation_score = cots_cycle_score(
        validation_sim.year, validation_sim.sim_cots_adult,
        validation_obs.year, validation_obs.cotsptow;
        fixed_observation_scale=training_score.observation_scale
    )
    peak_gate = validation_score.n_sim_peaks == validation_score.n_obs_peaks &&
        validation_score.n_matched_peaks == validation_score.n_obs_peaks &&
        validation_score.peak_timing_penalty <= 0.5 &&
        validation_score.peak_height_penalty <= 0.30 &&
        validation_score.peak_prominence_penalty <= 0.35

    push!(results, (
        sim_reef,
        split_year,
        split_year + 1,
        training_score.observation_scale,
        validation_score.total_loss,
        validation_score.matched_rmse,
        validation_score.percent_bias_abs,
        validation_score.peak_count_penalty,
        validation_score.peak_timing_penalty,
        validation_score.amplitude_penalty,
        validation_score.peak_height_penalty,
        validation_score.peak_prominence_penalty,
        validation_score.n_matched_peaks,
        join(validation_score.sim_peak_years, ";"),
        join(validation_score.obs_peak_years, ";"),
        peak_gate
    ))
end

output_path = joinpath(run_dir, "temporal_holdout_validation.csv")
CSV.write(output_path, results)

summary = Dict(
    "run_id" => run_id,
    "training_window" => "through $split_year",
    "validation_window" => "$(split_year + 1) onward",
    "reefs_evaluated" => nrow(results),
    "reefs_passing_peak_gate" => count(results.passed_peak_gate),
    "mean_validation_loss" => mean(results.validation_loss),
    "mean_validation_amplitude_penalty" => mean(results.validation_amplitude_penalty),
    "mean_validation_matched_peaks" => mean(results.validation_matched_peaks),
    "passed" => all(results.passed_peak_gate)
)
open(joinpath(run_dir, "temporal_holdout_summary.toml"), "w") do io
    TOML.print(io, summary; sorted=true)
end

show(stdout, MIME("text/plain"), results)
println()
println("Wrote $output_path")
