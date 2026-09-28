using Pkg
const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
Pkg.activate(joinpath(REPO_ROOT, "sandbox"))

using CSV
using DataFrames
using Statistics
using TOML

include(joinpath(@__DIR__, "cots_cycle_metrics.jl"))

coral_run_id = get(ENV, "BBO_CORAL_RUN_ID", "expanded_peak_pilot_seed20260930")
cots_run_id = get(ENV, "BBO_COTS_RUN_ID", "expanded_cotsconn_pilot_seed20260930")
runs_root = joinpath(REPO_ROOT, "sandbox", "calibration", "runs")
coral_dir = joinpath(runs_root, coral_run_id)
cots_dir = joinpath(runs_root, cots_run_id)

coral = CSV.read(joinpath(coral_dir, "evaluated_candidates.csv"), DataFrame)
cots = CSV.read(joinpath(cots_dir, "evaluated_candidates.csv"), DataFrame)
coral = coral[coral.status .== "success", :]
cots = cots[cots.status .== "success", :]
sort!(coral, :eval_id)
sort!(cots, :eval_id)

coral_metadata = TOML.parsefile(joinpath(coral_dir, "run_metadata.toml"))
cots_metadata = TOML.parsefile(joinpath(cots_dir, "run_metadata.toml"))
param_names = String.(coral_metadata["search"]["parameter_names"])
param_names == String.(cots_metadata["search"]["parameter_names"]) ||
    error("Connectivity runs used different parameter contracts")
coral.eval_id == cots.eval_id || error("Connectivity runs evaluated different candidate IDs")
for name in param_names
    all(isapprox.(Float64.(coral[!, name]), Float64.(cots[!, name]); atol=1e-12, rtol=1e-12)) ||
        error("Candidate values differ for $name")
end

metrics = [
    "loss",
    "mean_rmse",
    "mean_abs_percent_bias",
    "mean_peak_count_penalty",
    "mean_peak_timing_penalty",
    "mean_period_penalty",
    "mean_amplitude_penalty",
    "mean_peak_height_penalty",
    "mean_peak_prominence_penalty",
    "mean_matched_peaks"
]

paired = DataFrame(eval_id=Int.(coral.eval_id))
for metric in metrics
    coral_values = Float64.(coral[!, metric])
    cots_values = Float64.(cots[!, metric])
    paired[!, Symbol("coral_" * metric)] = coral_values
    paired[!, Symbol("cots_" * metric)] = cots_values
    paired[!, Symbol("delta_" * metric)] = cots_values .- coral_values
end
CSV.write(joinpath(cots_dir, "connectivity_factorial_comparison.csv"), paired)

outcomes = [
    "loss",
    "mean_peak_count_penalty",
    "mean_peak_timing_penalty",
    "mean_amplitude_penalty",
    "mean_peak_height_penalty",
    "mean_peak_prominence_penalty"
]
screen = DataFrame(
    parameter=String[],
    outcome=String[],
    coral_spearman=Float64[],
    cots_spearman=Float64[],
    direction_stable=Bool[],
    mean_abs_spearman=Float64[]
)
for parameter in param_names
    x = Float64.(coral[!, parameter])
    x_rank = rank_average(x)
    for outcome in outcomes
        coral_rho = safe_cor(x_rank, rank_average(Float64.(coral[!, outcome])))
        cots_rho = safe_cor(x_rank, rank_average(Float64.(cots[!, outcome])))
        stable = isfinite(coral_rho) && isfinite(cots_rho) &&
            sign(coral_rho) == sign(cots_rho)
        mean_abs = mean(abs.([coral_rho, cots_rho]))
        push!(screen, (parameter, outcome, coral_rho, cots_rho, stable, mean_abs))
    end
end
sort!(screen, [:outcome, order(:mean_abs_spearman, rev=true)])
CSV.write(joinpath(cots_dir, "exploratory_parameter_screen.csv"), screen)

coral_by_reef = CSV.read(joinpath(coral_dir, "evaluated_by_reef.csv"), DataFrame)
cots_by_reef = CSV.read(joinpath(cots_dir, "evaluated_by_reef.csv"), DataFrame)
peak_count(value) = isempty(strip(string(value))) ? 0 : length(split(string(value), ';'))
coral_second_peaks = count(>=(2), peak_count.(coral_by_reef.sim_peak_years))
cots_second_peaks = count(>=(2), peak_count.(cots_by_reef.sim_peak_years))

summary = Dict(
    "coral_run_id" => coral_run_id,
    "cots_run_id" => cots_run_id,
    "paired_candidates" => nrow(paired),
    "cots_candidates_with_lower_loss" => count(<(0.0), paired.delta_loss),
    "mean_loss_delta_cots_minus_coral" => mean(paired.delta_loss),
    "mean_peak_timing_delta_cots_minus_coral" =>
        mean(paired.delta_mean_peak_timing_penalty),
    "mean_amplitude_delta_cots_minus_coral" =>
        mean(paired.delta_mean_amplitude_penalty),
    "coral_candidate_reef_rows_with_two_or_more_peaks" => coral_second_peaks,
    "cots_candidate_reef_rows_with_two_or_more_peaks" => cots_second_peaks,
    "screening_note" =>
        "Exploratory rank correlations only (six candidates); use to prioritize, not infer importance."
)
open(joinpath(cots_dir, "connectivity_factorial_summary.toml"), "w") do io
    TOML.print(io, summary; sorted=true)
end

println(summary)
println("Top exploratory peak/amplitude associations:")
show(
    stdout,
    MIME("text/plain"),
    first(screen[in.(screen.outcome, Ref([
        "mean_peak_timing_penalty",
        "mean_amplitude_penalty"
    ])), :], 10)
)
println()
