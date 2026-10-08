using Pkg
const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
Pkg.activate(joinpath(REPO_ROOT, "sandbox"))

using CSV
using DataFrames

run_id = get(ENV, "BBO_RUN_ID", "expanded_cotsconn_pilot_seed20260930")
run_dir = joinpath(REPO_ROOT, "sandbox", "calibration", "runs", run_id)
trajectories = CSV.read(joinpath(run_dir, "best_calibrated_trajectories.csv"), DataFrame)
best = first(CSV.read(joinpath(run_dir, "best_summary.csv"), DataFrame))
split_year = parse(Int, get(ENV, "COTS_TEMPORAL_SPLIT_YEAR", "2005"))
allee_threshold = parse(Float64, get(ENV, "COTS_ALLEE_THRESHOLD", "3.0"))

diagnostics = DataFrame(
    reef_name=String[],
    first_peak_year=Int[],
    first_peak_adults=Float64[],
    post_split_max_year=Int[],
    post_split_max_adults=Float64[],
    final_year=Int[],
    final_recruits=Float64[],
    final_juveniles=Float64[],
    final_adults=Float64[],
    final_body_condition=Float64[],
    final_coral_cover=Float64[],
    allee_threshold=Float64[],
    final_adult_to_allee_threshold=Float64[],
    final_allee_multiplier=Float64[]
)

for reef in unique(trajectories.reef_name)
    reef_rows = sort(trajectories[trajectories.reef_name .== reef, :], :year)
    pre = reef_rows[reef_rows.year .<= split_year, :]
    post = reef_rows[reef_rows.year .> split_year, :]
    isempty(pre) && continue
    isempty(post) && continue
    pre_peak = pre[argmax(pre.sim_cots_adult), :]
    post_peak = post[argmax(post.sim_cots_adult), :]
    final = post[end, :]
    final_adults = Float64(final.sim_cots_adult)
    allee_multiplier = final_adults^2 / (allee_threshold^2 + final_adults^2)
    push!(diagnostics, (
        string(reef),
        Int(pre_peak.year),
        Float64(pre_peak.sim_cots_adult),
        Int(post_peak.year),
        Float64(post_peak.sim_cots_adult),
        Int(final.year),
        Float64(final.sim_cots_recruits),
        Float64(final.sim_cots_juveniles),
        final_adults,
        Float64(final.sim_cots_condition),
        Float64(final.sim_coral_cover),
        allee_threshold,
        final_adults / allee_threshold,
        allee_multiplier
    ))
end

sort!(diagnostics, :reef_name)
output_path = joinpath(run_dir, "recurrence_diagnostics.csv")
CSV.write(output_path, diagnostics)
show(stdout, MIME("text/plain"), diagnostics)
println()
println("Wrote $output_path")
