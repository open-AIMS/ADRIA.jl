"""Held-out-years re-score for promotion gate 4a (crash_mechanism_protocol.md, "Preregistered
promotion checks", 2026-10-19). Uses the frozen v2/v3 scoring functions on truncated windows,
exactly as score v3 truncates its manta part:

- fit window:   v2 loss (primary treatment per reef) on manta 1992-2012, model and observations
                truncated before smoothing;
- check window: v2 loss on manta 2013-2021 (same truncation) + score v3's diver term
                (1.5 x mean |log((A + 0.5) / (T + 0.5))|, 2020-26 diver/cull targets).

Linear observation model (A_half = 0, q 0.015), as the primary v3 rows. The frozen scores are not
changed; this script only adds windowed losses.

Usage: julia --project=sandbox rescore_holdout.jl <output_dir> <run_dir>:<arm,arm,...> ...
"""

using Pkg
const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
Pkg.activate(joinpath(ROOT, "sandbox"))

using CSV, DataFrames, Dates, TOML
include(joinpath(@__DIR__, "cots_cycle_metrics.jl"))
include(joinpath(@__DIR__, "calibration_run.jl"))
include(joinpath(@__DIR__, "score_v2.jl"))
include(joinpath(@__DIR__, "score_v3.jl"))

const Q = 0.015
const FIT = (1992, 2012)
const CHECK = (2013, SCORE_V3_LAST_MANTA_YEAR)
const OBS = joinpath(ROOT, "sandbox", "data", "reef_cots.csv")
const DIVER_TARGETS = joinpath(ROOT, "sandbox", "data", "derived", "density_targets_2020s", "targets_reef_year.csv")
const TARGETS = Dict("Lizard Island Reef" => "Lizard Isles", "MacGillivray Reef" => "Macgillivray Reef",
    "North Direction Reef" => "North Direction Island", "Eyrie Reef" => "Eyrie Reef")

function windowed_v2(reef, treatment, years, adults, obs_years, obs_values, window)
    keep_obs = (obs_years .>= window[1]) .& (obs_years .<= window[2])
    oy, ov = obs_years[keep_obs], obs_values[keep_obs]
    ov = treatment == "raw" ? ov : calendar_moving_average(oy, ov; window_years=3)
    keep_sim = (years .>= window[1]) .& (years .<= window[2])
    s = cots_cycle_score(years[keep_sim], adults[keep_sim], oy, ov; fixed_observation_scale=Q)
    return score_v2_loss(reef, s.total_loss, s.peak_count_penalty, s.sim_peak_years, s.obs_peak_years),
        s.n_matched_peaks, length(s.obs_peak_years)
end

function main(args)
    length(args) >= 2 || error("Usage: rescore_holdout.jl <output_dir> <run_dir>:<arms> ...")
    out = abspath(args[1])
    isfile(joinpath(out, "holdout_scores.csv")) && error("Refusing to overwrite $out")
    mkpath(out)
    raw_obs = CSV.read(OBS, DataFrame)
    table = CSV.read(DIVER_TARGETS, DataFrame)
    table = table[table.report_year .>= SCORE_V3_FIRST_DIVER_YEAR, :]
    diver = Dict(String(r) => (Int.(table[table.target_reef .== r, :report_year]),
        Float64.(table[table.target_reef .== r, :target_density_ha])) for r in unique(table.target_reef))
    rows = DataFrame(run=String[], treatment=String[], reef_name=String[], v2_treatment=String[],
        loss_fit=Float64[], matched_fit=Int[], observed_fit=Int[], loss_check_manta=Float64[],
        matched_check=Int[], observed_check=Int[], diver_error=Float64[], loss_check=Float64[],
        loss_1992_2021_selfcheck=Float64[])
    inputs = Dict{String,String}()
    for spec in args[2:end]
        run_dir, arm_list = rsplit(spec, ":"; limit=2)
        run_dir = abspath(run_dir)
        traj = CSV.read(joinpath(run_dir, "trajectories.csv"), DataFrame)
        inputs[joinpath(basename(run_dir), "trajectories.csv")] = file_sha256(joinpath(run_dir, "trajectories.csv"))
        for arm in split(arm_list, ","), (reef, observed) in TARGETS
            r = sort(traj[(traj.treatment .== arm) .& (traj.reef_name .== reef), :], :year)
            nrow(r) in (34, 36) || error("Missing trajectory rows for $arm / $reef")
            years, adults = Int.(r.year)[2:end], Float64.(r.adults_ha)[2:end]   # as the frozen score (from 1992)
            obs = sort(raw_obs[(raw_obs.reef_name .== observed) .& (raw_obs.year .>= 1992) .& (raw_obs.year .<= 2024), :], :year)
            oy, ov = Int.(obs.year), Float64.(obs.cotsptow)
            treatment = SCORE_V2_TREATMENT[reef]
            fit, mf, of = windowed_v2(reef, treatment, years, adults, oy, ov, FIT)
            chk, mc, oc = windowed_v2(reef, treatment, years, adults, oy, ov, CHECK)
            full = windowed_v2(reef, treatment, years, adults, oy, ov, (1992, SCORE_V3_LAST_MANTA_YEAR))[1]
            de = haskey(diver, reef) ? diver_log_error(Int.(r.year), Float64.(r.adults_ha), diver[reef]...)[1] : NaN
            push!(rows, (basename(run_dir), String(arm), reef, treatment, fit, mf, of, chk, mc, oc, de, score_v3_loss(chk, de), full))
        end
    end
    CSV.write(joinpath(out, "holdout_scores.csv"), rows)
    meta = Dict("status" => "promotion_gate_4a_holdout_rescore", "fit_window" => collect(FIT), "check_window_manta" => collect(CHECK),
        "check_diver" => "score v3 diver term (2020-26)", "q" => Q, "observation_model" => "linear (A_half 0)",
        "runs_and_arms" => collect(args[2:end]), "created" => string(now()),
        "adria_revision" => git_revision(ROOT), "inputs" => inputs, "observations" => file_sha256(OBS),
        "diver_targets" => file_sha256(DIVER_TARGETS), "script" => file_sha256(@__FILE__),
        "outputs" => Dict("holdout_scores.csv" => file_sha256(joinpath(out, "holdout_scores.csv"))))
    open(joinpath(out, "metadata.toml"), "w") do io
        TOML.print(io, meta)
    end
    println("Held-out re-score: ", nrow(rows), " reef-arm rows")
end

main(ARGS)
