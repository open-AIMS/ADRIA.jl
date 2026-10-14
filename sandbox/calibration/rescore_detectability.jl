"""Re-score stored 1991-run trajectories under the preregistered density-dependent
manta observation model (crash_mechanism_protocol.md, 2026-10-13):

    CPUE = q * A * A / (A + A_half),  q = 0.015 COTS/tow per adult/ha

Both observation treatments of the runner (raw, smoothed_3y) are scored.
A_half = 0 is the frozen linear model and is scored exactly as the runner does
(fixed_observation_scale = q on adults), so it must reproduce the stored losses.
Simulations are not rerun; only the observation mapping changes.

Usage: julia --project=sandbox rescore_detectability.jl [--q=<q>] <output_dir> <run_dir>:<arm,arm,...> ...

`--q` replaces the frozen q = 0.015 (e.g. the Pratchett 2005 paired fits, 2026-10-14); the
exact-replay check of the A_half = 0 control then does not apply and is recorded as such.

Every row also carries frozen score v2 (score_v2.jl, protocol 2026-10-14): `loss_v2` (Eyrie's
unobservable first wave dropped from the peak count; other reefs identical to `loss`) and
`score_v2_primary` (the row is the v2 treatment for its reef: raw for North Direction,
smoothed_3y otherwise). Summed v2 loss = sum of `loss_v2` over primary rows.
"""

using Pkg
const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
Pkg.activate(joinpath(ROOT, "sandbox"))

using CSV, DataFrames, Dates, TOML
include(joinpath(@__DIR__, "cots_cycle_metrics.jl"))
include(joinpath(@__DIR__, "calibration_run.jl"))
include(joinpath(@__DIR__, "score_v2.jl"))

const Q_FROZEN = 0.015
const Q = let flag = filter(a -> startswith(a, "--q="), ARGS)
    isempty(flag) ? Q_FROZEN : parse(Float64, split(only(flag), "=")[2])
end
0.0 < Q < 1.0 || error("q must lie in (0, 1)")
const A_HALF = (0.0, 1.0, 3.0, 10.0, 30.0)
const OBS = joinpath(ROOT, "sandbox", "data", "reef_cots.csv")
const TARGETS = Dict("Lizard Island Reef" => "Lizard Isles",
    "MacGillivray Reef" => "Macgillivray Reef",
    "North Direction Reef" => "North Direction Island",
    "Eyrie Reef" => "Eyrie Reef")

detected_cpue(adults, a_half) = a_half == 0.0 ? Q .* adults :
    Q .* adults .* adults ./ (adults .+ a_half)

function main(args)
    args = filter(a -> !startswith(a, "--q="), args)
    length(args) >= 2 || error("Usage: rescore_detectability.jl <output_dir> <run_dir>:<arms> ...")
    out = abspath(args[1])
    isfile(joinpath(out, "detectability_scores.csv")) && error("Refusing to overwrite $out")
    mkpath(out)
    raw_obs = CSV.read(OBS, DataFrame)
    scores = DataFrame(run=String[], treatment=String[], reef_name=String[],
        observation_treatment=String[], a_half=Float64[],
        loss=Float64[], stored_frozen_loss=Float64[], sim_peak_years=String[],
        sim_peak_heights_cpue=String[], matched_peaks=Int[], obs_peak_years=String[],
        peak_count_penalty=Float64[], loss_v2=Float64[], score_v2_primary=Bool[])
    series = DataFrame(run=String[], treatment=String[], reef_name=String[], a_half=Float64[],
        year=Int[], adults_ha=Float64[], detected_cpue=Float64[])
    inputs = Dict{String,String}()
    for spec in args[2:end]
        run_dir, arm_list = rsplit(spec, ":"; limit=2)
        run_dir = abspath(run_dir)
        arms = split(arm_list, ",")
        trajectories = CSV.read(joinpath(run_dir, "trajectories.csv"), DataFrame)
        stored = CSV.read(joinpath(run_dir, "scores.csv"), DataFrame)
        for name in ("trajectories.csv", "scores.csv", "metadata.toml")
            inputs[joinpath(basename(run_dir), name)] = file_sha256(joinpath(run_dir, name))
        end
        for arm in arms, (reef, observed) in TARGETS
            rows = sort(trajectories[(trajectories.treatment .== arm) .&
                (trajectories.reef_name .== reef), :], :year)
            nrow(rows) == 34 || error("Missing trajectory rows for $arm / $reef in $run_dir")
            years = Int.(rows.year)
            adults = Float64.(rows.adults_ha)
            obs = sort(raw_obs[(raw_obs.reef_name .== observed) .& (raw_obs.year .>= 1992) .&
                (raw_obs.year .<= 2024), :], :year)
            obs_years = Int.(obs.year)
            obs_values = Float64.(obs.cotsptow)
            smoothed = calendar_moving_average(obs_years, obs_values; window_years=3)
            for (treatment, values) in (("raw", obs_values), ("smoothed_3y", smoothed))
            frozen = stored[(stored.treatment .== arm) .& (stored.reef_name .== reef) .&
                (stored.observation_treatment .== treatment), :loss]
            length(frozen) == 1 || error("Missing stored score for $arm / $reef / $treatment")
            for a_half in A_HALF
                # Control: the runner's exact call (scale q on adults). Treatment: detected CPUE, scale 1.
                sim, scale = a_half == 0.0 ? (adults, Q) : (detected_cpue(adults, a_half), 1.0)
                all(isfinite, sim) || error("Nonfinite detected CPUE")
                score = cots_cycle_score(years[2:end], sim[2:end], obs_years, values;
                    fixed_observation_scale=scale)
                window = (years .>= score.scoring_start_year) .& (years .<= score.scoring_end_year)
                peaks = detect_cots_peaks(years[window], (sim .* scale)[window]; smooth_window_years=3)
                push!(scores, (basename(run_dir), arm, reef, treatment, a_half, score.total_loss, only(frozen),
                    join(score.sim_peak_years, ';'), join(round.([p.value for p in peaks]; digits=5), ';'),
                    score.n_matched_peaks, join(score.obs_peak_years, ';'), score.peak_count_penalty,
                    score_v2_loss(reef, score.total_loss, score.peak_count_penalty,
                        score.sim_peak_years, score.obs_peak_years),
                    score_v2_primary(reef, treatment)))
                treatment == "raw" || continue
                for (year, adult, value) in zip(years, adults, sim .* scale)
                    push!(series, (basename(run_dir), arm, reef, a_half, year, adult, value))
                end
            end
            end
        end
    end
    control = scores[scores.a_half .== 0.0, :]
    replay = maximum(abs.(control.loss .- control.stored_frozen_loss))
    Q != Q_FROZEN || replay == 0.0 || error("A_half = 0 control does not reproduce stored losses: $replay")
    CSV.write(joinpath(out, "detectability_scores.csv"), scores)
    CSV.write(joinpath(out, "detectability_series.csv"), series)
    metadata = Dict(
        "status" => "observation_sensitivity_not_promoted",
        "observation_model" => "CPUE = q * A * A / (A + A_half); A_half = 0 is the frozen linear model",
        "score_v2" => "loss_v2 / score_v2_primary columns (crash_mechanism_protocol.md, frozen score v2, 2026-10-14)",
        "q_cots_tow_per_adult_ha" => Q, "a_half_adults_ha" => collect(A_HALF),
        "control_max_abs_loss_delta" => Q == Q_FROZEN ? replay : "not applicable (q override $(Q))",
        "runs_and_arms" => collect(args[2:end]),
        "created" => string(now()),
        "adria_revision" => git_revision(ROOT),
        "cotsmod_revision" => git_revision(joinpath(ROOT, "..", "COTSMod.jl")),
        "inputs" => merge(inputs, Dict(
            "observations" => file_sha256(OBS),
            "score_definition" => file_sha256(joinpath(@__DIR__, "cots_cycle_metrics.jl")),
            "score_v2_definition" => file_sha256(joinpath(@__DIR__, "score_v2.jl")),
            "script" => file_sha256(@__FILE__))),
        "outputs" => Dict(name => file_sha256(joinpath(out, name)) for name in
            ("detectability_scores.csv", "detectability_series.csv")),
    )
    open(joinpath(out, "metadata.toml"), "w") do io
        TOML.print(io, metadata)
    end
    println("Rescored ", nrow(scores), " reef-arm-A_half rows; control replay delta ", replay)
end

main(ARGS)
