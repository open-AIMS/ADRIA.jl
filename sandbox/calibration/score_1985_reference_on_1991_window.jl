"""Re-score the archived 1985 Owen reference only on the common 1992-2024 window."""

using Pkg
const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
Pkg.activate(joinpath(ROOT, "sandbox"))
cd(ROOT)

using CSV, DataFrames, TOML
include(joinpath(@__DIR__, "cots_cycle_metrics.jl"))
include(joinpath(@__DIR__, "calibration_run.jl"))

const RUN = get(ENV, "OWEN_1991_RUN_DIR",
    joinpath(@__DIR__, "runs", "20261007T199105_owen_1991_initialization"))
const ARCHIVE = joinpath(@__DIR__, "runs", "20261007T103000_owen_seed_boundary_screen",
    "trajectories.csv")
const GATE = joinpath(@__DIR__, "runs", "20261001T113401_lizard_domain_gate",
    "scores.csv")
const OBS = joinpath(ROOT, "sandbox", "data", "reef_cots.csv")
const OUTPUT = joinpath(RUN, "historical_reference_common_window_scores.csv")
isfile(OUTPUT) && error("Refusing to overwrite $OUTPUT")
targets = [
    (model="Lizard Island Reef", observed="Lizard Isles"),
    (model="MacGillivray Reef", observed="Macgillivray Reef"),
    (model="North Direction Reef", observed="North Direction Island"),
    (model="Eyrie Reef", observed="Eyrie Reef"),
]
archive = CSV.read(ARCHIVE, DataFrame)
archive = archive[archive.treatment .== "reference_full", :]
gate = CSV.read(GATE, DataFrame)
gate = gate[(gate.seed .== 20260930) .& (gate.treatment .== "v1_legacy_closed") .&
    (gate.observation_treatment .== "raw"), :]
scales = Dict(String(row.reef_name) => Float64(row.observation_scale) for row in eachrow(gate))
observations = CSV.read(OBS, DataFrame)
scores = DataFrame(reef_name=String[], observation_treatment=String[], loss=Float64[],
    sim_peak_years=String[], sim_peak_heights_cpue=String[], matched_peaks=Int[])
for target in targets
    trajectory = sort(archive[(archive.reef_name .== target.model) .&
        (archive.year .>= 1992) .& (archive.year .<= 2024), :], :year)
    nrow(trajectory) == 33 || error("Incomplete archived trajectory: $(target.model)")
    obs = sort(observations[(observations.reef_name .== target.observed) .&
        (observations.year .>= 1992) .& (observations.year .<= 2024), :], :year)
    obs_years = Int.(obs.year)
    obs_raw = Float64.(obs.cotsptow)
    obs_smoothed = calendar_moving_average(obs_years, obs_raw; window_years=3)
    years = Int.(trajectory.year)
    adults = Float64.(trajectory.adults_ha)
    scale = scales[target.model]
    for (name, values) in (("raw", obs_raw), ("smoothed_3y", obs_smoothed))
        score = cots_cycle_score(years, adults, obs_years, values;
            fixed_observation_scale=scale)
        window = (years .>= score.scoring_start_year) .&
            (years .<= score.scoring_end_year)
        peaks = detect_cots_peaks(years[window], (adults .* scale)[window];
            smooth_window_years=3)
        push!(scores, (target.model, name, score.total_loss,
            join(score.sim_peak_years, ';'),
            join(round.([p.value for p in peaks]; digits=5), ';'),
            score.n_matched_peaks))
    end
end
CSV.write(OUTPUT, scores)
metadata = Dict(
    "status" => "common_window_reference_for_1991_screen",
    "score_years" => [1992, 2024],
    "archive" => file_sha256(ARCHIVE),
    "observations" => file_sha256(OBS),
    "frozen_scales" => file_sha256(GATE),
    "script" => file_sha256(@__FILE__),
    "output" => file_sha256(OUTPUT),
)
open(joinpath(RUN, "historical_reference_common_window_metadata.toml"), "w") do io
    TOML.print(io, metadata)
end
println(scores)
