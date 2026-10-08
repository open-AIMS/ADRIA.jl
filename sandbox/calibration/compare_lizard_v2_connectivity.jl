using Pkg
const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
Pkg.activate(joinpath(REPO_ROOT, "sandbox"))
cd(REPO_ROOT)

using ADRIA
using CSV, DataFrames, Dates, LinearAlgebra, Random, SHA, Statistics, TOML

include(joinpath(@__DIR__, "cots_cycle_metrics.jl"))
include(joinpath(@__DIR__, "calibration_run.jl"))

const SEED = parse(Int, get(ENV, "COTS_V2_PILOT_SEED", "20260930"))
const RUN_ID = get(
    ENV, "COTS_V2_PILOT_RUN_ID",
    Dates.format(now(), "yyyymmddTHHMMSS") * "_lizard_v2_connectivity_seed" * string(SEED)
)
occursin(r"^[A-Za-z0-9_.-]+$", RUN_ID) || error("Unsupported COTS_V2_PILOT_RUN_ID")
const OUT_DIR = joinpath(REPO_ROOT, "sandbox", "calibration", "runs", RUN_ID)
isdir(OUT_DIR) && !isempty(readdir(OUT_DIR)) && error("Run directory already contains outputs: $OUT_DIR")
mkpath(OUT_DIR)

const DOMAIN_PATH = joinpath(REPO_ROOT, "sandbox", "data", "Lizard_Historical_v0.2")
const CANDIDATE_PATH = joinpath(
    REPO_ROOT, "sandbox", "calibration", "runs",
    "expanded_cotsconn_pilot_seed20260930", "best_summary.csv"
)
isfile(CANDIDATE_PATH) || error("Missing frozen candidate: $CANDIDATE_PATH")
best = sort(CSV.read(CANDIDATE_PATH, DataFrame), :loss)[1, :]

ENV["ADRIA_COTS_CONNECTIVITY_MODE"] = "cots"
ENV["COTS_CONNECTIVITY_TEMPORAL_MODE"] = "cycle"
ENV["COTS_CONNECTIVITY_SEED"] = string(SEED)
ENV["COTS_CALENDAR_START_YEAR"] = "1985"
ENV["COTS_ALLEE_THRESHOLD"] = "3.0"
ENV["COTS_EXTERNAL_PULSE"] = "false"
ENV["COTS_EXTERNAL_SOURCE_DENSITY"] = "0.0"
ENV["COTS_JUVENILE_STORAGE"] = "false"
ENV["COTS_HABITAT_MEDIATION"] = "false"
ENV["COTS_SEPARATED_RECRUITMENT"] = "true"
ENV["COTS_APPLY_LARVAL_SURVIVAL"] = "true"
ENV["COTS_SETTLEMENT_PROBABILITY"] = "1.0"
ENV["COTS_INITIAL_MULTIPLIER"] = string(best.seed_mult)
ENV["COTS_SEED_FIRST_N"] = "10"
ENV["ADRIA_DEBUG_SEED_FIRST_N"] = "10"

site_map = CSV.read(joinpath(DOMAIN_PATH, "site_to_reef.csv"), DataFrame)
site_map.reef_name_clean = first.(split.(String.(site_map.reef_name), " ("))
@assert nrow(site_map) == 2895
@assert site_map.site_id == ["FNG_V2_$(lpad(i, 4, '0'))" for i in 1:2895]
all(site_map.area_m2 .> 0.0) || error("Site areas must be positive")
target_map = Dict(
    "Lizard Island Reef" => "Lizard Isles",
    "MacGillivray Reef" => "Macgillivray Reef",
    "North Direction Reef" => "North Direction Island",
    "Eyrie Reef" => "Eyrie Reef"
)
observations = CSV.read(joinpath(REPO_ROOT, "sandbox", "data", "reef_cots.csv"), DataFrame)

ENV["ADRIA_COTS_CONNECTIVITY_DATASET"] = "legacy"
control_domain = ADRIA.load_domain(ADRIA.LizardDomain, DOMAIN_PATH, "historical")
reference_survival = mean(control_domain.cots_forcing.source_survival)
reference_survival > 0.0 || error("Legacy mean larval survival must be positive")
larval_fecundity = Float64(best.a_ricker) / reference_survival
ENV["COTS_LARVAL_FECUNDITY"] = string(larval_fecundity)
allee_threshold = 3.0
b_ricker = Float64(best.b_ricker)
production(adults) = larval_fecundity * adults * exp(-b_ricker * adults) *
    adults^2 / (allee_threshold^2 + adults^2)
adult_grid = range(0.0, 100.0; length=10001)
boundary_production = maximum(production.(adult_grid))

scenario = ADRIA.sample(control_domain, 2)[1:1, :]
defaults = ADRIA.param_table(control_domain)
for name in names(scenario)
    scenario[!, name] .= defaults[1, name]
    if name != "allee_threshold" && hasproperty(best, Symbol(name))
        scenario[!, name] .= best[Symbol(name)]
    end
end
scenario[!, :allee_threshold] .= allee_threshold

treatments = [
    (name="legacy_closed", dataset="legacy", boundary=0.0),
    (name="legacy_evidence_boundary", dataset="legacy", boundary=boundary_production),
    (name="owen_closed", dataset="owen_global_2018_2023", boundary=0.0),
    (name="owen_evidence_boundary", dataset="owen_global_2018_2023", boundary=boundary_production)
]

scores = DataFrame(
    treatment=String[], dataset=String[], reef_name=String[], loss=Float64[],
    sim_peak_count=Int[], sim_peak_years=String[], sim_peak_heights_cpue=String[],
    sim_peak_prominences_cpue=String[], matched_peaks=Int[],
    peak_timing_penalty=Float64[], amplitude_penalty=Float64[],
    observation_scale=Float64[], mean_external_pelagic=Float64[]
)
trajectories = DataFrame(
    treatment=String[], reef_name=String[], year=Int[], adults_ha=Float64[],
    coral_cover=Float64[], external_pelagic_ha=Float64[], forcing_year=Int[]
)
control_scale = Dict{String,Float64}()

for treatment in treatments
    ENV["ADRIA_COTS_CONNECTIVITY_DATASET"] = treatment.dataset
    ENV["COTS_EXTERNAL_OUTBREAK_PRODUCTION"] = string(treatment.boundary)
    domain = treatment.dataset == "legacy" ? control_domain :
        ADRIA.load_domain(ADRIA.LizardDomain, DOMAIN_PATH, "historical")
    Random.seed!(SEED)
    result = ADRIA.run_scenario(domain, scenario[1, :])
    years = Int.(result.cots_calendar_year)
    years == collect(1985:2024) || error("Unexpected model calendar years")
    coral_cover = dropdims(sum(result.raw; dims=(2, 3)); dims=(2, 3))
    for reef_name in keys(target_map)
        indices = findall(site_map.reef_name_clean .== reef_name)
        isempty(indices) && error("Missing target reef $reef_name in V2 sites")
        weights = Float64.(site_map.area_m2[indices])
        weights ./= sum(weights)
        adults = [dot(weights, result.cots_log[t, 3, indices]) for t in eachindex(years)]
        covers = [dot(weights, coral_cover[t, indices]) for t in eachindex(years)]
        external = [dot(weights, result.cots_flow_log[t, 12, indices]) for t in eachindex(years)]
        for t in eachindex(years)
            push!(trajectories, (
                treatment.name, reef_name, years[t], adults[t], covers[t],
                external[t], result.cots_forcing_year[t]
            ))
        end
        observed = observations[observations.reef_name .== target_map[reef_name], :]
        if treatment.name == "legacy_closed"
            provisional = cots_cycle_score(
                years, adults, Int.(observed.year), Float64.(observed.cotsptow)
            )
            control_scale[reef_name] = provisional.observation_scale
        end
        score = cots_cycle_score(
            years, adults, Int.(observed.year), Float64.(observed.cotsptow);
            fixed_observation_scale=control_scale[reef_name]
        )
        scoring_window = (years .>= score.scoring_start_year) .&
            (years .<= score.scoring_end_year)
        peaks = detect_cots_peaks(
            years[scoring_window],
            (adults .* control_scale[reef_name])[scoring_window];
            smooth_window_years=3
        )
        push!(scores, (
            treatment.name, treatment.dataset, reef_name, score.total_loss,
            score.n_sim_peaks, join(score.sim_peak_years, ';'),
            join(round.([peak.value for peak in peaks]; digits=4), ';'),
            join(round.([peak.prominence for peak in peaks]; digits=4), ';'),
            score.n_matched_peaks, score.peak_timing_penalty,
            score.amplitude_penalty, control_scale[reef_name], mean(external)
        ))
    end
    println("Completed ", treatment.name)
end

CSV.write(joinpath(OUT_DIR, "scores.csv"), scores)
CSV.write(joinpath(OUT_DIR, "trajectories.csv"), trajectories)
metadata = Dict(
    "design" => "same-seed 2 x 2 Lizard v0.2 legacy/Owen by closed/evidence boundary pilot",
    "seed" => SEED,
    "domain" => "Lizard_Historical_v0.2",
    "area_weighting" => "site polygon area in m^2, normalized within each target reef",
    "observation_scale" => "fitted once on legacy_closed, then frozen for all treatments",
    "allee_threshold_cots_ha" => allee_threshold,
    "prepelagic_fecundity" => larval_fecundity,
    "boundary_production_recruits_ha_year" => boundary_production,
    "source_survival_reference" => reference_survival,
    "candidate_path" => relpath(CANDIDATE_PATH, REPO_ROOT),
    "candidate_sha256" => file_sha256(CANDIDATE_PATH),
    "domain_provenance_sha256" => file_sha256(joinpath(DOMAIN_PATH, "provenance.json")),
    "legacy_forcing_sha256" => file_sha256(joinpath(DOMAIN_PATH, "cots_connectivity", "annual_forcing", "provenance.toml")),
    "owen_forcing_sha256" => file_sha256(joinpath(DOMAIN_PATH, "cots_connectivity", "datasets", "owen_global_2018_2023", "provenance.json")),
    "adria_revision" => git_revision(REPO_ROOT),
    "cotsmod_revision" => git_revision(normpath(joinpath(REPO_ROOT, "..", "COTSMod.jl"))),
    "results" => Dict(
        "scores_sha256" => file_sha256(joinpath(OUT_DIR, "scores.csv")),
        "trajectories_sha256" => file_sha256(joinpath(OUT_DIR, "trajectories.csv"))
    )
)
open(joinpath(OUT_DIR, "metadata.toml"), "w") do io
    TOML.print(io, metadata; sorted=true)
end
println("Wrote paired V2 connectivity pilot to ", OUT_DIR)
