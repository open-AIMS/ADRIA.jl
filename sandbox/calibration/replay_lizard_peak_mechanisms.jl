"""Bounded V1 replay: distinguish old parameter, Allee, network and recruitment effects."""

using Pkg
const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
Pkg.activate(joinpath(ROOT, "sandbox"))
cd(ROOT)

using ADRIA
using CSV, DataFrames, Dates, Random, Statistics, TOML

include(joinpath(@__DIR__, "cots_cycle_metrics.jl"))
include(joinpath(@__DIR__, "reef_observation_mapping.jl"))
include(joinpath(@__DIR__, "calibration_run.jl"))

const SEED = 20260930
const RUN_ID = get(ENV, "LIZARD_REPLAY_RUN_ID",
    Dates.format(now(), "yyyymmddTHHMMSS") * "_lizard_peak_mechanism_replay")
occursin(r"^[A-Za-z0-9_.-]+$", RUN_ID) || error("Invalid run ID")
const OUT = joinpath(ROOT, "sandbox", "calibration", "runs", RUN_ID)
ispath(OUT) && error("Refusing to overwrite $OUT")
mkpath(OUT)
const DOMAIN_PATH = joinpath(ROOT, "sandbox", "data", "Lizard_Historical_v0.1")
const OLD_PATH = joinpath(ROOT, "sandbox", "data", "bbo_best_summary.csv")
const CURRENT_PATH = joinpath(ROOT, "sandbox", "calibration", "runs",
    "expanded_cotsconn_pilot_seed20260930", "best_summary.csv")
const GATE_DIR = joinpath(ROOT, "sandbox", "calibration", "runs",
    "20261001T113401_lizard_domain_gate")
const GATE_SCORES_PATH = joinpath(GATE_DIR, "scores.csv")
const OBS_PATH = joinpath(ROOT, "sandbox", "data", "reef_cots.csv")
const TARGETS = [
    (model="Lizard Island Reef", observed="Lizard Isles"),
    (model="MacGillivray Reef", observed="Macgillivray Reef"),
    (model="North Direction Reef", observed="North Direction Island"),
    (model="Eyrie Reef", observed="Eyrie Reef"),
]
const FLOW_NAMES = ["local_fecundity", "background_immigration", "internal_immigration",
    "external_immigration", "maturation", "retained_juveniles", "settlement_gate",
    "local_retention", "pelagic_survivors", "settled_recruits",
    "external_potential", "external_pelagic"]

old = sort(CSV.read(OLD_PATH, DataFrame), :loss)[1, :]
current = sort(CSV.read(CURRENT_PATH, DataFrame), :loss)[1, :]
gate_scores = CSV.read(GATE_SCORES_PATH, DataFrame)
scale_rows = gate_scores[(gate_scores.seed .== SEED) .&
    (gate_scores.treatment .== "v1_legacy_closed") .&
    (gate_scores.observation_treatment .== "raw"), :]
nrow(scale_rows) == length(TARGETS) || error("Frozen observation scales not found")
fixed_scales = Dict(String(row.reef_name) => Float64(row.observation_scale)
    for row in eachrow(scale_rows))

ENV["ADRIA_COTS_CONNECTIVITY_DATASET"] = "legacy"
ENV["COTS_CONNECTIVITY_TEMPORAL_MODE"] = "mean"
ENV["COTS_CONNECTIVITY_SEED"] = string(SEED)
ENV["COTS_CALENDAR_START_YEAR"] = "1985"
ENV["COTS_EXTERNAL_PULSE"] = "false"
ENV["COTS_EXTERNAL_SOURCE_DENSITY"] = "0.0"
ENV["COTS_EXTERNAL_OUTBREAK_PRODUCTION"] = "0.0"
ENV["COTS_JUVENILE_STORAGE"] = "false"
ENV["COTS_HABITAT_MEDIATION"] = "false"
ENV["COTS_SEED_FIRST_N"] = "10"
ENV["ADRIA_DEBUG_SEED_FIRST_N"] = "10"
ENV["COTS_SETTLEMENT_PROBABILITY"] = "1.0"

domains = Dict{String,Any}()
for mode in ("coral", "cots")
    ENV["ADRIA_COTS_CONNECTIVITY_MODE"] = mode
    domain = ADRIA.load_domain(ADRIA.LizardDomain, DOMAIN_PATH, "historical")
    domains[mode] = (domain=domain,
        mapping=load_reef_area_weights(domain, DOMAIN_PATH, [t.model for t in TARGETS]))
end

function scenario_from_candidate(domain, candidate, threshold)
    Random.seed!(SEED)
    scenario = ADRIA.sample(domain, 2)[1:1, :]
    defaults = ADRIA.param_table(domain)
    for name in names(scenario)
        scenario[!, name] .= defaults[1, name]
        if name != "allee_threshold" && hasproperty(candidate, Symbol(name))
            scenario[!, name] .= candidate[Symbol(name)]
        end
    end
    scenario[!, :allee_threshold] .= threshold
    return scenario
end

const TREATMENTS = [
    (name="old_coral_allee1_counterfactual", candidate="old", mode="coral",
     allee=1.0, separated=false, survival=false),
    (name="old_coral_allee3", candidate="old", mode="coral",
     allee=3.0, separated=false, survival=false),
    (name="old_cots_allee3", candidate="old", mode="cots",
     allee=3.0, separated=false, survival=false),
    (name="current_cots_legacy_allee3", candidate="current", mode="cots",
     allee=3.0, separated=false, survival=false),
    (name="current_cots_separated_nosurvival_allee3", candidate="current", mode="cots",
     allee=3.0, separated=true, survival=false),
    (name="current_cots_separated_survival_allee3", candidate="current", mode="cots",
     allee=3.0, separated=true, survival=true),
]

raw_obs = CSV.read(OBS_PATH, DataFrame)
observed = Dict{String,NamedTuple}()
for target in TARGETS
    rows = sort(raw_obs[(raw_obs.reef_name .== target.observed) .&
        (raw_obs.year .>= 1985) .& (raw_obs.year .<= 2024), :], :year)
    years = Int.(rows.year)
    observed[target.model] = (years=years, raw=Float64.(rows.cotsptow),
        smooth=calendar_moving_average(years, Float64.(rows.cotsptow); window_years=3))
end

scores = DataFrame(treatment=String[], reef_name=String[], observation_treatment=String[],
    loss=Float64[], sim_peak_count=Int[], sim_peak_years=String[],
    matched_peaks=Int[], amplitude_penalty=Float64[], peak_timing_penalty=Float64[],
    observation_scale=Float64[])
trajectories = DataFrame(treatment=String[], reef_name=String[], year=Int[],
    recruits_ha=Float64[], juveniles_ha=Float64[], adults_ha=Float64[],
    simulated_cpue=Float64[], coral_cover=Float64[], body_condition=Float64[])
flows = DataFrame(treatment=String[], reef_name=String[], year=Int[])
for name in FLOW_NAMES
    flows[!, name] = Float64[]
end
failures = DataFrame(treatment=String[], error=String[])

for treatment in TREATMENTS
    candidate = treatment.candidate == "old" ? old : current
    entry = domains[treatment.mode]
    ENV["ADRIA_COTS_CONNECTIVITY_MODE"] = treatment.mode
    ENV["COTS_ALLEE_THRESHOLD"] = string(treatment.allee)
    ENV["COTS_SEPARATED_RECRUITMENT"] = string(treatment.separated)
    ENV["COTS_APPLY_LARVAL_SURVIVAL"] = string(treatment.survival)
    ENV["COTS_INITIAL_MULTIPLIER"] = string(candidate.seed_mult)
    # Fecundity is an explicit pre-pelagic parameter in separated mode; keep
    # the same frozen comparison conversion used in the V1/V2 gate.
    reference_survival = mean(entry.domain.cots_forcing.source_survival)
    ENV["COTS_LARVAL_FECUNDITY"] = string(Float64(current.a_ricker) / reference_survival)
    scenario = scenario_from_candidate(entry.domain, candidate, treatment.allee)
    Random.seed!(SEED)
    try
        result = ADRIA.run_scenario(entry.domain, scenario[1, :])
        years = Int.(result.cots_calendar_year)
        years == collect(1985:2024) || error("Unexpected calendar years")
        cover = dropdims(sum(result.raw; dims=(2, 3)); dims=(2, 3))
        for target in TARGETS
            reef = target.model
            mapping = entry.mapping[reef]
            stages = [aggregate_reef_density(Matrix(result.cots_log[:, channel, :]), mapping)
                for channel in 1:3]
            reef_cover = aggregate_reef_density(cover, mapping)
            condition = aggregate_reef_density(Matrix(result.cots_condition_log), mapping)
            reef_flows = [aggregate_reef_density(Matrix(result.cots_flow_log[:, channel, :]),
                mapping) for channel in 1:length(FLOW_NAMES)]
            scale = fixed_scales[reef]
            for (index, year) in enumerate(years)
                push!(trajectories, (treatment.name, reef, year,
                    stages[1][index], stages[2][index], stages[3][index],
                    stages[3][index] * scale, reef_cover[index], condition[index]))
                push!(flows, (treatment.name, reef, year,
                    (values[index] for values in reef_flows)...))
            end
            obs = observed[reef]
            for (obs_name, obs_values) in (("raw", obs.raw), ("smoothed_3y", obs.smooth))
                score = cots_cycle_score(years, stages[3], obs.years, obs_values;
                    fixed_observation_scale=scale)
                push!(scores, (treatment.name, reef, obs_name, score.total_loss,
                    score.n_sim_peaks, join(score.sim_peak_years, ';'),
                    score.n_matched_peaks, score.amplitude_penalty,
                    score.peak_timing_penalty, scale))
            end
        end
        println("Completed ", treatment.name)
    catch err
        push!(failures, (treatment.name, sprint(showerror, err, catch_backtrace())))
        println(stderr, "FAILED ", treatment.name, ": ", err)
    end
    flush(stdout)
    CSV.write(joinpath(OUT, "scores.csv"), scores)
    CSV.write(joinpath(OUT, "trajectories.csv"), trajectories)
    CSV.write(joinpath(OUT, "flows.csv"), flows)
    CSV.write(joinpath(OUT, "failures.csv"), failures)
end

inputs = Dict(
    "old_candidate" => file_sha256(OLD_PATH),
    "current_candidate" => file_sha256(CURRENT_PATH),
    "frozen_gate_scores" => file_sha256(GATE_SCORES_PATH),
    "observations" => file_sha256(OBS_PATH),
    "v1_site_map" => file_sha256(joinpath(DOMAIN_PATH, "site_to_reef.csv")),
    "v1_forcing" => file_sha256(joinpath(DOMAIN_PATH, "cots_connectivity", "annual_forcing", "provenance.toml")),
    "script" => file_sha256(@__FILE__),
    "score_script" => file_sha256(joinpath(@__DIR__, "cots_cycle_metrics.jl")),
    "mapping_script" => file_sha256(joinpath(@__DIR__, "reef_observation_mapping.jl")),
    "scenario_code" => file_sha256(joinpath(ROOT, "ADRIA", "src", "scenario.jl")),
    "cots_core" => file_sha256(joinpath(ROOT, "..", "COTSMod.jl", "src", "COTSMod.jl")),
    "manifest" => file_sha256(joinpath(ROOT, "sandbox", "Manifest.toml")),
)
metadata = Dict(
    "status" => "exploratory_not_promoted",
    "seed" => SEED,
    "temporal_mode" => "mean_static_control",
    "legacy_allee1_label" => "counterfactual only; evidence-based threshold remains 3 COTS/ha",
    "observation_mapping" => "within-reef polygon-area-weighted COTS ha^-1 and coral fraction",
    "observation_scale" => "frozen per-reef raw-CPUE scales from V1 legacy closed reference seed in 20261001T113401 gate",
    "external_boundary" => "off in all treatments; V2 Owen boundary evaluated separately in paired gate",
    "old_replay_limit" => "archived simulation environment is not fully recorded, so this is a parameter-path replay, not an exact reproduction",
    "old_default_parameters" => "unset values from current V1 param_table; this may differ from September archived defaults",
    "julia_version" => string(VERSION),
    "adria_revision" => git_revision(ROOT),
    "cotsmod_revision" => git_revision(joinpath(ROOT, "..", "COTSMod.jl")),
    "git_status" => read(`git status --porcelain`, String),
    "treatments" => [Dict(String(k) => v for (k, v) in pairs(t)) for t in TREATMENTS],
    "inputs" => inputs,
    "outputs" => Dict(file => file_sha256(joinpath(OUT, file)) for file in
        ("scores.csv", "trajectories.csv", "flows.csv", "failures.csv")),
)
open(joinpath(OUT, "metadata.toml"), "w") do io
    TOML.print(io, metadata; sorted=true)
end
println("Saved ", OUT)
