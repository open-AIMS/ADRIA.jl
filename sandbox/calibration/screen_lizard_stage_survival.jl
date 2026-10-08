"""One-seed, fixed-Allee sensitivity of juvenile/adult survival with Owen boundary."""

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
const RUN_ID = get(ENV, "LIZARD_SURVIVAL_RUN_ID",
    Dates.format(now(), "yyyymmddTHHMMSS") * "_lizard_stage_survival_screen")
occursin(r"^[A-Za-z0-9_.-]+$", RUN_ID) || error("Invalid run ID")
const OUT = joinpath(ROOT, "sandbox", "calibration", "runs", RUN_ID)
ispath(OUT) && error("Refusing to overwrite $OUT")
mkpath(OUT)
const DOMAIN_PATH = joinpath(ROOT, "sandbox", "data", "Lizard_Historical_v0.2")
const CANDIDATE_PATH = joinpath(ROOT, "sandbox", "calibration", "runs",
    "expanded_cotsconn_pilot_seed20260930", "best_summary.csv")
const GATE_SCORES_PATH = joinpath(ROOT, "sandbox", "calibration", "runs",
    "20261001T113401_lizard_domain_gate", "scores.csv")
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

candidate = sort(CSV.read(CANDIDATE_PATH, DataFrame), :loss)[1, :]
gate_scores = CSV.read(GATE_SCORES_PATH, DataFrame)
scale_rows = gate_scores[(gate_scores.seed .== SEED) .&
    (gate_scores.treatment .== "v1_legacy_closed") .&
    (gate_scores.observation_treatment .== "raw"), :]
nrow(scale_rows) == length(TARGETS) || error("Frozen CPUE scales missing")
fixed_scales = Dict(String(row.reef_name) => Float64(row.observation_scale)
    for row in eachrow(scale_rows))

ENV["ADRIA_COTS_CONNECTIVITY_DATASET"] = "owen_global_2018_2023"
ENV["ADRIA_COTS_CONNECTIVITY_MODE"] = "cots"
ENV["COTS_CONNECTIVITY_TEMPORAL_MODE"] = "sample"
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
ENV["COTS_INITIAL_MULTIPLIER"] = string(candidate.seed_mult)
ENV["COTS_SEED_FIRST_N"] = "10"
ENV["ADRIA_DEBUG_SEED_FIRST_N"] = "10"

domain = ADRIA.load_domain(ADRIA.LizardDomain, DOMAIN_PATH, "historical")
mapping = load_reef_area_weights(domain, DOMAIN_PATH, [t.model for t in TARGETS])
survival_ref = mean(domain.cots_forcing.source_survival)
fecundity = Float64(candidate.a_ricker) / survival_ref
ENV["COTS_LARVAL_FECUNDITY"] = string(fecundity)
production(adults) = fecundity * adults * exp(-Float64(candidate.b_ricker) * adults) *
    adults^2 / (3.0^2 + adults^2)
boundary_ceiling = maximum(production.(range(0.0, 100.0; length=10001)))
ENV["COTS_EXTERNAL_OUTBREAK_PRODUCTION"] = string(boundary_ceiling)

Random.seed!(SEED)
base_scenario = ADRIA.sample(domain, 2)[1:1, :]
defaults = ADRIA.param_table(domain)
for name in names(base_scenario)
    base_scenario[!, name] .= defaults[1, name]
    if name != "allee_threshold" && hasproperty(candidate, Symbol(name))
        base_scenario[!, name] .= candidate[Symbol(name)]
    end
end
base_scenario[!, :allee_threshold] .= 3.0

# Legacy defaults are sensitivity anchors, not newly supported biological bounds.
const TREATMENTS = [
    (name="current_m2_m3", m2=Float64(candidate.m2), m3=Float64(candidate.m3)),
    (name="legacy_m2_only", m2=Float64(defaults[1, :m2]), m3=Float64(candidate.m3)),
    (name="legacy_m3_only", m2=Float64(candidate.m2), m3=Float64(defaults[1, :m3])),
    (name="legacy_m2_m3", m2=Float64(defaults[1, :m2]), m3=Float64(defaults[1, :m3])),
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
    scenario = copy(base_scenario)
    scenario[!, :m2] .= treatment.m2
    scenario[!, :m3] .= treatment.m3
    Random.seed!(SEED)
    try
        result = ADRIA.run_scenario(domain, scenario[1, :])
        years = Int.(result.cots_calendar_year)
        years == collect(1985:2024) || error("Unexpected calendar years")
        cover = dropdims(sum(result.raw; dims=(2, 3)); dims=(2, 3))
        for target in TARGETS
            reef = target.model
            stages = [aggregate_reef_density(Matrix(result.cots_log[:, channel, :]), mapping[reef])
                for channel in 1:3]
            reef_cover = aggregate_reef_density(cover, mapping[reef])
            condition = aggregate_reef_density(Matrix(result.cots_condition_log), mapping[reef])
            reef_flows = [aggregate_reef_density(Matrix(result.cots_flow_log[:, channel, :]),
                mapping[reef]) for channel in 1:length(FLOW_NAMES)]
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

metadata = Dict(
    "status" => "exploratory_not_promoted",
    "seed" => SEED,
    "fixed_allee_threshold_cots_ha" => 3.0,
    "temporal_mode" => "sampled hydrodynamic year; not historical hindcast",
    "source_assumptions" => "provisional Owen: 56 unmatched sources omitted, full nonfinite rows zeroed, 2017 survival reused after 2017",
    "stage_survival_note" => "legacy m2/m3 values are sensitivity anchors, not evidence-bounded replacements",
    "prepelagic_fecundity" => fecundity,
    "boundary_ceiling_recruits_ha_year" => boundary_ceiling,
    "observation_scale" => "frozen V1 legacy closed scale from 20261001T113401 gate",
    "julia_version" => string(VERSION),
    "adria_revision" => git_revision(ROOT),
    "cotsmod_revision" => git_revision(joinpath(ROOT, "..", "COTSMod.jl")),
    "git_status" => read(`git status --porcelain`, String),
    "treatments" => [Dict(String(k) => v for (k, v) in pairs(t)) for t in TREATMENTS],
    "inputs" => Dict(
        "candidate" => file_sha256(CANDIDATE_PATH),
        "gate_scores" => file_sha256(GATE_SCORES_PATH),
        "observations" => file_sha256(OBS_PATH),
        "v2_domain" => file_sha256(joinpath(DOMAIN_PATH, "provenance.json")),
        "owen_forcing" => file_sha256(joinpath(DOMAIN_PATH, "cots_connectivity", "datasets",
            "owen_global_2018_2023", "provenance.json")),
        "script" => file_sha256(@__FILE__),
        "mapping_script" => file_sha256(joinpath(@__DIR__, "reef_observation_mapping.jl")),
        "score_script" => file_sha256(joinpath(@__DIR__, "cots_cycle_metrics.jl")),
        "manifest" => file_sha256(joinpath(ROOT, "sandbox", "Manifest.toml")),
    ),
    "outputs" => Dict(file => file_sha256(joinpath(OUT, file)) for file in
        ("scores.csv", "trajectories.csv", "flows.csv", "failures.csv")),
)
open(joinpath(OUT, "metadata.toml"), "w") do io
    TOML.print(io, metadata; sorted=true)
end
println("Saved ", OUT)
