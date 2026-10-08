"""Bounded one-seed dynamic screen of opt-in low-coral-cover COTS mortality.

Replay the archived joint-fourfold Owen control; change only the new
low_cover_mortality switch/strength at fixed 0.10 pre-feeding prey cover.
The treatment is diagnostic, not an evidence-calibrated default.
"""

using Pkg
const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
Pkg.activate(joinpath(ROOT, "sandbox"))
cd(ROOT)

using ADRIA, CSV, DataFrames, Dates, Random, TOML
include(joinpath(@__DIR__, "cots_cycle_metrics.jl"))
include(joinpath(@__DIR__, "reef_observation_mapping.jl"))
include(joinpath(@__DIR__, "calibration_run.jl"))

const SEED = 20260930
const RUN_ID = get(ENV, "OWEN_LOW_COVER_RUN_ID",
    Dates.format(now(), "yyyymmddTHHMMSS") * "_owen_low_cover_mortality")
occursin(r"^[A-Za-z0-9_.-]+$", RUN_ID) || error("Invalid run ID")
const OUT = joinpath(@__DIR__, "runs", RUN_ID)
ispath(OUT) && error("Refusing to overwrite $OUT")
mkpath(OUT)
const V2 = joinpath(ROOT, "sandbox", "data", "Lizard_Historical_v0.2")
const SCREEN = joinpath(@__DIR__, "runs", "20261001T162000_owen_embedded_mortality_screen")
const REPLICATES = joinpath(@__DIR__, "runs", "20261001T163000_owen_embedded_replicates")
const BEST = joinpath(@__DIR__, "runs", "expanded_cotsconn_pilot_seed20260930",
    "best_summary.csv")
const GATE = joinpath(@__DIR__, "runs", "20261001T113401_lizard_domain_gate")
const OBS = joinpath(ROOT, "sandbox", "data", "reef_cots.csv")
const DIAGNOSTIC = joinpath(@__DIR__, "runs",
    "20261002T201000_owen_food_trigger_diagnostic", "trigger_audit", "metadata.json")
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
const TREATMENTS = [
    (name="joint4_reference", enabled=false, strength=0.0),
    (name="low_cover_0p5", enabled=true, strength=0.5),
    (name="low_cover_1p0", enabled=true, strength=1.0),
]

function recorded_revision(path::String, env_key::String)::String
    haskey(ENV, env_key) && return ENV[env_key]
    try
        return git_revision(path)
    catch
        error("Git revision unavailable; set $env_key from `git rev-parse HEAD`")
    end
end

screen_meta = TOML.parsefile(joinpath(SCREEN, "metadata.toml"))
fecundity = 4Float64(screen_meta["matched_fecundity_pre_pelagic_ha_year"])
boundary = 4Float64(screen_meta["matched_boundary_pre_pelagic_ha_year"])
candidate = sort(CSV.read(BEST, DataFrame), :loss)[1, :]
gate_scores = CSV.read(joinpath(GATE, "scores.csv"), DataFrame)
scale_rows = gate_scores[(gate_scores.seed .== SEED) .&
    (gate_scores.treatment .== "v1_legacy_closed") .&
    (gate_scores.observation_treatment .== "raw"), :]
fixed_scales = Dict(String(row.reef_name) => Float64(row.observation_scale)
    for row in eachrow(scale_rows))
length(fixed_scales) == length(TARGETS) || error("Missing fixed scales")

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
ENV["COTS_APPLY_LARVAL_SURVIVAL"] = "false"
ENV["COTS_SETTLEMENT_PROBABILITY"] = "1.0"
ENV["COTS_IMMIGRATION_SCALAR"] = "1.0"
ENV["COTS_INITIAL_MULTIPLIER"] = string(candidate.seed_mult)
ENV["COTS_SEED_FIRST_N"] = "10"
ENV["ADRIA_DEBUG_SEED_FIRST_N"] = "10"
ENV["COTS_LARVAL_FECUNDITY"] = string(fecundity)
ENV["COTS_EXTERNAL_OUTBREAK_PRODUCTION"] = string(boundary)
ENV["COTS_LOW_COVER_THRESHOLD"] = "0.10"
delete!(ENV, "COTS_INTERNAL_SETTLEMENT_BLACKOUT_START_YEAR")
delete!(ENV, "COTS_INTERNAL_SETTLEMENT_BLACKOUT_END_YEAR")

domain = ADRIA.load_domain(ADRIA.LizardDomain, V2, "historical")
mapping = load_reef_area_weights(domain, V2, [t.model for t in TARGETS])
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
raw_obs = CSV.read(OBS, DataFrame)
archived = CSV.read(joinpath(REPLICATES, "trajectories.csv"), DataFrame)
archived = archived[(archived.seed .== SEED) .&
    (archived.treatment .== "embedded_joint_4x"), :]
nrow(archived) == 40length(TARGETS) || error("Incomplete archived baseline")

scores = DataFrame(treatment=String[], reef_name=String[], observation_treatment=String[],
    loss=Float64[], sim_peak_count=Int[], sim_peak_years=String[],
    sim_peak_heights_cpue=String[], matched_peaks=Int[])
trajectories = DataFrame(treatment=String[], reef_name=String[], year=Int[],
    forcing_year=Int[], adults_ha=Float64[], simulated_cpue=Float64[], coral_cover=Float64[])
flows = DataFrame(treatment=String[], reef_name=String[], year=Int[])
for name in FLOW_NAMES
    flows[!, name] = Float64[]
end
food = DataFrame(treatment=String[], reef_name=String[], year=Int[],
    condition_area_mean=Float64[], adult_food_survival=Float64[],
    adult_food_deaths_ha=Float64[], juvenile_food_deaths_ha=Float64[])
failures = DataFrame(treatment=String[], error=String[])

for treatment in TREATMENTS
    ENV["COTS_LOW_COVER_MORTALITY"] = string(treatment.enabled)
    ENV["COTS_LOW_COVER_STRENGTH"] = string(treatment.strength)
    Random.seed!(SEED)
    try
        result = ADRIA.run_scenario(domain, base_scenario[1, :])
        years = Int.(result.cots_calendar_year)
        years == collect(1985:2024) || error("Unexpected model years")
        cover = dropdims(sum(result.raw; dims=(2, 3)); dims=(2, 3))
        for target in TARGETS
            reef = target.model
            m = mapping[reef]
            adults = aggregate_reef_density(Matrix(result.cots_log[:, 3, :]), m)
            coral = aggregate_reef_density(cover, m)
            maturation = aggregate_reef_density(Matrix(result.cots_flow_log[:, 5, :]), m)
            condition = aggregate_reef_density(Matrix(result.cots_condition_log), m)
            reef_flows = [aggregate_reef_density(Matrix(result.cots_flow_log[:, channel, :]), m)
                for channel in 1:length(FLOW_NAMES)]
            if treatment.name == "joint4_reference"
                reference = sort(archived[archived.reef_name .== reef, :], :year)
                years == Int.(reference.year) || error("Baseline years changed")
                Int.(result.cots_forcing_year) == Int.(reference.forcing_year) ||
                    error("Baseline hydrodynamic years changed")
                all(isapprox.(adults, Float64.(reference.adults_ha); atol=1e-12, rtol=0)) ||
                    error("Baseline adult state changed")
                all(isapprox.(coral, Float64.(reference.coral_cover); atol=1e-12, rtol=0)) ||
                    error("Baseline coral state changed")
            end
            scale = fixed_scales[reef]
            for (idx, year) in enumerate(years)
                push!(trajectories, (treatment.name, reef, year,
                    Int(result.cots_forcing_year[idx]), adults[idx], adults[idx] * scale,
                    coral[idx]))
                push!(flows, (treatment.name, reef, year,
                    (values[idx] for values in reef_flows)...))
                if idx >= 2
                    survivor = max(0.0, adults[idx] - maturation[idx])
                    available = max(0.0, (1 - Float64(candidate.m3)) * adults[idx - 1])
                    effective = available > 1e-10 ? clamp(survivor / available, 0.0, 1.0) : 1.0
                    adult_food_deaths = max(0.0, available - survivor)
                    juvenile_food_deaths = 0.0
                    for (site, weight) in zip(m.indices, m.weights)
                        before = Float64(result.cots_log[idx - 1, 3, site])
                        after = Float64(result.cots_log[idx, 3, site])
                        matured = Float64(result.cots_flow_log[idx, 5, site])
                        site_available = (1 - Float64(candidate.m3)) * before
                        site_f = site_available > 1e-10 ?
                            clamp((after - matured) / site_available, 0.0, 1.0) : 1.0
                        juvenile_food_deaths += weight *
                            Float64(result.cots_log[idx - 1, 2, site]) *
                            (1 - Float64(candidate.m2)) * (1 - site_f)
                    end
                    push!(food, (treatment.name, reef, year, condition[idx], effective,
                        adult_food_deaths, juvenile_food_deaths))
                end
            end
            obs = sort(raw_obs[(raw_obs.reef_name .== target.observed) .&
                (raw_obs.year .>= 1985) .& (raw_obs.year .<= 2024), :], :year)
            obs_years = Int.(obs.year)
            obs_values = Float64.(obs.cotsptow)
            smoothed = calendar_moving_average(obs_years, obs_values; window_years=3)
            for (name, values) in (("raw", obs_values), ("smoothed_3y", smoothed))
                score = cots_cycle_score(years, adults, obs_years, values;
                    fixed_observation_scale=scale)
                window = (years .>= score.scoring_start_year) .&
                    (years .<= score.scoring_end_year)
                peaks = detect_cots_peaks(years[window], (adults .* scale)[window];
                    smooth_window_years=3)
                push!(scores, (treatment.name, reef, name, score.total_loss,
                    score.n_sim_peaks, join(score.sim_peak_years, ';'),
                    join(round.([p.value for p in peaks]; digits=5), ';'),
                    score.n_matched_peaks))
            end
        end
        println("Completed ", treatment.name)
    catch err
        push!(failures, (treatment.name, sprint(showerror, err, catch_backtrace())))
        println(stderr, "FAILED ", treatment.name, ": ", err)
    end
    flush(stdout)
    for (file, data) in (("scores.csv", scores), ("trajectories.csv", trajectories),
        ("flows.csv", flows), ("food_diagnostics.csv", food), ("failures.csv", failures))
        CSV.write(joinpath(OUT, file), data)
    end
end

metadata = Dict(
    "status" => "bounded_low_cover_mortality_screen_not_promoted",
    "seed" => SEED,
    "hypothesis" => "independent cover-only food mortality acts more strongly in the trough than both waves",
    "fixed_allee_threshold_cots_ha" => 3.0,
    "cover_threshold_fraction" => 0.10,
    "bound_note" => "strengths 0.5 and 1.0 are diagnostic, not independently estimated field mortality",
    "fecundity_pre_pelagic_recruits_per_adult_year" => fecundity,
    "boundary_pre_pelagic_recruits_per_ha_year" => boundary,
    "reference_replay" => "adult and coral checked against archived seed-20260930 joint4 to 1e-12; hydrodynamic year exact",
    "trough_gate" => "three-year-smoothed minimum / smaller adjacent peak <=0.10 on all three two-wave reefs with seven observed peak matches",
    "julia_version" => string(VERSION),
    "adria_revision" => recorded_revision(ROOT, "OWEN_ADRIA_REVISION"),
    "cotsmod_revision" => recorded_revision(joinpath(ROOT, "..", "COTSMod.jl"),
        "OWEN_COTSMOD_REVISION"),
    "treatments" => [Dict(String(k) => v for (k, v) in pairs(t)) for t in TREATMENTS],
    "inputs" => Dict(
        "diagnostic_gate" => file_sha256(DIAGNOSTIC),
        "screen_metadata" => file_sha256(joinpath(SCREEN, "metadata.toml")),
        "reference_trajectories" => file_sha256(joinpath(REPLICATES, "trajectories.csv")),
        "candidate" => file_sha256(BEST),
        "gate_scores" => file_sha256(joinpath(GATE, "scores.csv")),
        "observations" => file_sha256(OBS),
        "script" => file_sha256(@__FILE__),
        "adapter" => file_sha256(joinpath(ROOT, "ADRIA", "src", "ecosystem", "cots.jl")),
        "scenario" => file_sha256(joinpath(ROOT, "ADRIA", "src", "scenario.jl")),
        "core" => file_sha256(joinpath(ROOT, "..", "COTSMod.jl", "src", "COTSMod.jl")),
        "owen_forcing" => file_sha256(joinpath(V2, "cots_connectivity", "datasets",
            "owen_global_2018_2023", "provenance.json")),
    ),
    "outputs" => Dict(file => file_sha256(joinpath(OUT, file)) for file in
        ("scores.csv", "trajectories.csv", "flows.csv", "food_diagnostics.csv", "failures.csv")),
)
open(joinpath(OUT, "metadata.toml"), "w") do io
    TOML.print(io, metadata; sorted=true)
end
println("Saved ", OUT)
