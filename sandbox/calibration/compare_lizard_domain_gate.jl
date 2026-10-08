"""Replicated Lizard V1/V2 and legacy/Owen peak comparison, without fitting ecology."""

using Pkg
const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
Pkg.activate(joinpath(REPO_ROOT, "sandbox"))
cd(REPO_ROOT)

using ADRIA
using CSV, DataFrames, Dates, Random, SHA, Statistics, TOML

include(joinpath(@__DIR__, "cots_cycle_metrics.jl"))
include(joinpath(@__DIR__, "reef_observation_mapping.jl"))
include(joinpath(@__DIR__, "calibration_run.jl"))

const SEEDS = parse.(Int, split(get(ENV, "LIZARD_GATE_SEEDS", "20260930,20261001,20261002"), ','))
length(SEEDS) >= 3 && length(unique(SEEDS)) == length(SEEDS) ||
    error("LIZARD_GATE_SEEDS requires at least three distinct recorded seeds")
const REFERENCE_SEED = 20260930
REFERENCE_SEED in SEEDS || error("Reference seed 20260930 must be in LIZARD_GATE_SEEDS")
const TEMPORAL_MODE = lowercase(get(ENV, "LIZARD_GATE_TEMPORAL_MODE", "cycle"))
TEMPORAL_MODE in ("cycle", "sample") || error("LIZARD_GATE_TEMPORAL_MODE must be cycle or sample")
const RUN_ID = get(ENV, "LIZARD_GATE_RUN_ID", Dates.format(now(), "yyyymmddTHHMMSS") * "_lizard_domain_gate")
occursin(r"^[A-Za-z0-9_.-]+$", RUN_ID) || error("Unsupported LIZARD_GATE_RUN_ID")
const OUT_DIR = joinpath(REPO_ROOT, "sandbox", "calibration", "runs", RUN_ID)
ispath(OUT_DIR) && error("Refusing to overwrite $OUT_DIR")
mkpath(OUT_DIR)

const V1_PATH = joinpath(REPO_ROOT, "sandbox", "data", "Lizard_Historical_v0.1")
const V2_PATH = joinpath(REPO_ROOT, "sandbox", "data", "Lizard_Historical_v0.2")
const CANDIDATE_PATH = joinpath(REPO_ROOT, "sandbox", "calibration", "runs",
    "expanded_cotsconn_pilot_seed20260930", "best_summary.csv")
isfile(CANDIDATE_PATH) || error("Missing frozen biological candidate")
best = sort(CSV.read(CANDIDATE_PATH, DataFrame), :loss)[1, :]

# The old scientific control remains untouched.  Only this named experiment
# enables separated production, annual forcing, and the evidence boundary.
ENV["ADRIA_COTS_CONNECTIVITY_MODE"] = "cots"
ENV["COTS_CONNECTIVITY_TEMPORAL_MODE"] = TEMPORAL_MODE
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

const TARGETS = [
    (model="Lizard Island Reef", observed="Lizard Isles"),
    (model="MacGillivray Reef", observed="Macgillivray Reef"),
    (model="North Direction Reef", observed="North Direction Island"),
    (model="Eyrie Reef", observed="Eyrie Reef"),
]
const MODEL_REEFS = [target.model for target in TARGETS]
const OBSERVATIONS_PATH = joinpath(REPO_ROOT, "sandbox", "data", "reef_cots.csv")
observations = CSV.read(OBSERVATIONS_PATH, DataFrame)
observed = Dict{String,NamedTuple}()
observed_rows = DataFrame(reef_name=String[], year=Int[], raw_cpue=Float64[],
                          smoothed_3y_cpue=Float64[], tows=Int[])
for target in TARGETS
    rows = sort(observations[(observations.reef_name .== target.observed) .&
                              (observations.year .>= 1985) .&
                              (observations.year .<= 2024), :], :year)
    length(unique(rows.year)) == nrow(rows) || error("Duplicate survey-year rows: $(target.observed)")
    years = Int.(rows.year)
    raw = Float64.(rows.cotsptow)
    all(isfinite, raw) && all(raw .>= 0.0) || error("Invalid observed CPUE")
    smooth = calendar_moving_average(years, raw; window_years=3)
    observed[target.model] = (years=years, raw=raw, smooth=smooth)
    for index in eachindex(years)
        push!(observed_rows, (target.model, years[index], raw[index], smooth[index], Int(rows.tows[index])))
    end
end
CSV.write(joinpath(OUT_DIR, "observations.csv"), observed_rows)

function load_treatment_domain(domain_path, dataset)
    ENV["ADRIA_COTS_CONNECTIVITY_DATASET"] = dataset
    domain = ADRIA.load_domain(ADRIA.LizardDomain, domain_path, "historical")
    mapping = load_reef_area_weights(domain, domain_path, MODEL_REEFS)
    return (domain=domain, mapping=mapping)
end

domains = Dict(
    "v1_legacy" => load_treatment_domain(V1_PATH, "legacy"),
    "v2_legacy" => load_treatment_domain(V2_PATH, "legacy"),
    "v2_owen" => load_treatment_domain(V2_PATH, "owen_global_2018_2023"),
)

weight_rows = DataFrame(domain=String[], reef_name=String[], site_id=String[],
                        site_index=Int[], area_m2=Float64[], within_reef_weight=Float64[])
area_rows = DataFrame(domain=String[], reef_name=String[], site_count=Int[], total_area_m2=Float64[])
for (key, entry) in sort(collect(domains); by=first)
    for reef_name in MODEL_REEFS
        mapping = entry.mapping[reef_name]
        push!(area_rows, (key, reef_name, length(mapping.indices), mapping.area_m2))
        for (index, weight) in zip(mapping.indices, mapping.weights)
            push!(weight_rows, (key, reef_name, String(entry.domain.loc_ids[index]),
                                index, Float64(entry.domain.loc_data.area[index]), weight))
        end
    end
end
CSV.write(joinpath(OUT_DIR, "reef_observation_weights.csv"), weight_rows)
CSV.write(joinpath(OUT_DIR, "reef_area_summary.csv"), area_rows)

reference_survival = mean(domains["v2_legacy"].domain.cots_forcing.source_survival)
reference_survival > 0.0 || error("Invalid V2 legacy reference survival")
fecundity = Float64(best.a_ricker) / reference_survival
ENV["COTS_LARVAL_FECUNDITY"] = string(fecundity)
allee_threshold = 3.0
b_ricker = Float64(best.b_ricker)
production(adults) = fecundity * adults * exp(-b_ricker * adults) *
    adults^2 / (allee_threshold^2 + adults^2)
boundary_ceiling = maximum(production.(range(0.0, 100.0; length=10001)))

Random.seed!(REFERENCE_SEED)
scenario = ADRIA.sample(domains["v1_legacy"].domain, 2)[1:1, :]
defaults = ADRIA.param_table(domains["v1_legacy"].domain)
for name in names(scenario)
    scenario[!, name] .= defaults[1, name]
    if name != "allee_threshold" && hasproperty(best, Symbol(name))
        scenario[!, name] .= best[Symbol(name)]
    end
end
scenario[!, :allee_threshold] .= allee_threshold

const TREATMENTS = [
    (name="v1_legacy_closed", domain="v1_legacy", dataset="legacy", boundary=false),
    (name="v2_legacy_closed", domain="v2_legacy", dataset="legacy", boundary=false),
    (name="v1_legacy_boundary", domain="v1_legacy", dataset="legacy", boundary=true),
    (name="v2_legacy_boundary", domain="v2_legacy", dataset="legacy", boundary=true),
    (name="v2_owen_closed", domain="v2_owen", dataset="owen_global_2018_2023", boundary=false),
    (name="v2_owen_boundary", domain="v2_owen", dataset="owen_global_2018_2023", boundary=true),
]

scores = DataFrame(seed=Int[], treatment=String[], domain=String[], dataset=String[],
    boundary=Bool[], reef_name=String[], observation_treatment=String[],
    loss=Float64[], sim_peak_count=Int[], sim_peak_years=String[],
    sim_peak_heights_cpue=String[], sim_peak_prominences_cpue=String[],
    obs_peak_count=Int[], obs_peak_years=String[], obs_peak_heights_cpue=String[],
    matched_peaks=Int[], peak_timing_penalty=Float64[], amplitude_penalty=Float64[],
    period_penalty=Float64[], observation_scale=Float64[],
    mean_external_pelagic_ha=Float64[])
trajectories = DataFrame(seed=Int[], treatment=String[], domain=String[],
    dataset=String[], boundary=Bool[], reef_name=String[], year=Int[],
    adults_ha=Float64[], simulated_cpue=Float64[], coral_cover=Float64[],
    external_pelagic_ha=Float64[], forcing_year=Int[])
failures = DataFrame(seed=Int[], treatment=String[], error=String[])
fixed_scales = Dict{String,Float64}()

for seed in vcat(REFERENCE_SEED, filter(!=(REFERENCE_SEED), SEEDS))
    for treatment in TREATMENTS
        entry = domains[treatment.domain]
        ENV["ADRIA_COTS_CONNECTIVITY_DATASET"] = treatment.dataset
        ENV["COTS_EXTERNAL_OUTBREAK_PRODUCTION"] = string(treatment.boundary ? boundary_ceiling : 0.0)
        ENV["COTS_CONNECTIVITY_SEED"] = string(seed)
        Random.seed!(seed)
        try
            result = ADRIA.run_scenario(entry.domain, scenario[1, :])
            years = Int.(result.cots_calendar_year)
            years == collect(1985:2024) || error("Unexpected model years")
            coral = dropdims(sum(result.raw; dims=(2, 3)); dims=(2, 3))
            adults_by_site = Matrix(result.cots_log[:, 3, :])
            external_by_site = Matrix(result.cots_flow_log[:, 12, :])
            all(isfinite, adults_by_site) && minimum(adults_by_site) >= 0.0 ||
                error("Adult COTS state is nonfinite or negative")
            for target in TARGETS
                reef_name = target.model
                mapping = entry.mapping[reef_name]
                adults = aggregate_reef_density(adults_by_site, mapping)
                covers = aggregate_reef_density(coral, mapping)
                external = aggregate_reef_density(external_by_site, mapping)
                obs = observed[reef_name]
                if seed == REFERENCE_SEED && treatment.name == "v1_legacy_closed"
                    provisional = cots_cycle_score(years, adults, obs.years, obs.raw)
                    fixed_scales[reef_name] = provisional.observation_scale
                end
                scale = fixed_scales[reef_name]
                for index in eachindex(years)
                    push!(trajectories, (seed, treatment.name, treatment.domain,
                        treatment.dataset, treatment.boundary, reef_name, years[index],
                        adults[index], adults[index] * scale, covers[index],
                        external[index], Int(result.cots_forcing_year[index])))
                end
                for (obs_name, obs_values) in (("raw", obs.raw), ("smoothed_3y", obs.smooth))
                    score = cots_cycle_score(
                        years, adults, obs.years, obs_values;
                        fixed_observation_scale=scale
                    )
                    scoring_window = (years .>= score.scoring_start_year) .&
                        (years .<= score.scoring_end_year)
                    sim_peaks = detect_cots_peaks(
                        years[scoring_window], (adults .* scale)[scoring_window];
                        smooth_window_years=3
                    )
                    obs_peaks = detect_cots_peaks(obs.years, obs_values)
                    length(sim_peaks) == score.n_sim_peaks || error("Peak-report count mismatch")
                    push!(scores, (seed, treatment.name, treatment.domain,
                        treatment.dataset, treatment.boundary, reef_name, obs_name,
                        score.total_loss, score.n_sim_peaks,
                        join(score.sim_peak_years, ';'),
                        join(round.([p.value for p in sim_peaks]; digits=5), ';'),
                        join(round.([p.prominence for p in sim_peaks]; digits=5), ';'),
                        score.n_obs_peaks, join(score.obs_peak_years, ';'),
                        join(round.([p.value for p in obs_peaks]; digits=5), ';'),
                        score.n_matched_peaks, score.peak_timing_penalty,
                        score.amplitude_penalty, score.period_penalty,
                        scale, mean(external)))
                end
            end
            println("Completed seed=", seed, " treatment=", treatment.name)
            flush(stdout)
        catch err
            push!(failures, (seed, treatment.name, sprint(showerror, err, catch_backtrace())))
            println(stderr, "FAILED seed=", seed, " treatment=", treatment.name, ": ", err)
            if seed == REFERENCE_SEED && treatment.name == "v1_legacy_closed"
                CSV.write(joinpath(OUT_DIR, "failures.csv"), failures)
                rethrow()
            end
        end
        CSV.write(joinpath(OUT_DIR, "scores.csv"), scores)
        CSV.write(joinpath(OUT_DIR, "trajectories.csv"), trajectories)
        CSV.write(joinpath(OUT_DIR, "failures.csv"), failures)
    end
end

function sha(path)
    isfile(path) || error("Missing provenance input: $path")
    return file_sha256(path)
end

project_path = joinpath(REPO_ROOT, "sandbox", "Project.toml")
manifest_path = joinpath(REPO_ROOT, "sandbox", "Manifest.toml")
metadata = Dict(
    "design" => "3-seed paired V1/V2 legacy and V2 Owen, each closed/evidence boundary",
    "temporal_mode" => TEMPORAL_MODE,
    "seed_interpretation" => TEMPORAL_MODE == "cycle" ?
        "cycle is deterministic; seed replicates are a reproducibility check, not independent outcomes" :
        "sample selects a hydrodynamic year for each model year using the recorded seed; this is a transport sensitivity, not a historical hindcast",
    "status" => "exploratory_not_promoted",
    "seeds" => SEEDS,
    "reference_seed_for_observation_scale" => REFERENCE_SEED,
    "observation_aggregation" => "within-reef polygon-area-weighted site COTS ha^-1",
    "observation_conversion" => "per-reef least-squares scale fitted only on V1 legacy closed reference seed, then frozen",
    "observation_treatments" => ["raw", "smoothed_3y"],
    "smoothing" => "symmetric calendar +/-1 year mean of available survey CPUE; no interpolation across missing years",
    "scoring_years" => "common observation/model window ending 2024",
    "allee_threshold_cots_ha" => allee_threshold,
    "prepelagic_fecundity" => fecundity,
    "boundary_ceiling_recruits_ha_year" => boundary_ceiling,
    "owen_policy_status" => "provisional: zero full nonfinite rows, exclude 56 unmatched sources, reuse 2017 survival",
    "julia_version" => string(VERSION),
    "adria_revision" => git_revision(REPO_ROOT),
    "cotsmod_revision" => git_revision(normpath(joinpath(REPO_ROOT, "..", "COTSMod.jl"))),
    "git_status" => read(`git status --porcelain`, String),
    "inputs" => Dict(
        "candidate" => sha(CANDIDATE_PATH),
        "observations" => sha(OBSERVATIONS_PATH),
        "v1_site_map" => sha(joinpath(V1_PATH, "site_to_reef.csv")),
        "v2_domain_provenance" => sha(joinpath(V2_PATH, "provenance.json")),
        "v1_legacy_forcing" => sha(joinpath(V1_PATH, "cots_connectivity", "annual_forcing", "provenance.toml")),
        "v2_legacy_forcing" => sha(joinpath(V2_PATH, "cots_connectivity", "annual_forcing", "provenance.toml")),
        "v2_owen_forcing" => sha(joinpath(V2_PATH, "cots_connectivity", "datasets", "owen_global_2018_2023", "provenance.json")),
        "project" => sha(project_path),
        "manifest" => sha(manifest_path),
        "comparison_script" => sha(@__FILE__),
        "observation_mapping_script" => sha(joinpath(@__DIR__, "reef_observation_mapping.jl")),
        "score_script" => sha(joinpath(@__DIR__, "cots_cycle_metrics.jl")),
    ),
    "outputs" => Dict(filename => sha(joinpath(OUT_DIR, filename)) for filename in
        ["observations.csv", "reef_observation_weights.csv", "reef_area_summary.csv",
         "scores.csv", "trajectories.csv", "failures.csv"]),
)
open(joinpath(OUT_DIR, "metadata.toml"), "w") do io
    TOML.print(io, metadata; sorted=true)
end
println("Wrote replicated Lizard domain gate run to ", OUT_DIR)
