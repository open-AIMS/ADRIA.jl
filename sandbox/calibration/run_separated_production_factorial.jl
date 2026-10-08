using Pkg
const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
Pkg.activate(joinpath(REPO_ROOT, "sandbox"))
cd(REPO_ROOT)

using ADRIA
using CSV, DataFrames, Dates, Random, SHA, Statistics, TOML

include(joinpath(@__DIR__, "cots_cycle_metrics.jl"))
include(joinpath(@__DIR__, "calibration_run.jl"))

seed = parse(Int, get(ENV, "COTS_FACTORIAL_SEED", "20260929"))
run_id = get(
    ENV,
    "COTS_FACTORIAL_RUN_ID",
    Dates.format(now(), "yyyymmddTHHMMSS") * "_separated_production_seed" * string(seed)
)
occursin(r"^[A-Za-z0-9_.-]+$", run_id) || error("Unsupported COTS_FACTORIAL_RUN_ID")
out_dir = joinpath(REPO_ROOT, "sandbox", "calibration", "runs", run_id)
isdir(out_dir) && !isempty(readdir(out_dir)) && error(
    "Run directory already contains outputs; choose a new COTS_FACTORIAL_RUN_ID: $out_dir"
)
mkpath(out_dir)

ENV["ADRIA_COTS_CONNECTIVITY_MODE"] = "cots"
ENV["COTS_CONNECTIVITY_TEMPORAL_MODE"] = "cycle"
ENV["COTS_CONNECTIVITY_SEED"] = string(seed)
ENV["COTS_ALLEE_THRESHOLD"] = "3.0"
ENV["COTS_CALENDAR_START_YEAR"] = "1985"
ENV["COTS_EXTERNAL_PULSE"] = "false"
ENV["COTS_EXTERNAL_SOURCE_DENSITY"] = "0.0"
ENV["COTS_JUVENILE_STORAGE"] = "false"
ENV["COTS_HABITAT_MEDIATION"] = "false"
ENV["COTS_SEED_FIRST_N"] = get(ENV, "COTS_SEED_FIRST_N", "10")
ENV["ADRIA_DEBUG_SEED_FIRST_N"] = ENV["COTS_SEED_FIRST_N"]

domain_version = get(ENV, "LIZARD_DOMAIN_VERSION", "Lizard_Historical_v0.1")
occursin(r"^[A-Za-z0-9_.-]+$", domain_version) || error("Unsupported LIZARD_DOMAIN_VERSION")
domain_path = joinpath(REPO_ROOT, "sandbox", "data", domain_version)
cots_dataset = get(ENV, "ADRIA_COTS_CONNECTIVITY_DATASET", "legacy")
occursin(r"^[A-Za-z0-9_.-]+$", cots_dataset) || error("Unsupported ADRIA_COTS_CONNECTIVITY_DATASET")
cots_data_path = cots_dataset == "legacy" ?
    joinpath(domain_path, "cots_connectivity") :
    joinpath(domain_path, "cots_connectivity", "datasets", cots_dataset)
dom = ADRIA.load_domain(ADRIA.LizardDomain, domain_path, "historical")
isnothing(dom.cots_forcing) && error(
    "Annual COTS forcing is absent; run sandbox/domain_building/build_lizard_cots_connectivity.jl"
)

requested_run_dir = get(ENV, "BBO_RUN_DIR", "")
if isempty(requested_run_dir) && haskey(ENV, "BBO_RUN_ID")
    requested_run_dir = joinpath(
        REPO_ROOT, "sandbox", "calibration", "runs", ENV["BBO_RUN_ID"]
    )
end
artifact_dir = isempty(requested_run_dir) ?
    joinpath(REPO_ROOT, "sandbox", "data") : normpath(requested_run_dir)
summary_path = joinpath(
    artifact_dir, isempty(requested_run_dir) ? "bbo_best_summary.csv" : "best_summary.csv"
)
isfile(summary_path) || (summary_path = joinpath(
    artifact_dir,
    isempty(requested_run_dir) ? "bbo_evaluated_candidates.csv" : "evaluated_candidates.csv"
))
isfile(summary_path) || error(
    "No calibrated candidate was found; set BBO_RUN_DIR or BBO_RUN_ID"
)
best = sort(CSV.read(summary_path, DataFrame), :loss)[1, :]

scenario = ADRIA.sample(dom, 2)[1:1, :]
defaults = ADRIA.param_table(dom)
for name in names(scenario)
    scenario[!, name] .= defaults[1, name]
    if name != "allee_threshold" && hasproperty(best, Symbol(name))
        scenario[!, name] .= best[Symbol(name)]
    end
end
hasproperty(best, :seed_mult) &&
    (ENV["COTS_INITIAL_MULTIPLIER"] = string(best.seed_mult))

reference_survival = mean(dom.cots_forcing.source_survival)
reference_survival > 0.0 || error("Mean source survival must be positive")
legacy_a_ricker = Float64(best.a_ricker)
larval_fecundity = legacy_a_ricker / reference_survival
allee_threshold = 3.0
b_ricker = Float64(best.b_ricker)
production(adults) = larval_fecundity * adults * exp(-b_ricker * adults) *
    adults^2 / (allee_threshold^2 + adults^2)
allee_source_production = production(allee_threshold)
adult_grid = range(0.0, 100.0; length=10001)
ricker_peak_source_production = maximum(production.(adult_grid))
ricker_peak_adult_density = adult_grid[argmax(production.(adult_grid))]

treatments = [
    (
        name="legacy_annual_internal",
        temporal="cycle",
        separated=false,
        larval_fecundity=legacy_a_ricker,
        survival=false,
        external_production=0.0,
        settlement_probability=1.0
    ),
    (
        name="separated_scaled_mean_closed",
        temporal="mean",
        separated=true,
        larval_fecundity=larval_fecundity,
        survival=true,
        external_production=0.0,
        settlement_probability=1.0
    ),
    (
        name="separated_evidence_mean_allee_source",
        temporal="mean",
        separated=true,
        larval_fecundity=larval_fecundity,
        survival=true,
        external_production=allee_source_production,
        settlement_probability=1.0
    ),
    (
        name="separated_evidence_mean_ricker_peak_source",
        temporal="mean",
        separated=true,
        larval_fecundity=larval_fecundity,
        survival=true,
        external_production=ricker_peak_source_production,
        settlement_probability=1.0
    ),
    (
        name="separated_evidence_cycle_ricker_peak_source",
        temporal="cycle",
        separated=true,
        larval_fecundity=larval_fecundity,
        survival=true,
        external_production=ricker_peak_source_production,
        settlement_probability=1.0
    )
]

site_to_reef = CSV.read(joinpath(domain_path, "site_to_reef.csv"), DataFrame)
site_to_reef.reef_name_clean = [split(name, " (")[1] for name in site_to_reef.reef_name]
target_map = Dict(
    "Lizard Island Reef" => "Lizard Isles",
    "MacGillivray Reef" => "Macgillivray Reef",
    "North Direction Reef" => "North Direction Island",
    "Eyrie Reef" => "Eyrie Reef"
)
observations = CSV.read(joinpath(REPO_ROOT, "sandbox", "data", "reef_cots.csv"), DataFrame)

trajectories = DataFrame(
    treatment=String[], reef_name=String[], year=Int[], recruits=Float64[],
    juveniles=Float64[], adults=Float64[], coral_cover=Float64[]
)
flows = DataFrame(
    treatment=String[], reef_name=String[], year=Int[], forcing_year=Int[],
    local_fecundity=Float64[], background_immigration=Float64[],
    internal_immigration=Float64[], external_immigration=Float64[],
    maturation=Float64[], retained_juveniles=Float64[], settlement_gate=Float64[],
    local_retention=Float64[], pelagic_survivors=Float64[],
    settled_recruits=Float64[], external_potential=Float64[],
    external_pelagic=Float64[]
)
scores = DataFrame(
    treatment=String[], reef_name=String[], loss=Float64[], sim_peak_count=Int[],
    sim_peak_years=String[], sim_peak_heights_cpue=String[],
    sim_peak_prominences_cpue=String[], observed_peak_years=String[],
    observed_peak_heights_cpue=String[], observed_peak_prominences_cpue=String[],
    observed_peak_count=Int[], matched_peaks=Int[], observation_scale=Float64[],
    peak_timing_penalty=Float64[], amplitude_penalty=Float64[]
)

for treatment in treatments
    ENV["COTS_CONNECTIVITY_TEMPORAL_MODE"] = treatment.temporal
    ENV["COTS_SEPARATED_RECRUITMENT"] = string(treatment.separated)
    ENV["COTS_LARVAL_FECUNDITY"] = string(treatment.larval_fecundity)
    ENV["COTS_SETTLEMENT_PROBABILITY"] = string(treatment.settlement_probability)
    ENV["COTS_APPLY_LARVAL_SURVIVAL"] = string(treatment.survival)
    ENV["COTS_EXTERNAL_OUTBREAK_PRODUCTION"] = string(treatment.external_production)
    Random.seed!(seed)
    result = ADRIA.run_scenario(dom, scenario[1, :])
    years = result.cots_calendar_year
    first(years) == 1985 && last(years) == 2024 ||
        error("Unexpected historical calendar mapping: $(first(years))-$(last(years))")
    coral_cover = dropdims(sum(result.raw, dims=(2, 3)), dims=(2, 3))

    for reef in keys(target_map)
        site_indices = findall(site_to_reef.reef_name_clean .== reef)
        isempty(site_indices) && continue
        reef_adults = [mean(result.cots_log[t, 3, site_indices]) for t in eachindex(years)]
        for t in eachindex(years)
            push!(trajectories, (
                treatment.name, reef, years[t],
                mean(result.cots_log[t, 1, site_indices]),
                mean(result.cots_log[t, 2, site_indices]), reef_adults[t],
                mean(coral_cover[t, site_indices])
            ))
            flow_values = [
                mean(result.cots_flow_log[t, channel, site_indices]) for channel in 1:12
            ]
            push!(flows, (
                treatment.name, reef, years[t], result.cots_forcing_year[t],
                flow_values...
            ))
        end

        observed = observations[observations.reef_name .== target_map[reef], :]
        score = cots_cycle_score(
            years, reef_adults, Int.(observed.year), Float64.(observed.cotsptow)
        )
        sim_yearly, sim_raw = yearly_series(years, reef_adults)
        obs_yearly, obs_raw = yearly_series(
            Int.(observed.year), Float64.(observed.cotsptow)
        )
        sim_yearly, sim_raw, obs_yearly, obs_raw, _, _ = _common_window(
            sim_yearly, sim_raw, obs_yearly, obs_raw
        )
        sim_peaks = detect_cots_peaks(
            sim_yearly, score.observation_scale .* sim_raw; smooth_window_years=3
        )
        obs_peaks = detect_cots_peaks(obs_yearly, obs_raw; smooth_window_years=1)
        push!(scores, (
            treatment.name, reef, score.total_loss, score.n_sim_peaks,
            join(score.sim_peak_years, ';'),
            join(round.(getfield.(sim_peaks, :value); digits=6), ';'),
            join(round.(getfield.(sim_peaks, :prominence); digits=6), ';'),
            join(score.obs_peak_years, ';'),
            join(round.(getfield.(obs_peaks, :value); digits=6), ';'),
            join(round.(getfield.(obs_peaks, :prominence); digits=6), ';'),
            score.n_obs_peaks, score.n_matched_peaks, score.observation_scale,
            score.peak_timing_penalty, score.amplitude_penalty
        ))
    end
end

CSV.write(joinpath(out_dir, "mechanism_trajectories.csv"), trajectories)
CSV.write(joinpath(out_dir, "mechanism_flows.csv"), flows)
CSV.write(joinpath(out_dir, "mechanism_scores.csv"), scores)

score_summary = combine(
    groupby(scores, :treatment),
    :loss => mean => :mean_loss,
    :sim_peak_count => sum => :detected_peaks,
    :matched_peaks => sum => :matched_peaks,
    :peak_timing_penalty => mean => :mean_peak_timing_penalty,
    :amplitude_penalty => mean => :mean_amplitude_penalty
)
flow_summary = combine(
    groupby(flows, :treatment),
    :local_fecundity => mean => :mean_potential_fecundity,
    :background_immigration => mean => :mean_background_immigration,
    :local_retention => mean => :mean_local_retention,
    :internal_immigration => mean => :mean_internal_immigration,
    :external_potential => mean => :mean_external_potential,
    :external_pelagic => mean => :mean_external_pelagic,
    :external_immigration => mean => :mean_external_settled,
    :settled_recruits => mean => :mean_total_settled_recruits
)
mechanism_summary = leftjoin(score_summary, flow_summary; on=:treatment)
CSV.write(joinpath(out_dir, "mechanism_summary.csv"), mechanism_summary)

forcing_dir = joinpath(cots_data_path, "annual_forcing")
metadata = Dict(
    "run_id" => run_id,
    "seed" => seed,
    "julia_version" => string(VERSION),
    "manifest_sha256" => file_sha256(joinpath(REPO_ROOT, "sandbox", "Manifest.toml")),
    "fixed_allee_threshold_cots_per_ha" => allee_threshold,
    "reference_larval_survival" => reference_survival,
    "legacy_a_ricker" => legacy_a_ricker,
    "pre_pelagic_larval_fecundity" => larval_fecundity,
    "scale_preserving_definition" => "larval_fecundity = legacy a_ricker / mean q3baseline Lizard source survival",
    "outbreak_event_definition" => "probability COTS CPUE exceeds 0.22 COTS per tow",
    "outbreak_probability_interpretation" => "continuous expected source activity, not a binary probability cutoff",
    "allee_source_production_cots_per_ha_per_year" => allee_source_production,
    "ricker_peak_source_production_cots_per_ha_per_year" =>
        ricker_peak_source_production,
    "ricker_peak_adult_density_cots_per_ha" => ricker_peak_adult_density,
    "candidate_source" => summary_path,
    "adria_revision" => git_revision(REPO_ROOT),
    "cotsmod_revision" => git_revision(normpath(joinpath(REPO_ROOT, "..", "COTSMod.jl"))),
    "inputs" => Dict(
        "candidate" => file_sha256(summary_path),
        "observations" => file_sha256(joinpath(REPO_ROOT, "sandbox", "data", "reef_cots.csv")),
        "outbreak_predictions" => file_sha256(
            joinpath(REPO_ROOT, "sandbox", "data", "gbrPredsAdj_20262408.csv")
        ),
        "forcing_provenance" => file_sha256(joinpath(forcing_dir, "provenance.toml")),
        "outbreak_boundary" => file_sha256(
            joinpath(forcing_dir, "external_outbreak_boundary.csv")
        ),
        "factorial_script" => file_sha256(@__FILE__),
        "cycle_metrics" => file_sha256(joinpath(@__DIR__, "cots_cycle_metrics.jl")),
        "adria_scenario" => file_sha256(joinpath(REPO_ROOT, "ADRIA", "src", "scenario.jl")),
        "adria_cots_adapter" => file_sha256(
            joinpath(REPO_ROOT, "ADRIA", "src", "ecosystem", "cots.jl")
        ),
        "cotsmod_core" => file_sha256(
            joinpath(REPO_ROOT, "..", "COTSMod.jl", "src", "COTSMod.jl")
        )
    ),
    "treatments" => [Dict(
        "name" => treatment.name,
        "temporal_connectivity" => treatment.temporal,
        "separated_recruitment" => treatment.separated,
        "larval_fecundity" => treatment.larval_fecundity,
        "source_larval_survival" => treatment.survival,
        "external_outbreak_production" => treatment.external_production,
        "settlement_probability" => treatment.settlement_probability
    ) for treatment in treatments],
    "promotion_gate" =>
        "at least two plausible simulated peaks with improved timing and amplitude"
)
open(joinpath(out_dir, "run_metadata.toml"), "w") do io
    TOML.print(io, metadata; sorted=true)
end

println("Separated-production factorial written to $out_dir")
show(
    stdout,
    MIME("text/plain"),
    combine(groupby(scores, :treatment), :loss => mean, :sim_peak_count => mean)
)
println()
