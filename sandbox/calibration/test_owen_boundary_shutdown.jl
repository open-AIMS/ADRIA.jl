"""Counterfactual: zero only 2000-10 external supply to test inter-wave persistence.

This is an oracle diagnostic, not an evidence-derived historical boundary.
All ecological parameters, Owen matrices, source-survival assumption, Allee
threshold, seed, and observation scale match embedded_joint_4x.
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
const RUN_ID = get(ENV, "OWEN_SHUTDOWN_RUN_ID",
    Dates.format(now(), "yyyymmddTHHMMSS") * "_owen_boundary_shutdown")
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
ENV["COTS_INITIAL_MULTIPLIER"] = string(candidate.seed_mult)
ENV["COTS_SEED_FIRST_N"] = "10"
ENV["ADRIA_DEBUG_SEED_FIRST_N"] = "10"
ENV["COTS_LARVAL_FECUNDITY"] = string(fecundity)
ENV["COTS_EXTERNAL_OUTBREAK_PRODUCTION"] = string(boundary)

domain = ADRIA.load_domain(ADRIA.LizardDomain, V2, "historical")
mapping = load_reef_area_weights(domain, V2, [t.model for t in TARGETS])
forcing = domain.cots_forcing
for year in 2000:2010
    index = findfirst(==(year), forcing.outbreak_years)
    isnothing(index) && error("Missing outbreak boundary year $year")
    forcing.external_outbreak_potential[:, index, :] .= 0.0
    forcing.external_outbreak_pelagic[:, index, :] .= 0.0
end

Random.seed!(SEED)
scenario = ADRIA.sample(domain, 2)[1:1, :]
defaults = ADRIA.param_table(domain)
for name in names(scenario)
    scenario[!, name] .= defaults[1, name]
    if name != "allee_threshold" && hasproperty(candidate, Symbol(name))
        scenario[!, name] .= candidate[Symbol(name)]
    end
end
scenario[!, :allee_threshold] .= 3.0
Random.seed!(SEED)
result = ADRIA.run_scenario(domain, scenario[1, :])
years = Int.(result.cots_calendar_year)
years == collect(1985:2024) || error("Unexpected model years")
for t in findall(y -> 2000 <= y <= 2010, years)
    all(iszero, result.cots_flow_log[t, 11, :]) || error("Boundary potential not zero")
    all(iszero, result.cots_flow_log[t, 12, :]) || error("Boundary pelagic not zero")
end

raw_obs = CSV.read(OBS, DataFrame)
cover = dropdims(sum(result.raw; dims=(2, 3)); dims=(2, 3))
scores = DataFrame(reef_name=String[], observation_treatment=String[],
    loss=Float64[], sim_peak_count=Int[], sim_peak_years=String[],
    sim_peak_heights_cpue=String[], matched_peaks=Int[])
trajectories = DataFrame(reef_name=String[], year=Int[], forcing_year=Int[],
    adults_ha=Float64[], simulated_cpue=Float64[], coral_cover=Float64[])
flows = DataFrame(reef_name=String[], year=Int[])
for name in FLOW_NAMES
    flows[!, name] = Float64[]
end
for target in TARGETS
    reef = target.model
    m = mapping[reef]
    adults = aggregate_reef_density(Matrix(result.cots_log[:, 3, :]), m)
    coral = aggregate_reef_density(cover, m)
    reef_flows = [aggregate_reef_density(Matrix(result.cots_flow_log[:, channel, :]), m)
        for channel in 1:length(FLOW_NAMES)]
    scale = fixed_scales[reef]
    for (i, year) in enumerate(years)
        push!(trajectories, (reef, year, Int(result.cots_forcing_year[i]),
            adults[i], adults[i] * scale, coral[i]))
        push!(flows, (reef, year, (values[i] for values in reef_flows)...))
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
        push!(scores, (reef, name, score.total_loss, score.n_sim_peaks,
            join(score.sim_peak_years, ';'),
            join(round.([p.value for p in peaks]; digits=5), ';'),
            score.n_matched_peaks))
    end
end
CSV.write(joinpath(OUT, "scores.csv"), scores)
CSV.write(joinpath(OUT, "trajectories.csv"), trajectories)
CSV.write(joinpath(OUT, "flows.csv"), flows)

metadata = Dict(
    "status" => "oracle_counterfactual_not_promoted",
    "seed" => SEED,
    "treatment" => "embedded_joint_4x_boundary_zero_2000_2010",
    "boundary_zero_years" => collect(2000:2010),
    "intervention" => "set potential and pelagic external coefficients to zero for all sinks and hydrodynamic years in the window",
    "fixed_allee_threshold_cots_ha" => 3.0,
    "effective_fecundity_ha_year" => fecundity,
    "effective_boundary_ha_year" => boundary,
    "julia_version" => string(VERSION),
    "adria_revision" => git_revision(ROOT),
    "cotsmod_revision" => git_revision(joinpath(ROOT, "..", "COTSMod.jl")),
    "inputs" => Dict(
        "screen_metadata" => file_sha256(joinpath(SCREEN, "metadata.toml")),
        "replicate_trajectories" => file_sha256(joinpath(REPLICATES, "trajectories.csv")),
        "candidate" => file_sha256(BEST),
        "gate_scores" => file_sha256(joinpath(GATE, "scores.csv")),
        "observations" => file_sha256(OBS),
        "script" => file_sha256(@__FILE__),
        "owen_forcing" => file_sha256(joinpath(V2, "cots_connectivity", "datasets",
            "owen_global_2018_2023", "provenance.json")),
    ),
    "outputs" => Dict(file => file_sha256(joinpath(OUT, file)) for file in
        ("scores.csv", "trajectories.csv", "flows.csv")),
)
open(joinpath(OUT, "metadata.toml"), "w") do io
    TOML.print(io, metadata; sorted=true)
end
println("Saved ", OUT)
