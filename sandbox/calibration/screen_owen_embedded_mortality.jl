"""Bounded peak-capacity screen assuming Owen transport embeds larval mortality.

Hypothesis: removing the separate RME source multiplier, while matching the
control's mean effective transport intensity, may permit distinct multi-reef
outbreak peaks under existing production and boundary parameters. No core
equation, Allee threshold, observation scale, or mortality parameter changes.
The common rescaling is a diagnostic normalization, not a biological estimate.
Failure: higher production only makes one broad peak, wrong timing, or coral
collapse. Rollback: leave this opt-in assumption off; no default changes.
"""

using Pkg
const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
Pkg.activate(joinpath(ROOT, "sandbox"))
cd(ROOT)

using ADRIA, CSV, DataFrames, Dates, Random, Statistics, TOML
include(joinpath(@__DIR__, "cots_cycle_metrics.jl"))
include(joinpath(@__DIR__, "reef_observation_mapping.jl"))
include(joinpath(@__DIR__, "calibration_run.jl"))

const SEED = 20260930
const RUN_ID = get(ENV, "OWEN_EMBEDDED_RUN_ID",
    Dates.format(now(), "yyyymmddTHHMMSS") * "_owen_embedded_mortality_screen")
occursin(r"^[A-Za-z0-9_.-]+$", RUN_ID) || error("Invalid run ID")
const OUT = joinpath(@__DIR__, "runs", RUN_ID)
ispath(OUT) && error("Refusing to overwrite $OUT")
mkpath(OUT)
const V2 = joinpath(ROOT, "sandbox", "data", "Lizard_Historical_v0.2")
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

candidate = sort(CSV.read(BEST, DataFrame), :loss)[1, :]
gate_scores = CSV.read(joinpath(GATE, "scores.csv"), DataFrame)
scale_rows = gate_scores[(gate_scores.seed .== SEED) .&
    (gate_scores.treatment .== "v1_legacy_closed") .&
    (gate_scores.observation_treatment .== "raw"), :]
nrow(scale_rows) == length(TARGETS) || error("Missing frozen CPUE scales")
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
ENV["COTS_SETTLEMENT_PROBABILITY"] = "1.0"
ENV["COTS_INITIAL_MULTIPLIER"] = string(candidate.seed_mult)
ENV["COTS_SEED_FIRST_N"] = "10"
ENV["ADRIA_DEBUG_SEED_FIRST_N"] = "10"

domain = ADRIA.load_domain(ADRIA.LizardDomain, V2, "historical")
forcing = domain.cots_forcing
mapping = load_reef_area_weights(domain, V2, [t.model for t in TARGETS])
legacy_forcing = ADRIA.load_cots_connectivity_forcing(
    joinpath(V2, "cots_connectivity", "annual_forcing"), String.(domain.loc_ids))
legacy_mean_survival = mean(legacy_forcing.source_survival)
owen_mean_survival = mean(forcing.source_survival)
control_fecundity = Float64(candidate.a_ricker) / legacy_mean_survival
production(adults) = control_fecundity * adults *
    exp(-Float64(candidate.b_ricker) * adults) * adults^2 / (3.0^2 + adults^2)
control_boundary = maximum(production.(range(0.0, 100.0; length=10001)))

# Match only the four-reef, 2008-18 mean external pelagic flux of the archived
# control. No per-reef or per-year fitting; this scalar is held across treatments.
potential_sum = 0.0
pelagic_sum = 0.0
for calendar_year in 2008:2018
    tstep = calendar_year - 1985 + 1
    hydro = ADRIA.cots_forcing_index(forcing, tstep - 1, "sample", SEED)
    calendar_index = findfirst(==(calendar_year), forcing.outbreak_years)
    isnothing(calendar_index) && error("Missing boundary year $calendar_year")
    for target in TARGETS
        m = mapping[target.model]
        reefs = forcing.site_to_reef[m.indices]
        global potential_sum += sum(m.weights .* forcing.external_outbreak_potential[
            reefs, calendar_index, hydro])
        global pelagic_sum += sum(m.weights .* forcing.external_outbreak_pelagic[
            reefs, calendar_index, hydro])
    end
end
potential_sum > 0.0 && pelagic_sum > 0.0 || error("Invalid external normalization")
boundary_effective_ratio = pelagic_sum / potential_sum
matched_boundary = control_boundary * boundary_effective_ratio
matched_fecundity = control_fecundity * owen_mean_survival

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

const TREATMENTS = [
    (name="embedded_candidate", fecundity=Float64(candidate.a_ricker), boundary=matched_boundary),
    (name="embedded_mean_matched", fecundity=matched_fecundity, boundary=matched_boundary),
    (name="embedded_local_2x", fecundity=2matched_fecundity, boundary=matched_boundary),
    (name="embedded_boundary_2x", fecundity=matched_fecundity, boundary=2matched_boundary),
    (name="embedded_joint_2x", fecundity=2matched_fecundity, boundary=2matched_boundary),
    (name="embedded_joint_4x", fecundity=4matched_fecundity, boundary=4matched_boundary),
]

raw_obs = CSV.read(OBS, DataFrame)
observed = Dict{String,NamedTuple}()
for target in TARGETS
    rows = sort(raw_obs[(raw_obs.reef_name .== target.observed) .&
        (raw_obs.year .>= 1985) .& (raw_obs.year .<= 2024), :], :year)
    years = Int.(rows.year)
    values = Float64.(rows.cotsptow)
    observed[target.model] = (years=years, raw=values,
        smooth=calendar_moving_average(years, values; window_years=3))
end

scores = DataFrame(treatment=String[], reef_name=String[], observation_treatment=String[],
    loss=Float64[], sim_peak_count=Int[], sim_peak_years=String[],
    sim_peak_heights_cpue=String[], sim_peak_prominences_cpue=String[],
    matched_peaks=Int[], amplitude_penalty=Float64[], peak_timing_penalty=Float64[],
    observation_scale=Float64[])
trajectories = DataFrame(treatment=String[], reef_name=String[], year=Int[],
    forcing_year=Int[], recruits_ha=Float64[], juveniles_ha=Float64[],
    adults_ha=Float64[], simulated_cpue=Float64[], coral_cover=Float64[])
flows = DataFrame(treatment=String[], reef_name=String[], year=Int[])
for name in FLOW_NAMES
    flows[!, name] = Float64[]
end
failures = DataFrame(treatment=String[], error=String[])

for treatment in TREATMENTS
    ENV["COTS_APPLY_LARVAL_SURVIVAL"] = "false"
    ENV["COTS_LARVAL_FECUNDITY"] = string(treatment.fecundity)
    ENV["COTS_EXTERNAL_OUTBREAK_PRODUCTION"] = string(treatment.boundary)
    Random.seed!(SEED)
    try
        result = ADRIA.run_scenario(domain, scenario[1, :])
        years = Int.(result.cots_calendar_year)
        years == collect(1985:2024) || error("Unexpected model years")
        cover = dropdims(sum(result.raw; dims=(2, 3)); dims=(2, 3))
        for target in TARGETS
            reef = target.model
            m = mapping[reef]
            stages = [aggregate_reef_density(Matrix(result.cots_log[:, channel, :]), m)
                for channel in 1:3]
            reef_cover = aggregate_reef_density(cover, m)
            reef_flows = [aggregate_reef_density(Matrix(result.cots_flow_log[:, channel, :]), m)
                for channel in 1:length(FLOW_NAMES)]
            scale = fixed_scales[reef]
            for (index, year) in enumerate(years)
                push!(trajectories, (treatment.name, reef, year,
                    Int(result.cots_forcing_year[index]), stages[1][index],
                    stages[2][index], stages[3][index], stages[3][index] * scale,
                    reef_cover[index]))
                push!(flows, (treatment.name, reef, year,
                    (values[index] for values in reef_flows)...))
            end
            obs = observed[reef]
            for (obs_name, obs_values) in (("raw", obs.raw), ("smoothed_3y", obs.smooth))
                score = cots_cycle_score(years, stages[3], obs.years, obs_values;
                    fixed_observation_scale=scale)
                scoring_window = (years .>= score.scoring_start_year) .&
                    (years .<= score.scoring_end_year)
                peaks = detect_cots_peaks(years[scoring_window],
                    (stages[3] .* scale)[scoring_window]; smooth_window_years=3)
                length(peaks) == score.n_sim_peaks || error("Peak count mismatch")
                push!(scores, (treatment.name, reef, obs_name, score.total_loss,
                    score.n_sim_peaks, join(score.sim_peak_years, ';'),
                    join(round.([p.value for p in peaks]; digits=5), ';'),
                    join(round.([p.prominence for p in peaks]; digits=5), ';'),
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
    "hypothesis" => "Owen matrix embeds larval mortality; no separate RME multiplier",
    "seed" => SEED,
    "temporal_mode" => "sample",
    "fixed_allee_threshold_cots_ha" => 3.0,
    "source_survival_applied_separately" => false,
    "control_fecundity_pre_pelagic_ha_year" => control_fecundity,
    "owen_local_mean_survival" => owen_mean_survival,
    "matched_fecundity_pre_pelagic_ha_year" => matched_fecundity,
    "control_boundary_pre_pelagic_ha_year" => control_boundary,
    "boundary_effective_ratio_2008_2018_four_reef" => boundary_effective_ratio,
    "matched_boundary_pre_pelagic_ha_year" => matched_boundary,
    "normalization_note" => "Single four-reef 2008-18 mean-flux match to archived control; not an evidence estimate",
    "provisional_owen_policy" => "56 unmatched sources omitted; nonfinite full rows zeroed",
    "coral_observation_status" => "9 m LTMP versus whole-reef model diagnostic only",
    "julia_version" => string(VERSION),
    "adria_revision" => git_revision(ROOT),
    "cotsmod_revision" => git_revision(joinpath(ROOT, "..", "COTSMod.jl")),
    "git_status" => read(`git status --porcelain`, String),
    "treatments" => [Dict(String(k) => v for (k, v) in pairs(t)) for t in TREATMENTS],
    "inputs" => Dict(
        "candidate" => file_sha256(BEST),
        "gate_scores" => file_sha256(joinpath(GATE, "scores.csv")),
        "observations" => file_sha256(OBS),
        "script" => file_sha256(@__FILE__),
        "mapping_script" => file_sha256(joinpath(@__DIR__, "reef_observation_mapping.jl")),
        "score_script" => file_sha256(joinpath(@__DIR__, "cots_cycle_metrics.jl")),
        "owen_forcing" => file_sha256(joinpath(V2, "cots_connectivity", "datasets",
            "owen_global_2018_2023", "provenance.json")),
    ),
    "outputs" => Dict(file => file_sha256(joinpath(OUT, file)) for file in
        ("scores.csv", "trajectories.csv", "flows.csv", "failures.csv")),
)
open(joinpath(OUT, "metadata.toml"), "w") do io
    TOML.print(io, metadata; sorted=true)
end
println("Saved ", OUT)
