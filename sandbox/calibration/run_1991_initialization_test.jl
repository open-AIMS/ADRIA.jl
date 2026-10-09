"""Bounded 1991 initialization screen: COTS conversion x coral initial cover.

The 1985-2024 V2/Owen baseline remains untouched. This experimental run slices
historical environmental forcing to 1991-2024, replaces only opt-in initial
states, and retains the archived ecological parameters and observation score.
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
const EXPERIMENT = get(ENV, "OWEN_1991_EXPERIMENT", "initialization")
EXPERIMENT in ("initialization", "cohort_history", "gap_grazing",
    "crash_factorial", "lagged_hazard", "cohort_audit", "structural_crash") ||
    error("Unknown 1991 experiment")
const RUN_ID = get(ENV, "OWEN_1991_RUN_ID",
    Dates.format(now(), "yyyymmddTHHMMSS") * "_owen_1991_" * EXPERIMENT)
occursin(r"^[A-Za-z0-9_.-]+$", RUN_ID) || error("Invalid run ID")
const OUT = joinpath(@__DIR__, "runs", RUN_ID)
ispath(OUT) && error("Refusing to overwrite $OUT")
const V2 = joinpath(ROOT, "sandbox", "data", "Lizard_Historical_v0.2")
const SURFACE_DIR = get(ENV, "OWEN_1991_SURFACE_DIR",
    joinpath(@__DIR__, "runs", "20261007T199104_1991_initial_surfaces"))
const SCREEN = joinpath(@__DIR__, "runs", "20261001T162000_owen_embedded_mortality_screen")
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

surface_path = joinpath(SURFACE_DIR, "surface.csv")
surface = CSV.read(surface_path, DataFrame)
required = (:site_id, :cots_cpue_1989_1991_idw,
    :coral_manta_9m_fraction_1989_1991_idw, :outbreak_probability_1991)
all(in(propertynames(surface)), required) || error("1991 surface has missing columns")
nrow(surface) == 2895 || error("Unexpected 1991 surface site count")
all(isfinite, Matrix{Float64}(surface[:, collect(required[2:end])])) ||
    error("1991 surface contains missing/nonfinite values")
all(x -> 0.0 <= x <= 1.0, surface.coral_manta_9m_fraction_1989_1991_idw) ||
    error("Invalid coral cover surface")
mkpath(OUT)

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
length(fixed_scales) == length(TARGETS) || error("Missing frozen observation scales")

ENV["ADRIA_COTS_CONNECTIVITY_DATASET"] = "owen_global_2018_2023"
ENV["ADRIA_COTS_CONNECTIVITY_MODE"] = "cots"
ENV["COTS_CONNECTIVITY_TEMPORAL_MODE"] = "sample"
ENV["COTS_CONNECTIVITY_SEED"] = string(SEED)
ENV["COTS_ALLEE_THRESHOLD"] = "3.0"
ENV["COTS_EXTERNAL_PULSE"] = "false"
ENV["COTS_EXTERNAL_SOURCE_DENSITY"] = "0.0"
ENV["COTS_JUVENILE_STORAGE"] = "false"
ENV["COTS_SIZE_WEIGHTED_FECUNDITY"] = "false"
ENV["COTS_LAGGED_ADULT_HAZARD"] = "false"
ENV["COTS_HAZARD_FOOD_COUPLED"] = "false"
ENV["COTS_HAZARD_RHO"] = "0.5"
ENV["COTS_HAZARD_THRESHOLD_HA"] = "1.5"
ENV["COTS_HAZARD_WIDTH_HA"] = "0.5"
ENV["COTS_HAZARD_CORAL_THRESHOLD"] = "0.15"
ENV["COTS_HAZARD_CORAL_WIDTH"] = "0.03"
ENV["COTS_HAZARD_MAX"] = "0.0"
ENV["COTS_HABITAT_MEDIATION"] = "false"
ENV["COTS_SEPARATED_RECRUITMENT"] = "true"
ENV["COTS_APPLY_LARVAL_SURVIVAL"] = "false"
ENV["COTS_SETTLEMENT_PROBABILITY"] = "1.0"
ENV["COTS_IMMIGRATION_SCALAR"] = "1.0"
ENV["COTS_LARVAL_FECUNDITY"] = string(fecundity)
ENV["COTS_EXTERNAL_OUTBREAK_PRODUCTION"] = string(boundary)
ENV["COTS_INITIAL_MULTIPLIER"] = string(Float64(candidate.seed_mult))
ENV["COTS_LOW_COVER_MORTALITY"] = "false"
ENV["COTS_LOW_COVER_STRENGTH"] = "0.0"
delete!(ENV, "COTS_INTERNAL_SETTLEMENT_BLACKOUT_START_YEAR")
delete!(ENV, "COTS_INTERNAL_SETTLEMENT_BLACKOUT_END_YEAR")
delete!(ENV, "COTS_INITIAL_STATE_CSV")

domain = ADRIA.load_domain(ADRIA.LizardDomain, V2, "historical")
String.(surface.site_id) == String.(domain.loc_ids) ||
    error("1991 surface order does not match model site order")
mapping = load_reef_area_weights(domain, V2, [t.model for t in TARGETS])
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

# The input NetCDF has absolute calendar years: do not relabel 1985 forcing.
all_years = Int.(domain.env_layer_md.timeframe)
all_years == collect(1985:2024) || error("Unexpected V2 historical years")
keep = findall(>=(1991), all_years)
years = all_years[keep]
domain.dhw_scens = ADRIA.DataCube(Float64.(Array(domain.dhw_scens)[keep, :, :]);
    timesteps=years, sites=domain.loc_ids, scenarios=1:size(domain.dhw_scens, 3))
domain.wave_scens = ADRIA.ZeroDataCube(; T=Float64, timesteps=years,
    locs=domain.loc_ids, scenarios=[1])
domain.cyclone_mortality_scens = ADRIA.ZeroDataCube(; T=Float64,
    timesteps=years, locs=domain.loc_ids,
    species=ADRIA.functional_group_names(), scenarios=[1])
md = domain.env_layer_md
domain.env_layer_md = ADRIA.EnvLayer(md.dpkg_path, md.loc_data_fn, md.loc_id_col,
    md.cluster_id_col, md.init_coral_cov_fn, md.connectivity_fn, md.DHW_fn,
    md.wave_fn, years)
baseline_cover = Float64.(Array(domain.init_coral_cover))
size(baseline_cover, 2) == nrow(surface) || error("Unexpected initial coral dimensions")
baseline_total = vec(sum(baseline_cover; dims=1))
all(>(0.0), baseline_total) || error("Cannot rescale zero initial coral cover")
interpolated_cover = Float64.(surface.coral_manta_9m_fraction_1989_1991_idw)
new_cover = baseline_cover .* reshape(interpolated_cover ./ baseline_total, 1, :)
all(isapprox.(vec(sum(new_cover; dims=1)), interpolated_cover; atol=1e-12)) ||
    error("Interpolated coral totals not preserved")

function set_initial_cover!(values)
    domain.init_coral_cover = ADRIA.DataCube(values;
        species=collect(domain.init_coral_cover.species), locations=domain.loc_ids)
end

function make_initial_state(alpha::Float64, name::String; cohort_multiplier::Float64=0.0)::String
    # alpha is an explicitly provisional COTS/tow per adult COTS/ha conversion.
    cpue = Float64.(surface.cots_cpue_1989_1991_idw)
    adults = cpue ./ alpha
    state = DataFrame(site_id=String.(surface.site_id),
        juvenile_ha=6.0 .* cohort_multiplier .* adults,
        subadult_ha=3.0 .* cohort_multiplier .* adults,
        adult_ha=adults)
    path = joinpath(OUT, "initial_state_$(name).csv")
    CSV.write(path, state)
    return path
end

treatments = EXPERIMENT == "initialization" ? [
    (name="inherited_1991", alpha=NaN, coral="inherited_2023", initial="legacy"),
    (name="idw_alpha015_old_coral", alpha=0.015, coral="inherited_2023", initial="adult_only"),
    (name="idw_alpha015_idw_coral", alpha=0.015, coral="manta_1989_1991", initial="adult_only"),
    (name="idw_alpha150_old_coral", alpha=0.150, coral="inherited_2023", initial="adult_only"),
    (name="idw_alpha150_idw_coral", alpha=0.150, coral="manta_1989_1991", initial="adult_only"),
] : EXPERIMENT == "cohort_history" ? [
    (name="inherited_1991", alpha=NaN, coral="inherited_2023", initial="legacy"),
    (name="idw_alpha015_cohort_half", alpha=0.015, coral="inherited_2023", initial="cohort_half"),
    (name="idw_alpha015_cohort_full", alpha=0.015, coral="inherited_2023", initial="cohort_full"),
    (name="idw_alpha150_cohort_full", alpha=0.150, coral="inherited_2023", initial="cohort_full"),
    (name="coral_only_no_cots", alpha=NaN, coral="inherited_2023", initial="no_cots"),
] : EXPERIMENT == "crash_factorial" ? [
    (name="cohort_full_control", alpha=0.015, coral="inherited_2023", initial="cohort_full"),
    (name="cohort_full_stage_gate", alpha=0.015, coral="inherited_2023", initial="cohort_full"),
    (name="cohort_full_mortality030", alpha=0.015, coral="inherited_2023", initial="cohort_full"),
    (name="cohort_full_stage_gate_mortality030", alpha=0.015, coral="inherited_2023", initial="cohort_full"),
] : EXPERIMENT == "lagged_hazard" ? [
    (name="cohort_full_control", alpha=0.015, coral="inherited_2023", initial="cohort_full"),
    (name="density_h075", alpha=0.015, coral="inherited_2023", initial="cohort_full"),
    (name="density_h150", alpha=0.015, coral="inherited_2023", initial="cohort_full"),
    (name="density_food_h075", alpha=0.015, coral="inherited_2023", initial="cohort_full"),
    (name="density_food_h150", alpha=0.015, coral="inherited_2023", initial="cohort_full"),
] : EXPERIMENT == "cohort_audit" ? [
    (name="cohort_full_control", alpha=0.015, coral="inherited_2023", initial="cohort_full"),
] : EXPERIMENT == "structural_crash" ? [
    (name="cohort_full_control", alpha=0.015, coral="inherited_2023", initial="cohort_full"),
    (name="preferred_starvation", alpha=0.015, coral="inherited_2023", initial="cohort_full"),
    (name="senescence", alpha=0.015, coral="inherited_2023", initial="cohort_full"),
    (name="preferred_starvation_senescence", alpha=0.015, coral="inherited_2023", initial="cohort_full"),
] : [
    (name="inherited_1991", alpha=NaN, coral="inherited_2023", initial="legacy"),
    (name="inherited_grazing_off_2004_2012", alpha=NaN,
        coral="inherited_2023", initial="legacy"),
]
initial_paths = EXPERIMENT == "initialization" ? Dict(
    "idw_alpha015" => make_initial_state(0.015, "alpha015"),
    "idw_alpha150" => make_initial_state(0.150, "alpha150"),
) : EXPERIMENT == "cohort_history" ? Dict(
    "idw_alpha015_cohort_half" => make_initial_state(0.015, "alpha015_cohort_half"; cohort_multiplier=0.5),
    "idw_alpha015_cohort_full" => make_initial_state(0.015, "alpha015_cohort_full"; cohort_multiplier=1.0),
    "idw_alpha150_cohort_full" => make_initial_state(0.150, "alpha150_cohort_full"; cohort_multiplier=1.0),
) : EXPERIMENT in ("crash_factorial", "lagged_hazard", "cohort_audit", "structural_crash") ? let
    state_path = make_initial_state(0.015, "alpha015_cohort_full"; cohort_multiplier=1.0)
    Dict(t.name => state_path for t in treatments)
end : Dict{String,String}()
archived_m3 = Float64(scenario[1, :m3])
raw_obs = CSV.read(OBS, DataFrame)
scores = DataFrame(treatment=String[], reef_name=String[], observation_treatment=String[],
    loss=Float64[], sim_peak_count=Int[], sim_peak_years=String[],
    sim_peak_heights_cpue=String[], matched_peaks=Int[])
trajectories = DataFrame(treatment=String[], reef_name=String[], year=Int[],
    forcing_year=Int[], adults_ha=Float64[], simulated_cpue=Float64[], coral_cover=Float64[])
flows = DataFrame(treatment=String[], reef_name=String[], year=Int[])
for name in FLOW_NAMES
    flows[!, name] = Float64[]
end
hazards = DataFrame(treatment=String[], reef_name=String[], year=Int[],
    burden_ha=Float64[], hazard=Float64[], adult_deaths_ha=Float64[])
cohorts = DataFrame(treatment=String[], reef_name=String[], year=Int[],
    recruits_ha=Float64[], subadults_ha=Float64[], adults_ha=Float64[],
    opening_recruits_ha=Float64[], opening_subadults_ha=Float64[],
    opening_adults_ha=Float64[], maturation_ha=Float64[],
    adult_carryover_ha=Float64[], adult_losses_ha=Float64[],
    subadult_food_factor=Float64[], adult_food_factor=Float64[],
    recruits_to_subadults_ha=Float64[])
# Preregistered structural screen (crash_mechanism_protocol.md, 2026-10-09).
const STRUCTURAL_THRESHOLD = 0.05     # preferred cover, fraction of habitable area
const STRUCTURAL_MEMORY_YEARS = 0.5
const STRUCTURAL_SENESCENCE_AGE = 6   # years
const STRUCTURAL_SENESCENCE_MORTALITY = 0.8  # annual fraction
mechanisms = DataFrame(treatment=String[], reef_name=String[], year=Int[],
    recruits_ha=Float64[], subadults_ha=Float64[], adults_ha=Float64[],
    maturation_ha=Float64[], fast_coral=Float64[], total_coral=Float64[],
    preferred_food_memory=Float64[], food_survival=Float64[],
    senescent_adults_ha=Float64[], senescent_deaths_ha=Float64[])
initial_reef = DataFrame(treatment=String[], reef_name=String[],
    seed_adults_ha=Float64[], seed_corals_fraction=Float64[],
    idw_cpue=Float64[], idw_manta_coral_fraction=Float64[])
failures = DataFrame(treatment=String[], error=String[])

for t in treatments
    ENV["COTS_LAGGED_ADULT_HAZARD"] = EXPERIMENT == "lagged_hazard" &&
        t.name != "cohort_full_control" ? "true" : "false"
    ENV["COTS_HAZARD_FOOD_COUPLED"] = EXPERIMENT == "lagged_hazard" &&
        occursin("density_food", t.name) ? "true" : "false"
    ENV["COTS_HAZARD_MAX"] = EXPERIMENT == "lagged_hazard" &&
        occursin("h150", t.name) ? "1.5" :
        EXPERIMENT == "lagged_hazard" && occursin("h075", t.name) ? "0.75" : "0.0"
    ENV["COTS_JUVENILE_STORAGE"] = EXPERIMENT == "crash_factorial" &&
        occursin("stage_gate", t.name) ? "true" : "false"
    structural = EXPERIMENT == "structural_crash"
    ENV["COTS_PREFERRED_PREY_STARVATION"] = structural &&
        occursin("preferred_starvation", t.name) ? "true" : "false"
    ENV["COTS_STARVATION_PREFERRED_THRESHOLD"] = string(STRUCTURAL_THRESHOLD)
    ENV["COTS_STARVATION_MEMORY_YEARS"] = string(STRUCTURAL_MEMORY_YEARS)
    ENV["COTS_ADULT_SENESCENCE"] = structural &&
        occursin("senescence", t.name) ? "true" : "false"
    ENV["COTS_SENESCENCE_AGE"] = string(STRUCTURAL_SENESCENCE_AGE)
    ENV["COTS_SENESCENCE_MORTALITY"] = string(STRUCTURAL_SENESCENCE_MORTALITY)
    ENV["COTS_PER_CAPITA_CONSUMPTION"] = "false"
    scenario[!, :m3] .= EXPERIMENT == "crash_factorial" &&
        occursin("mortality030", t.name) ? 0.30 : archived_m3
    set_initial_cover!(t.coral == "manta_1989_1991" ? new_cover : baseline_cover)
    if t.name == "inherited_grazing_off_2004_2012"
        ENV["COTS_DIAGNOSTIC_GRAZING_OFF_START_YEAR"] = "2004"
        ENV["COTS_DIAGNOSTIC_GRAZING_OFF_END_YEAR"] = "2012"
    else
        delete!(ENV, "COTS_DIAGNOSTIC_GRAZING_OFF_START_YEAR")
        delete!(ENV, "COTS_DIAGNOSTIC_GRAZING_OFF_END_YEAR")
    end
    ENV["ADRIA_COTS_ENABLED"] = t.initial == "no_cots" ? "false" : "true"
    ENV["COTS_LARVAL_FECUNDITY"] = t.initial == "no_cots" ? "0.0" : string(fecundity)
    ENV["COTS_EXTERNAL_OUTBREAK_PRODUCTION"] = t.initial == "no_cots" ? "0.0" : string(boundary)
    ENV["COTS_IMMIGRATION_SCALAR"] = t.initial == "no_cots" ? "0.0" : "1.0"
    if t.initial in ("legacy", "no_cots")
        delete!(ENV, "COTS_INITIAL_STATE_CSV")
    else
        key = EXPERIMENT == "initialization" ?
            (t.alpha == 0.015 ? "idw_alpha015" : "idw_alpha150") : t.name
        ENV["COTS_INITIAL_STATE_CSV"] = initial_paths[key]
    end
    Random.seed!(SEED)
    try
        result = ADRIA.run_scenario(domain, scenario[1, :])
        Int.(result.cots_calendar_year) == years || error("1991 calendar alignment failed")
        if t.initial == "no_cots"
            maximum(abs, result.cots_log) == 0.0 ||
                error("No-COTS diagnostic received nonzero COTS state")
        end
        coral = dropdims(sum(result.raw; dims=(2, 3)); dims=(2, 3))
        fast_coral = dropdims(sum(Array(result.raw)[:, 1:3, :, :]; dims=(2, 3)); dims=(2, 3))
        if EXPERIMENT == "lagged_hazard"
            size(result.cots_hazard_log) == (length(years), 3, length(domain.loc_ids)) ||
                error("Unexpected site hazard log dimensions")
            site_table = DataFrame(
                treatment=fill(t.name, length(years) * length(domain.loc_ids)),
                site_index=repeat(collect(1:length(domain.loc_ids)); inner=length(years)),
                year=repeat(years, length(domain.loc_ids)),
                adults_ha=vec(Matrix(result.cots_log[:, 3, :])),
                coral_cover=vec(coral),
                burden_ha=vec(Matrix(result.cots_hazard_log[:, 1, :])),
                hazard=vec(Matrix(result.cots_hazard_log[:, 2, :])),
                adult_deaths_ha=vec(Matrix(result.cots_hazard_log[:, 3, :]))
            )
            path = joinpath(OUT, "site_hazard.csv")
            CSV.write(path, site_table; append=isfile(path), writeheader=!isfile(path))
        end
        for target in TARGETS
            reef = target.model
            m = mapping[reef]
            adults = aggregate_reef_density(Matrix(result.cots_log[:, 3, :]), m)
            coral_reef = aggregate_reef_density(coral, m)
            state_seed = t.initial == "no_cots" ? zeros(length(domain.loc_ids)) :
                t.initial == "legacy" ?
                0.1 .* [p > 0.1 ? min(p / 0.5, 1.0) : 0.0
                    for p in domain.cots_init_density] .* Float64(candidate.seed_mult) :
                Float64.(surface.cots_cpue_1989_1991_idw) ./ t.alpha
            push!(initial_reef, (t.name, reef,
                sum(m.weights .* state_seed[m.indices]),
                sum(m.weights .* vec(sum(Array(domain.init_coral_cover); dims=1))[m.indices]),
                sum(m.weights .* Float64.(surface.cots_cpue_1989_1991_idw[m.indices])),
                sum(m.weights .* interpolated_cover[m.indices])))
            reef_flows = [aggregate_reef_density(Matrix(result.cots_flow_log[:, channel, :]), m)
                for channel in 1:length(FLOW_NAMES)]
            reef_stages = EXPERIMENT in ("cohort_audit", "structural_crash") ?
                [aggregate_reef_density(Matrix(result.cots_log[:, stage, :]), m)
                    for stage in 1:3] : Vector{Float64}[]
            reef_hazards = EXPERIMENT == "lagged_hazard" ?
                [aggregate_reef_density(Matrix(result.cots_hazard_log[:, channel, :]), m)
                    for channel in 1:3] : Vector{Float64}[]
            reef_mechanisms = EXPERIMENT == "structural_crash" ?
                [aggregate_reef_density(Matrix(result.cots_mechanism_log[:, channel, :]), m)
                    for channel in 1:4] : Vector{Float64}[]
            reef_fast_coral = aggregate_reef_density(fast_coral, m)
            scale = fixed_scales[reef]
            for (idx, year) in enumerate(years)
                if EXPERIMENT == "structural_crash"
                    push!(mechanisms, (t.name, reef, year,
                        reef_stages[1][idx], reef_stages[2][idx], reef_stages[3][idx],
                        reef_flows[5][idx], reef_fast_coral[idx], coral_reef[idx],
                        (values[idx] for values in reef_mechanisms)...))
                end
                push!(trajectories, (t.name, reef, year,
                    Int(result.cots_forcing_year[idx]), adults[idx], adults[idx] * scale,
                    coral_reef[idx]))
                push!(flows, (t.name, reef, year,
                    (values[idx] for values in reef_flows)...))
                if EXPERIMENT == "lagged_hazard"
                    push!(hazards, (t.name, reef, year,
                        (values[idx] for values in reef_hazards)...))
                end
                if EXPERIMENT == "cohort_audit" && idx > 1
                    initial_adults = Float64.(surface.cots_cpue_1989_1991_idw) ./ t.alpha
                    opening = idx == 2 ?
                        (6.0 * sum(m.weights .* initial_adults[m.indices]),
                         3.0 * sum(m.weights .* initial_adults[m.indices]),
                         sum(m.weights .* initial_adults[m.indices])) :
                        (reef_stages[1][idx - 1], reef_stages[2][idx - 1],
                         reef_stages[3][idx - 1])
                    maturation = reef_flows[5][idx]
                    carryover = reef_stages[3][idx] - maturation
                    carryover >= -1e-9 || error("Negative adult carryover")
                    abs(reef_stages[2][idx] - opening[1] * (1.0 - scenario[1, :m1])) < 1e-8 ||
                        error("Unexpected subadult transition with storage off")
                    subadult_denominator = opening[2] * (1.0 - scenario[1, :m2])
                    adult_denominator = opening[3] * (1.0 - scenario[1, :m3])
                    push!(cohorts, (t.name, reef, year,
                        reef_stages[1][idx], reef_stages[2][idx], reef_stages[3][idx],
                        opening[1], opening[2], opening[3], maturation, carryover,
                        opening[3] - carryover,
                        subadult_denominator > 0 ? maturation / subadult_denominator : NaN,
                        adult_denominator > 0 ? carryover / adult_denominator : NaN,
                        reef_stages[2][idx]))
                end
            end
            t.initial == "no_cots" && continue
            obs = sort(raw_obs[(raw_obs.reef_name .== target.observed) .&
                (raw_obs.year .>= 1992) .& (raw_obs.year .<= 2024), :], :year)
            obs_years = Int.(obs.year)
            obs_values = Float64.(obs.cotsptow)
            smoothed = calendar_moving_average(obs_years, obs_values; window_years=3)
            for (name, values) in (("raw", obs_values), ("smoothed_3y", smoothed))
                score = cots_cycle_score(years[2:end], adults[2:end], obs_years, values;
                    fixed_observation_scale=scale)
                window = (years .>= score.scoring_start_year) .&
                    (years .<= score.scoring_end_year)
                peaks = detect_cots_peaks(years[window], (adults .* scale)[window];
                    smooth_window_years=3)
                push!(scores, (t.name, reef, name, score.total_loss,
                    score.n_sim_peaks, join(score.sim_peak_years, ';'),
                    join(round.([p.value for p in peaks]; digits=5), ';'),
                    score.n_matched_peaks))
            end
        end
        println("Completed ", t.name)
    catch err
        push!(failures, (t.name, sprint(showerror, err, catch_backtrace())))
        println(stderr, "FAILED ", t.name, ": ", err)
    end
    flush(stdout)
    for (file, table) in (("scores.csv", scores), ("trajectories.csv", trajectories),
        ("flows.csv", flows), ("initial_reef.csv", initial_reef),
        ("failures.csv", failures))
        CSV.write(joinpath(OUT, file), table)
    end
    EXPERIMENT == "lagged_hazard" && CSV.write(joinpath(OUT, "hazards.csv"), hazards)
    EXPERIMENT == "cohort_audit" && CSV.write(joinpath(OUT, "cohorts.csv"), cohorts)
    EXPERIMENT == "structural_crash" && CSV.write(joinpath(OUT, "mechanisms.csv"), mechanisms)
end
delete!(ENV, "COTS_INITIAL_STATE_CSV")
delete!(ENV, "ADRIA_COTS_ENABLED")
delete!(ENV, "COTS_DIAGNOSTIC_GRAZING_OFF_START_YEAR")
delete!(ENV, "COTS_DIAGNOSTIC_GRAZING_OFF_END_YEAR")
delete!(ENV, "COTS_LAGGED_ADULT_HAZARD")
delete!(ENV, "COTS_HAZARD_FOOD_COUPLED")
delete!(ENV, "COTS_HAZARD_MAX")
for name in ("COTS_PREFERRED_PREY_STARVATION", "COTS_STARVATION_PREFERRED_THRESHOLD",
        "COTS_STARVATION_MEMORY_YEARS", "COTS_ADULT_SENESCENCE", "COTS_SENESCENCE_AGE",
        "COTS_SENESCENCE_MORTALITY", "COTS_PER_CAPITA_CONSUMPTION")
    delete!(ENV, name)
end

metadata = Dict(
    "status" => "bounded_1991_$(EXPERIMENT)_screen_not_promoted",
    "experiment" => EXPERIMENT,
    "grazing_off_years" => EXPERIMENT == "gap_grazing" ? [2004, 2012] : Int[],
    "year_window" => [1991, 2024], "score_start_year" => 1992,
    "seed" => SEED, "allee_threshold_adults_ha" => 3.0,
    "cots_initial_profile" => EXPERIMENT == "initialization" ?
        "adult-only for IDW treatments; inherited 6:3:1 for control" :
        EXPERIMENT == "cohort_history" ?
        "IDW adult CPUE with 0.5 or 1.0 times inherited 6:3:1 immature:adult ratio; no-COTS coral counterfactual" :
        EXPERIMENT in ("crash_factorial", "lagged_hazard", "cohort_audit", "structural_crash") ?
        "all arms: IDW adult CPUE at alpha=0.015 and assumed 6:3:1 immature:adult ratio" :
        "inherited 6:3:1 control, with a 2004-2012 experimental grazing bypass",
    "archived_adult_mortality_fraction_per_year" => archived_m3,
    "crash_factorial_adult_mortality_fraction_per_year" => EXPERIMENT == "crash_factorial" ? 0.30 : archived_m3,
    "crash_factorial_stage_gate" => EXPERIMENT == "crash_factorial" ?
        "existing juvenile_storage switch; min=0.05 max=1 cover=0.4" : "not tested",
    "lagged_hazard_factorial" => EXPERIMENT == "lagged_hazard" ?
        "control plus density-only and density-times-food at hazard_max=0.75 and 1.5; rho=0.5 threshold=1.5 COTS/ha width=0.5 COTS/ha coral threshold=0.15 fraction coral width=0.03 fraction" : "not tested",
    "structural_crash_factorial" => EXPERIMENT == "structural_crash" ?
        "control, preferred_starvation, senescence and both. preferred threshold=$(STRUCTURAL_THRESHOLD) fast-cover fraction of habitable area; memory=$(STRUCTURAL_MEMORY_YEARS) yr; senescence age=$(STRUCTURAL_SENESCENCE_AGE) yr; senescent mortality=$(STRUCTURAL_SENESCENCE_MORTALITY) annual fraction; ages 2-5 keep archived m3; initial adults on stable age profile; per-capita consumption off" : "not tested",
    "alpha_015_meaning" => "0.015 COTS/tow per adult/ha; provisional threshold-equivalence sensitivity, not validated conversion",
    "alpha_150_meaning" => "0.150 COTS/tow per adult/ha; provisional high conversion sensitivity, not validated conversion",
    "coral_initialization" => "9 m manta fraction used as whole-reef proxy; preserve inherited species/size proportions",
    "outbreak_probability" => "1991 values recorded on initial surface; existing Owen external boundary uses them, not initial CPUE",
    "source_to_sink_orientation" => "Owen source rows, sink columns",
    "external_production_pre_pelagic_recruits_ha_year" => boundary,
    "fecundity_pre_pelagic_recruits_adult_year" => fecundity,
    "hydrodynamic_mode" => "sample", "hydrodynamic_seed" => SEED,
    "julia_version" => string(VERSION),
    "active_project" => string(Base.active_project()),
    "adria_revision" => get(ENV, "COTS_RUN_ADRIA_REVISION", git_revision(ROOT)),
    "cotsmod_revision" => get(ENV, "COTS_RUN_COTSMOD_REVISION",
        git_revision(joinpath(ROOT, "..", "COTSMod.jl"))),
    "inputs" => Dict(
        "surface" => file_sha256(surface_path),
        "surface_provenance" => file_sha256(joinpath(SURFACE_DIR, "metadata.json")),
        "observations" => file_sha256(OBS),
        "candidate" => file_sha256(BEST),
        "frozen_observation_scale" => file_sha256(joinpath(GATE, "scores.csv")),
        "historical_dhw" => file_sha256(joinpath(V2, "DHWs", "dhw_RCPhistorical.nc")),
        "initial_coral" => file_sha256(joinpath(V2, "initial_cover.nc")),
        "owen_forcing" => file_sha256(joinpath(V2, "cots_connectivity", "datasets",
            "owen_global_2018_2023", "provenance.json")),
        "script" => file_sha256(@__FILE__),
        "archived_cohort_trajectories" => EXPERIMENT in ("crash_factorial", "lagged_hazard", "cohort_audit", "structural_crash") ?
            file_sha256(joinpath(@__DIR__, "runs", "20261007T199111_owen_1991_cohort_history_final", "trajectories.csv")) : "not_applicable",
        "archived_cohort_initial_state" => EXPERIMENT in ("crash_factorial", "lagged_hazard", "cohort_audit", "structural_crash") ?
            file_sha256(joinpath(@__DIR__, "runs", "20261007T199111_owen_1991_cohort_history_final", "initial_state_alpha015_cohort_full.csv")) : "not_applicable",
        "adapter" => file_sha256(joinpath(ROOT, "ADRIA", "src", "ecosystem", "cots.jl")),
        "scenario" => file_sha256(joinpath(ROOT, "ADRIA", "src", "scenario.jl")),
        "core" => file_sha256(joinpath(ROOT, "..", "COTSMod.jl", "src", "COTSMod.jl")),
        "score_definition" => file_sha256(joinpath(@__DIR__, "cots_cycle_metrics.jl")),
        "reef_mapping" => file_sha256(joinpath(@__DIR__, "reef_observation_mapping.jl")),
        "project" => file_sha256(joinpath(ROOT, "sandbox", "Project.toml")),
        "manifest" => file_sha256(joinpath(ROOT, "sandbox", "Manifest.toml")),
    ),
    "outputs" => Dict(file => file_sha256(joinpath(OUT, file)) for file in
        vcat(["scores.csv", "trajectories.csv", "flows.csv", "initial_reef.csv",
              "failures.csv"], EXPERIMENT == "lagged_hazard" ?
              ["hazards.csv", "site_hazard.csv"] :
              EXPERIMENT == "cohort_audit" ? ["cohorts.csv"] :
              EXPERIMENT == "structural_crash" ? ["mechanisms.csv"] : String[],
              basename.(collect(values(initial_paths))))),
)
open(joinpath(OUT, "metadata.toml"), "w") do io
    TOML.print(io, metadata)
end
isempty(failures.treatment) || error("1991 screen has failed treatments; see failures.csv")
