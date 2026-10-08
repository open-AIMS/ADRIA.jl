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
EXPERIMENT in ("initialization", "cohort_history", "gap_grazing") ||
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
) : Dict{String,String}()
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
initial_reef = DataFrame(treatment=String[], reef_name=String[],
    seed_adults_ha=Float64[], seed_corals_fraction=Float64[],
    idw_cpue=Float64[], idw_manta_coral_fraction=Float64[])
failures = DataFrame(treatment=String[], error=String[])

for t in treatments
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
            scale = fixed_scales[reef]
            for (idx, year) in enumerate(years)
                push!(trajectories, (t.name, reef, year,
                    Int(result.cots_forcing_year[idx]), adults[idx], adults[idx] * scale,
                    coral_reef[idx]))
                push!(flows, (t.name, reef, year,
                    (values[idx] for values in reef_flows)...))
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
end
delete!(ENV, "COTS_INITIAL_STATE_CSV")
delete!(ENV, "ADRIA_COTS_ENABLED")
delete!(ENV, "COTS_DIAGNOSTIC_GRAZING_OFF_START_YEAR")
delete!(ENV, "COTS_DIAGNOSTIC_GRAZING_OFF_END_YEAR")

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
        "inherited 6:3:1 control, with a 2004-2012 experimental grazing bypass",
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
    "adria_revision" => git_revision(ROOT),
    "cotsmod_revision" => git_revision(joinpath(ROOT, "..", "COTSMod.jl")),
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
        "adapter" => file_sha256(joinpath(ROOT, "ADRIA", "src", "ecosystem", "cots.jl")),
        "scenario" => file_sha256(joinpath(ROOT, "ADRIA", "src", "scenario.jl")),
    ),
    "outputs" => Dict(file => file_sha256(joinpath(OUT, file)) for file in
        vcat(["scores.csv", "trajectories.csv", "flows.csv", "initial_reef.csv",
              "failures.csv"], basename.(collect(values(initial_paths))))),
)
open(joinpath(OUT, "metadata.toml"), "w") do io
    TOML.print(io, metadata)
end
isempty(failures.treatment) || error("1991 screen has failed treatments; see failures.csv")
