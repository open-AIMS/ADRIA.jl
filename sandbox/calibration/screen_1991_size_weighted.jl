"""Bounded, opt-in 1991 Owen test of size-weighted adult COTS production.

The only ecological changes among arms are breeder weighting and small-to-large
growth. The Allee half-saturation is always 3 adult COTS/ha, with the same
provisional global tow conversion, initial stage state, coral and forcing seed.
"""

using Pkg
const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
Pkg.activate(joinpath(ROOT, "sandbox"))
cd(ROOT)

using ADRIA, CSV, DataFrames, Dates, Random, TOML, Statistics
include(joinpath(@__DIR__, "cots_cycle_metrics.jl"))
include(joinpath(@__DIR__, "reef_observation_mapping.jl"))
include(joinpath(@__DIR__, "calibration_run.jl"))

const SEED = 20260930
const RUN_ID = get(ENV, "OWEN_SIZE_RUN_ID",
    Dates.format(now(), "yyyymmddTHHMMSS") * "_owen_1991_size_weighted")
occursin(r"^[A-Za-z0-9_.-]+$", RUN_ID) || error("Invalid run ID")
const OUT = joinpath(@__DIR__, "runs", RUN_ID)
ispath(OUT) && error("Refusing to overwrite $OUT")
const V2 = joinpath(ROOT, "sandbox", "data", "Lizard_Historical_v0.2")
const CAPACITY = joinpath(@__DIR__, "runs", "20261008Tcycle_capacity_24pairs")
const SCREEN = joinpath(@__DIR__, "runs", "20261001T162000_owen_embedded_mortality_screen")
const BEST = joinpath(@__DIR__, "runs", "expanded_cotsconn_pilot_seed20260930", "best_summary.csv")
const OBS = joinpath(ROOT, "sandbox", "data", "reef_cots.csv")
const PAIR = 20
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
const SIZE_NAMES = ["small_adults", "large_adults", "effective_breeders",
    "small_adult_growth"]

plan = DataFrame(treatment=["control", "size_fixed_growth", "size_food_growth",
        "size_food_growth_large75", "size_food_growth_d350"],
    enabled=[false, true, true, true, true],
    food_mediated=[false, false, true, true, true],
    initial_large_fraction=[0.0, 0.25, 0.25, 0.75, 0.25],
    large_diameter_mm=[300.0, 300.0, 300.0, 300.0, 350.0])
mkpath(OUT)
CSV.write(joinpath(OUT, "design.csv"), plan)

pair_rows = CSV.read(joinpath(CAPACITY, "design.csv"), DataFrame)
pair = only(eachrow(pair_rows[pair_rows.pair .== PAIR, :]))
q = Float64(pair.q_tow_per_adult_ha)
initial_path = joinpath(CAPACITY, "initial_state_pair_020.csv")
reference = CSV.read(joinpath(CAPACITY, "trajectories.csv"), DataFrame)
reference = reference[reference.treatment .== "pair_020_theta_3", :]
nrow(reference) == 34length(TARGETS) || error("Incomplete archived pair 20 reference")
screen_meta = TOML.parsefile(joinpath(SCREEN, "metadata.toml"))
fecundity = 4Float64(screen_meta["matched_fecundity_pre_pelagic_ha_year"])
boundary = 4Float64(screen_meta["matched_boundary_pre_pelagic_ha_year"])
candidate = sort(CSV.read(BEST, DataFrame), :loss)[1, :]

ENV["ADRIA_COTS_CONNECTIVITY_DATASET"] = "owen_global_2018_2023"
ENV["ADRIA_COTS_CONNECTIVITY_MODE"] = "cots"
ENV["COTS_CONNECTIVITY_TEMPORAL_MODE"] = "sample"
ENV["COTS_CONNECTIVITY_SEED"] = string(SEED)
ENV["COTS_EXTERNAL_PULSE"] = "false"
ENV["COTS_EXTERNAL_SOURCE_DENSITY"] = "0.0"
ENV["COTS_JUVENILE_STORAGE"] = "false"
ENV["COTS_HABITAT_MEDIATION"] = "false"
ENV["COTS_SEPARATED_RECRUITMENT"] = "true"
ENV["COTS_APPLY_LARVAL_SURVIVAL"] = "false" # Embedded in Owen transport.
ENV["COTS_SETTLEMENT_PROBABILITY"] = "1.0"
ENV["COTS_IMMIGRATION_SCALAR"] = "1.0"
ENV["COTS_INITIAL_MULTIPLIER"] = string(Float64(candidate.seed_mult))
ENV["COTS_LOW_COVER_MORTALITY"] = "false"
ENV["COTS_LOW_COVER_STRENGTH"] = "0.0"
ENV["COTS_INITIAL_STATE_CSV"] = initial_path
ENV["COTS_LARVAL_FECUNDITY"] = string(fecundity * pair.fecundity_multiplier)
ENV["COTS_EXTERNAL_OUTBREAK_PRODUCTION"] = string(boundary * pair.boundary_multiplier)
ENV["COTS_ALLEE_THRESHOLD"] = "3.0"
ENV["COTS_SIZE_SMALL_DIAMETER_MM"] = "200.0"
ENV["COTS_SIZE_FECUNDITY_SLOPE_PER_MM"] = "0.0115"
ENV["COTS_SIZE_GROWTH_MAX"] = "0.5"
ENV["COTS_SIZE_GROWTH_COVER"] = "0.4"
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
base_scenario[!, :m3] .= pair.adult_mortality
base_scenario[!, :a_F] .= Float64(candidate.a_F) * pair.feeding_multiplier
base_scenario[!, :a_S] .= Float64(candidate.a_S) * pair.feeding_multiplier
all_years = Int.(domain.env_layer_md.timeframe)
all_years == collect(1985:2024) || error("Unexpected V2 historical years")
keep = findall(>=(1991), all_years)
years = all_years[keep]
domain.dhw_scens = ADRIA.DataCube(Float64.(Array(domain.dhw_scens)[keep, :, :]);
    timesteps=years, sites=domain.loc_ids, scenarios=1:size(domain.dhw_scens, 3))
domain.wave_scens = ADRIA.ZeroDataCube(; T=Float64, timesteps=years,
    locs=domain.loc_ids, scenarios=[1])
domain.cyclone_mortality_scens = ADRIA.ZeroDataCube(; T=Float64, timesteps=years,
    locs=domain.loc_ids, species=ADRIA.functional_group_names(), scenarios=[1])
md = domain.env_layer_md
domain.env_layer_md = ADRIA.EnvLayer(md.dpkg_path, md.loc_data_fn, md.loc_id_col,
    md.cluster_id_col, md.init_coral_cov_fn, md.connectivity_fn, md.DHW_fn,
    md.wave_fn, years)
raw_obs = CSV.read(OBS, DataFrame)

scores = DataFrame(treatment=String[], reef_name=String[], loss=Float64[],
    sim_peak_count=Int[], sim_peak_years=String[], sim_peak_heights_cpue=String[],
    matched_peaks=Int[])
trajectories = DataFrame(treatment=String[], reef_name=String[], year=Int[],
    forcing_year=Int[], adults_ha=Float64[], simulated_cpue=Float64[],
    coral_cover=Float64[])
flows = DataFrame(treatment=String[], reef_name=String[], year=Int[])
for name in vcat(FLOW_NAMES, SIZE_NAMES)
    flows[!, name] = Float64[]
end
failures = DataFrame(treatment=String[], error=String[])

for arm in eachrow(plan)
    ENV["COTS_SIZE_WEIGHTED_FECUNDITY"] = string(arm.enabled)
    ENV["COTS_SIZE_GROWTH_FOOD_MEDIATED"] = string(arm.food_mediated)
    ENV["COTS_SIZE_INITIAL_LARGE_FRACTION"] = string(arm.initial_large_fraction)
    ENV["COTS_SIZE_LARGE_DIAMETER_MM"] = string(arm.large_diameter_mm)
    Random.seed!(SEED)
    try
        result = ADRIA.run_scenario(domain, base_scenario[1, :])
        Int.(result.cots_calendar_year) == years || error("Calendar alignment failed")
        coral = dropdims(sum(result.raw; dims=(2, 3)); dims=(2, 3))
        for target in TARGETS
            reef = target.model
            m = mapping[reef]
            adult = aggregate_reef_density(Matrix(result.cots_log[:, 3, :]), m)
            reef_coral = aggregate_reef_density(coral, m)
            reef_flows = [aggregate_reef_density(Matrix(result.cots_flow_log[:, channel, :]), m)
                for channel in 1:length(FLOW_NAMES)]
            reef_sizes = [aggregate_reef_density(Matrix(result.cots_size_log[:, channel, :]), m)
                for channel in 1:length(SIZE_NAMES)]
            all(isapprox.(reef_sizes[1] .+ reef_sizes[2], adult; atol=1e-9, rtol=1e-9)) ||
                error("Small/large adult accounting differs from total adults")
            if arm.treatment == "control"
                old = sort(reference[reference.reef_name .== reef, :], :year)
                years == Int.(old.year) || error("Reference year mismatch")
                Int.(result.cots_forcing_year) == Int.(old.forcing_year) ||
                    error("Reference forcing-year mismatch")
                all(isapprox.(adult, Float64.(old.adults_ha); atol=1e-12, rtol=0)) ||
                    error("Legacy adult baseline changed")
                all(isapprox.(reef_coral, Float64.(old.coral_cover); atol=1e-12, rtol=0)) ||
                    error("Legacy coral baseline changed")
            end
            for (idx, year) in enumerate(years)
                push!(trajectories, (arm.treatment, reef, year,
                    Int(result.cots_forcing_year[idx]), adult[idx], adult[idx] * q,
                    reef_coral[idx]))
                push!(flows, (arm.treatment, reef, year,
                    (v[idx] for v in vcat(reef_flows, reef_sizes))...))
            end
            obs = sort(raw_obs[(raw_obs.reef_name .== target.observed) .&
                (raw_obs.year .>= 1992) .& (raw_obs.year .<= 2024), :], :year)
            obs_years = Int.(obs.year)
            smoothed = calendar_moving_average(obs_years, Float64.(obs.cotsptow);
                window_years=3)
            score = cots_cycle_score(years[2:end], adult[2:end], obs_years,
                smoothed; fixed_observation_scale=q)
            window = (years .>= score.scoring_start_year) .&
                (years .<= score.scoring_end_year)
            peaks = detect_cots_peaks(years[window], (adult .* q)[window];
                smooth_window_years=3)
            push!(scores, (arm.treatment, reef, score.total_loss, score.n_sim_peaks,
                join(score.sim_peak_years, ';'),
                join(round.([p.value for p in peaks]; digits=6), ';'),
                score.n_matched_peaks))
        end
        println("Completed ", arm.treatment)
    catch err
        push!(failures, (arm.treatment, sprint(showerror, err, catch_backtrace())))
        println(stderr, "FAILED ", arm.treatment, ": ", err)
    end
    flush(stdout)
    for (name, table) in (("scores.csv", scores), ("trajectories.csv", trajectories),
        ("flows.csv", flows), ("failures.csv", failures))
        CSV.write(joinpath(OUT, name), table)
    end
end

inputs = Dict(
    "capacity_design" => file_sha256(joinpath(CAPACITY, "design.csv")),
    "capacity_initial_state" => file_sha256(initial_path),
    "reference_trajectories" => file_sha256(joinpath(CAPACITY, "trajectories.csv")),
    "reference_parameters" => file_sha256(BEST),
    "screen_metadata" => file_sha256(joinpath(SCREEN, "metadata.toml")),
    "observations" => file_sha256(OBS),
    "score_definition" => file_sha256(joinpath(@__DIR__, "cots_cycle_metrics.jl")),
    "area_weighted_mapping" => file_sha256(joinpath(@__DIR__, "reef_observation_mapping.jl")),
    "sandbox_project" => file_sha256(joinpath(ROOT, "sandbox", "Project.toml")),
    "sandbox_manifest" => file_sha256(joinpath(ROOT, "sandbox", "Manifest.toml")),
    "script" => file_sha256(@__FILE__),
    "adapter" => file_sha256(joinpath(ROOT, "ADRIA", "src", "ecosystem", "cots.jl")),
    "scenario" => file_sha256(joinpath(ROOT, "ADRIA", "src", "scenario.jl")),
    "core" => file_sha256(joinpath(ROOT, "..", "COTSMod.jl", "src", "COTSMod.jl")),
    "owen_forcing" => file_sha256(joinpath(V2, "cots_connectivity", "datasets",
        "owen_global_2018_2023", "provenance.json")),
)
outputs = Dict(file => file_sha256(joinpath(OUT, file)) for file in
    ("design.csv", "scores.csv", "trajectories.csv", "flows.csv", "failures.csv"))
metadata = Dict(
    "status" => "opt_in_size_weighted_qualitative_screen_not_promoted",
    "supported_allee_half_saturation_adults_ha" => 3.0,
    "initial_stage_state" => initial_path,
    "capacity_pair" => PAIR,
    "site_count" => length(domain.loc_ids),
    "year_count" => length(years),
    "reef_count" => length(TARGETS),
    "q_tow_per_adult_ha" => q,
    "q_status" => "provisional coherent global conversion",
    "hydrodynamic_seed" => SEED,
    "larval_survival_assumption" => "Owen connectivity includes survival",
    "external_boundary" => "evidence-derived outbreak prediction boundary, unchanged across arms",
    "connectivity_orientation" => "Owen source-by-sink reef matrices; 2018-2023 hydrodynamic years sampled with recorded seed",
    "state_units" => "small, large, effective breeders, total adults in COTS ha^-1; annual growth and recruitment fluxes in COTS ha^-1 yr^-1",
    "female_gonad_slope_per_mm" => 0.0115,
    "female_gonad_reference" => "Babcock et al. 2016 DOI:10.1007/s00227-016-3009-5",
    "class_diameter_and_growth_status" => "diagnostic assumptions, not observed historical sizes or rates",
    "julia_version" => string(VERSION),
    "active_project" => string(Base.active_project()),
    "adria_revision" => git_revision(ROOT),
    "cotsmod_revision" => git_revision(joinpath(ROOT, "..", "COTSMod.jl")),
    "inputs" => inputs,
    "outputs" => outputs,
)
open(joinpath(OUT, "metadata.toml"), "w") do io
    TOML.print(io, metadata; sorted=true)
end
println("Saved ", OUT)
isempty(failures.treatment) || error("Size screen had failed treatments; see failures.csv")
