"""Paired diagnostic capacity screen on the 1991 Owen/Lizard domain.

This is an opt-in, one-seed, space-filling screen, not an optimizer or a change
to any model default. Each design point is run at the supported 3 COTS/ha
Allee half-saturation and at 1 COTS/ha as a labelled counterfactual. A single
provisional tow-per-adult/ha scale is used both to initialize IDW adult density
and to map modeled adults back to CPUE. No conversion is claimed validated.
"""

using Pkg
const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
Pkg.activate(joinpath(ROOT, "sandbox"))
cd(ROOT)

using ADRIA, CSV, DataFrames, Dates, Random, TOML, Statistics
include(joinpath(@__DIR__, "cots_cycle_metrics.jl"))
include(joinpath(@__DIR__, "reef_observation_mapping.jl"))
include(joinpath(@__DIR__, "calibration_run.jl"))

const HYDRO_SEED = 20260930
const DESIGN_SEED = 20261008
const N_PAIRS = parse(Int, get(ENV, "OWEN_CAPACITY_PAIRS", "24"))
1 <= N_PAIRS <= 256 || error("OWEN_CAPACITY_PAIRS must be in 1:256")
const RUN_ID = get(ENV, "OWEN_CAPACITY_RUN_ID",
    Dates.format(now(), "yyyymmddTHHMMSS") * "_owen_1991_cycle_capacity")
occursin(r"^[A-Za-z0-9_.-]+$", RUN_ID) || error("Invalid run ID")
const OUT = joinpath(@__DIR__, "runs", RUN_ID)
ispath(OUT) && error("Refusing to overwrite $OUT")
const V2 = joinpath(ROOT, "sandbox", "data", "Lizard_Historical_v0.2")
const SURFACE_DIR = joinpath(@__DIR__, "runs", "20261007T199104_1991_initial_surfaces")
const SCREEN = joinpath(@__DIR__, "runs", "20261001T162000_owen_embedded_mortality_screen")
const BEST = joinpath(@__DIR__, "runs", "expanded_cotsconn_pilot_seed20260930", "best_summary.csv")
const ARCHIVE = joinpath(@__DIR__, "runs", "20261007T199111_owen_1991_cohort_history_final")
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

function stratified(rng::AbstractRNG, n::Int, low::Float64, high::Float64; log_scale=false)
    points = (collect(0:n-1) .+ rand(rng, n)) ./ n
    shuffle!(rng, points)
    return log_scale ? exp.(log(low) .+ points .* (log(high) - log(low))) :
        low .+ points .* (high - low)
end

function design(n::Int)
    rng = MersenneTwister(DESIGN_SEED)
    q = stratified(rng, n, 0.005, 0.30; log_scale=true)
    fec = stratified(rng, n, 0.5, 2.0)
    bound = stratified(rng, n, 0.0, 1.5)
    mort = stratified(rng, n, 0.08, 0.30)
    feed = stratified(rng, n, 0.6, 1.5)
    cohort = shuffle!(rng, [Float64((i - 1) % 3) / 2 for i in 1:n])
    coral = shuffle!(rng, [isodd(i) ? "inherited" : "idw_1991" for i in 1:n])
    rows = DataFrame(pair=Int[], q_tow_per_adult_ha=Float64[],
        immature_profile_multiplier=Float64[], coral_initialization=String[],
        fecundity_multiplier=Float64[], boundary_multiplier=Float64[],
        adult_mortality=Float64[], feeding_multiplier=Float64[])
    # A paired, archived-state replay validates that the ecological baseline has
    # not changed. The new CPUE scale is intentionally q=0.015 in both directions.
    push!(rows, (0, 0.015, 1.0, "inherited", 1.0, 1.0,
        Float64(candidate.m3), 1.0))
    for i in 1:n
        push!(rows, (i, q[i], cohort[i], coral[i], fec[i], bound[i], mort[i], feed[i]))
    end
    return rows
end

surface_path = joinpath(SURFACE_DIR, "surface.csv")
surface = CSV.read(surface_path, DataFrame)
nrow(surface) == 2895 || error("Unexpected 1991 surface site count")
all(isfinite, Float64.(surface.cots_cpue_1989_1991_idw)) ||
    error("Nonfinite initial COTS CPUE")
all(>=(0.0), Float64.(surface.cots_cpue_1989_1991_idw)) ||
    error("Negative initial COTS CPUE")
screen_meta = TOML.parsefile(joinpath(SCREEN, "metadata.toml"))
fecundity = 4Float64(screen_meta["matched_fecundity_pre_pelagic_ha_year"])
boundary = 4Float64(screen_meta["matched_boundary_pre_pelagic_ha_year"])
candidate = sort(CSV.read(BEST, DataFrame), :loss)[1, :]
plan = design(N_PAIRS)
mkpath(OUT)
CSV.write(joinpath(OUT, "design.csv"), plan)

ENV["ADRIA_COTS_CONNECTIVITY_DATASET"] = "owen_global_2018_2023"
ENV["ADRIA_COTS_CONNECTIVITY_MODE"] = "cots"
ENV["COTS_CONNECTIVITY_TEMPORAL_MODE"] = "sample"
ENV["COTS_CONNECTIVITY_SEED"] = string(HYDRO_SEED)
ENV["COTS_EXTERNAL_PULSE"] = "false"
ENV["COTS_EXTERNAL_SOURCE_DENSITY"] = "0.0"
ENV["COTS_JUVENILE_STORAGE"] = "false"
ENV["COTS_HABITAT_MEDIATION"] = "false"
ENV["COTS_SEPARATED_RECRUITMENT"] = "true"
ENV["COTS_APPLY_LARVAL_SURVIVAL"] = "false" # Owen transport assumed to include it.
ENV["COTS_SETTLEMENT_PROBABILITY"] = "1.0"
ENV["COTS_IMMIGRATION_SCALAR"] = "1.0"
ENV["COTS_INITIAL_MULTIPLIER"] = string(Float64(candidate.seed_mult))
ENV["COTS_LOW_COVER_MORTALITY"] = "false"
ENV["COTS_LOW_COVER_STRENGTH"] = "0.0"
delete!(ENV, "COTS_INTERNAL_SETTLEMENT_BLACKOUT_START_YEAR")
delete!(ENV, "COTS_INTERNAL_SETTLEMENT_BLACKOUT_END_YEAR")
delete!(ENV, "COTS_INITIAL_STATE_CSV")

domain = ADRIA.load_domain(ADRIA.LizardDomain, V2, "historical")
String.(surface.site_id) == String.(domain.loc_ids) || error("Initial surface/site mismatch")
mapping = load_reef_area_weights(domain, V2, [t.model for t in TARGETS])
Random.seed!(HYDRO_SEED)
base_scenario = ADRIA.sample(domain, 2)[1:1, :]
defaults = ADRIA.param_table(domain)
for name in names(base_scenario)
    base_scenario[!, name] .= defaults[1, name]
    if name != "allee_threshold" && hasproperty(candidate, Symbol(name))
        base_scenario[!, name] .= candidate[Symbol(name)]
    end
end

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
baseline_cover = Float64.(Array(domain.init_coral_cover))
baseline_total = vec(sum(baseline_cover; dims=1))
all(>(0.0), baseline_total) || error("Zero baseline coral cover")
idw_cover = Float64.(surface.coral_manta_9m_fraction_1989_1991_idw)
new_cover = baseline_cover .* reshape(idw_cover ./ baseline_total, 1, :)
raw_obs = CSV.read(OBS, DataFrame)
archive = CSV.read(joinpath(ARCHIVE, "trajectories.csv"), DataFrame)
archive = archive[archive.treatment .== "idw_alpha015_cohort_full", :]
nrow(archive) == 34length(TARGETS) || error("Incomplete archived reference")

scores = DataFrame(treatment=String[], pair=Int[], theta_ha=Float64[],
    reef_name=String[], loss=Float64[], sim_peak_count=Int[],
    sim_peak_years=String[], sim_peak_heights_cpue=String[], matched_peaks=Int[])
trajectories = DataFrame(treatment=String[], pair=Int[], theta_ha=Float64[],
    reef_name=String[], year=Int[], forcing_year=Int[], adults_ha=Float64[],
    simulated_cpue=Float64[], coral_cover=Float64[])
flows = DataFrame(treatment=String[], pair=Int[], theta_ha=Float64[],
    reef_name=String[], year=Int[])
for name in FLOW_NAMES
    flows[!, name] = Float64[]
end
failures = DataFrame(treatment=String[], pair=Int[], theta_ha=Float64[], error=String[])

for row in eachrow(plan)
    q = Float64(row.q_tow_per_adult_ha)
    adults0 = Float64.(surface.cots_cpue_1989_1991_idw) ./ q
    mult = Float64(row.immature_profile_multiplier)
    initial = DataFrame(site_id=String.(surface.site_id),
        juvenile_ha=6.0 .* mult .* adults0,
        subadult_ha=3.0 .* mult .* adults0,
        adult_ha=adults0)
    initial_path = joinpath(OUT, "initial_state_pair_$(lpad(row.pair, 3, '0')).csv")
    CSV.write(initial_path, initial)
    ENV["COTS_INITIAL_STATE_CSV"] = initial_path
    cover = row.coral_initialization == "inherited" ? baseline_cover : new_cover
    domain.init_coral_cover = ADRIA.DataCube(cover;
        species=collect(domain.init_coral_cover.species), locations=domain.loc_ids)
    ENV["COTS_LARVAL_FECUNDITY"] = string(fecundity * row.fecundity_multiplier)
    ENV["COTS_EXTERNAL_OUTBREAK_PRODUCTION"] = string(boundary * row.boundary_multiplier)
    for theta in (3.0, 1.0)
        treatment = "pair_$(lpad(row.pair, 3, '0'))_theta_$(Int(theta))"
        scenario = copy(base_scenario)
        scenario[!, :allee_threshold] .= theta
        scenario[!, :m3] .= row.adult_mortality
        scenario[!, :a_F] .= Float64(candidate.a_F) * row.feeding_multiplier
        scenario[!, :a_S] .= Float64(candidate.a_S) * row.feeding_multiplier
        ENV["COTS_ALLEE_THRESHOLD"] = string(theta)
        Random.seed!(HYDRO_SEED)
        try
            result = ADRIA.run_scenario(domain, scenario[1, :])
            Int.(result.cots_calendar_year) == years || error("Calendar alignment failed")
            coral = dropdims(sum(result.raw; dims=(2, 3)); dims=(2, 3))
            for target in TARGETS
                reef = target.model
                m = mapping[reef]
                adult = aggregate_reef_density(Matrix(result.cots_log[:, 3, :]), m)
                reef_coral = aggregate_reef_density(coral, m)
                reef_flows = [aggregate_reef_density(
                    Matrix(result.cots_flow_log[:, channel, :]), m)
                    for channel in 1:length(FLOW_NAMES)]
                if row.pair == 0 && theta == 3.0
                    reference = sort(archive[archive.reef_name .== reef, :], :year)
                    years == Int.(reference.year) || error("Reference year mismatch")
                    Int.(result.cots_forcing_year) == Int.(reference.forcing_year) ||
                        error("Reference hydrodynamic year mismatch")
                    all(isapprox.(adult, Float64.(reference.adults_ha);
                        atol=1e-12, rtol=0)) || error("Reference adult mismatch")
                    all(isapprox.(reef_coral, Float64.(reference.coral_cover);
                        atol=1e-12, rtol=0)) || error("Reference coral mismatch")
                end
                for (idx, year) in enumerate(years)
                    push!(trajectories, (treatment, row.pair, theta, reef, year,
                        Int(result.cots_forcing_year[idx]), adult[idx], adult[idx] * q,
                        reef_coral[idx]))
                    push!(flows, (treatment, row.pair, theta, reef, year,
                        (values[idx] for values in reef_flows)...))
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
                push!(scores, (treatment, row.pair, theta, reef, score.total_loss,
                    score.n_sim_peaks, join(score.sim_peak_years, ';'),
                    join(round.([p.value for p in peaks]; digits=6), ';'),
                    score.n_matched_peaks))
            end
            println("Completed ", treatment)
        catch err
            push!(failures, (treatment, row.pair, theta,
                sprint(showerror, err, catch_backtrace())))
            println(stderr, "FAILED ", treatment, ": ", err)
        end
        flush(stdout)
        for (name, table) in (("scores.csv", scores),
            ("trajectories.csv", trajectories), ("flows.csv", flows),
            ("failures.csv", failures))
            CSV.write(joinpath(OUT, name), table)
        end
    end
end
delete!(ENV, "COTS_INITIAL_STATE_CSV")

inputs = Dict(
    "surface" => file_sha256(surface_path),
    "surface_provenance" => file_sha256(joinpath(SURFACE_DIR, "metadata.json")),
    "reference_trajectories" => file_sha256(joinpath(ARCHIVE, "trajectories.csv")),
    "reference_parameters" => file_sha256(BEST),
    "screen_metadata" => file_sha256(joinpath(SCREEN, "metadata.toml")),
    "observations" => file_sha256(OBS),
    "script" => file_sha256(@__FILE__),
    "adapter" => file_sha256(joinpath(ROOT, "ADRIA", "src", "ecosystem", "cots.jl")),
    "core" => file_sha256(joinpath(ROOT, "..", "COTSMod.jl", "src", "COTSMod.jl")),
    "owen_forcing" => file_sha256(joinpath(V2, "cots_connectivity", "datasets",
        "owen_global_2018_2023", "provenance.json")),
)
outputs = Dict(file => file_sha256(joinpath(OUT, file)) for file in
    ("design.csv", "scores.csv", "trajectories.csv", "flows.csv", "failures.csv"))
initial_state_files = ["initial_state_pair_$(lpad(pair, 3, '0')).csv"
    for pair in plan.pair]
merge!(outputs, Dict(file => file_sha256(joinpath(OUT, file))
    for file in initial_state_files))
metadata = Dict(
    "status" => "paired_capacity_screen_not_promoted",
    "hypothesis" => "current parameters plus coherent provisional tow conversion can generate two peaks with deep troughs; theta=1 is counterfactual only",
    "n_pairs" => N_PAIRS,
    "design_seed" => DESIGN_SEED,
    "hydrodynamic_seed" => HYDRO_SEED,
    "year_window" => [1991, 2024],
    "score_start_year" => 1992,
    "supported_allee_half_saturation_adults_ha" => 3.0,
    "counterfactual_allee_half_saturation_adults_ha" => 1.0,
    "q_units" => "COTS/tow per adult COTS/ha; same q used for initialization and score, provisional global factor",
    "q_range" => [0.005, 0.30],
    "fecundity_multiplier_range" => [0.5, 2.0],
    "boundary_multiplier_range" => [0.0, 1.5],
    "adult_mortality_range" => [0.08, 0.30],
    "feeding_multiplier_range" => [0.6, 1.5],
    "initial_state_rows_per_pair" => nrow(surface),
    "fecundity_reference_pre_pelagic_recruits_per_adult_year" => fecundity,
    "boundary_reference_pre_pelagic_recruits_per_ha_year" => boundary,
    "larval_survival_assumption" => "Owen connectivity already includes survival",
    "immature_initialization" => "6:3:1 inherited juvenile:subadult:adult ratio multiplied by 0, 0.5, or 1; no cohort observations exist",
    "coral_initialization" => "inherited surface or IDW 1989-1991 9m manta cover as whole-reef proxy",
    "design_bounds_status" => "diagnostic, not independent biological bounds; no automatic promotion",
    "observation_status" => "EOTR and LTMP manta assumed comparable; SALAD/manta conversion remains unvalidated",
    "size_fecundity_reference" => "Babcock et al. 2016 DOI:10.1007/s00227-016-3009-5; female gonad mass 3.384 exp(0.0115 D_mm) g, R2 0.55; not implemented because adult size is untracked",
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
isempty(failures.treatment) || error("Capacity screen had failed treatments; see failures.csv")
