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
include(joinpath(@__DIR__, "connectivity_retention.jl"))

const SEED = 20260930
const EXPERIMENT = get(ENV, "OWEN_1991_EXPERIMENT", "initialization")
EXPERIMENT in ("initialization", "cohort_history", "gap_grazing",
    "crash_factorial", "lagged_hazard", "cohort_audit", "structural_crash",
    "unit_migration_q015", "ulp_noise_floor", "unit_compat_k2",
    "second_wave_audit", "switching_condition", "switching_repair",
    "habitat_disturbance", "benthic_cyclone", "trough_mechanisms", "initiation",
    "predation_capacity", "outbreak_collapse", "growth_starvation",
    "initial_composition", "food_grace", "coral_calibration", "joint_coral",
    "mellin_coral", "bleaching_sensitivity", "cots_rebaseline", "composition_rebaseline",
    "initial_state_rebaseline", "local_retention", "bleaching_floor", "coral_refit",
    "lizard_coral", "coral_refit_v2", "initial_share", "timing", "candidates", "culling", "coral_mix") ||
    error("Unknown 1991 experiment")
"""
Insert the next per-date sequence number into an unnumbered run ID,
`YYYYMMDDT_<name>` -> `YYYYMMDDT_NN_<name>`, so runs list in the order tried
(see renumber_runs.py). Numbered or time-stamped IDs are returned unchanged.
"""
function sequenced_run_id(run_id::AbstractString)::String
    m = match(r"^(\d{8})T_(?!\d{2}_)(.+)$", run_id)
    isnothing(m) && return String(run_id)
    runs = joinpath(@__DIR__, "runs")
    taken = [parse(Int, n.captures[1]) for n in
        (match(Regex("^$(m.captures[1])T_(\\d{2})_"), name) for name in
            (isdir(runs) ? readdir(runs) : String[])) if !isnothing(n)]
    return "$(m.captures[1])T_$(lpad(isempty(taken) ? 1 : maximum(taken) + 1, 2, '0'))_$(m.captures[2])"
end
const RUN_ID = sequenced_run_id(get(ENV, "OWEN_1991_RUN_ID",
    Dates.format(now(), "yyyymmdd") * "T_" * EXPERIMENT * "_1991_v1"))
occursin(r"^[A-Za-z0-9_.-]+$", RUN_ID) || error("Invalid run ID")
println("Run ID: ", RUN_ID)
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
# Unit migration (crash_mechanism_protocol.md, 2026-10-09): q is a fixed
# COTS/tow per adult COTS/ha convention (CoCoNet epsilon; Moran & De'ath 1992
# threshold pairing). k converts the archived compressed density unit to COTS/ha
# using Lizard's frozen scale; it is exact for Lizard only.
const MIGRATION = EXPERIMENT in ("unit_migration_q015", "unit_compat_k2", "second_wave_audit",
    "switching_condition", "switching_repair", "habitat_disturbance", "benthic_cyclone",
    "trough_mechanisms", "initiation", "predation_capacity", "outbreak_collapse",
    "growth_starvation", "initial_composition", "food_grace", "coral_calibration",
    "joint_coral", "mellin_coral", "bleaching_sensitivity", "cots_rebaseline",
    "composition_rebaseline", "initial_state_rebaseline", "local_retention", "bleaching_floor",
    "coral_refit", "lizard_coral", "coral_refit_v2", "initial_share", "timing", "candidates", "culling", "coral_mix")
# Preregistered first-wave initiation screen on the adopted v0.3 domain (2026-10-12):
# external supply x5 and/or 1991 juvenile preload x5, crossed with none / A / P+A.
const INITIATION = EXPERIMENT == "initiation"
# Preregistered predator-pit capacity tuning grid (2026-10-12): External x5 (working
# assumption) + gradual Allee, with P off or Q on a fixed log-spaced grid.
const PIT = EXPERIMENT == "predation_capacity"
const PIT_CAPACITIES = (20, 30, 45, 70, 100, 150, 220, 330, 500)
# Preregistered post-outbreak collapse x fast-coral growth screen (2026-10-12):
# leading candidate (Ext x5 + P(Q=30) + A) x ReefMod-style collapse x coral
# calibration group (Lizard default 1 vs per-reef RME layer, +/- fast extension x1.2).
const COLLAPSE = EXPERIMENT == "outbreak_collapse"
const COLLAPSE_Q = 30.0
const COLLAPSE_THRESHOLD_HA = 67.5   # ReefMod 2.7 per 400 m^2
const COLLAPSE_BACKGROUND_HA = 0.25  # ReefMod 0.01 per 400 m^2
const FAST_LINEAR_EXTENSION_MULTIPLIER = 1.2
const RME_GBR_LAYER = joinpath(ROOT, "sandbox", "data", "rme_ml_2025_06_05", "data_files",
    "region", "reefmod_gbr.gpkg")
# Preregistered fast-coral growth fit x faster starvation screen (2026-10-12): RME
# calibration groups (project lead: group 5 for the calibration reefs), fast-group
# linear-extension multiplier grid x starvation {ref, memory 0, post-grazing}.
const GROWTH = EXPERIMENT == "growth_starvation"
const GROWTH_MULTIPLIERS = (1.0, 1.05, 1.10, 1.15, 1.20)
growth_arm_multiplier(name) = parse(Int, split(name, "_")[1][2:end]) / 100
# Preregistered data-derived 1991 coral composition screen (2026-10-12): fast:slow
# split per reef from 1993-95 LTMP photo transects, totals and within-block structure
# kept; Q re-checked; fast-growth multiplier sensitivity.
const COMPOSITION = EXPERIMENT == "initial_composition"
const COMPOSITION_MULTIPLIERS = (1.0, 1.10, 1.20)
const PHOTO_TRANSECT_CSV = joinpath(ROOT, "sandbox", "data", "reef_photo_transect.csv")
const PHOTO_REEFS = Dict("Lizard Island Reef" => "Lizard Isles",
    "MacGillivray Reef" => "Macgillivray Reef", "North Direction Reef" => "North Direction Island")
const PHOTO_FAST_TAXA = ("Acropora", "Pocilloporidae")
const PHOTO_YEARS = 1993:1995
const FAST_GROUP_NAMES = ("tabular_Acropora", "corymbose_Acropora", "corymbose_non_Acropora")
composition_arm_q(name) = parse(Float64, split(name, "_q")[end])
# Preregistered food-scale screen (2026-10-12): plausible 1991 fast share (photo
# transect Acropora + Pocilloporidae + w * "Other", w = 0, 0.5, 1) x starvation grace
# {0, 1} year, on the run-14 best configuration (fast extension x1.0, Q = 30).
const GRACE = EXPERIMENT == "food_grace"
const OTHER_WEIGHTS = Dict("acrpoc" => 0.0, "halfother" => 0.5, "other" => 1.0)
# Preregistered coral-side calibration (2026-10-12): COTS grazing forced from observed
# CPUE / q (COTS dynamics off); slow (massive) x fast linear-extension multipliers on
# every RME calibration group; observed 1993-95 composition.
const CORAL = EXPERIMENT == "coral_calibration"
const SLOW_MULTIPLIERS = (0.40, 0.55, 0.70, 0.85, 1.00)
const FAST_MULTIPLIERS = (1.0, 1.2, 1.4)
const FORCED_ADULTS_CSV = joinpath(ROOT, "sandbox", "data", "derived",
    "observed_cots_forcing_1991_2024", "observed_adults_ha.csv")
coral_arm_name(slow, fast) =
    "forced_s$(lpad(round(Int, 100 * slow), 3, '0'))_f$(round(Int, 100 * fast))"
# Preregistered joint coral fit (2026-10-13): growth (s, f) x cyclone strength kappa x
# background mortality mu, COTS forced from observed CPUE / q. Mortality per site-year:
# m' = 1 - (1 - m)^kappa * (1 - mu) through the cyclone mortality step.
const JOINT = EXPERIMENT == "joint_coral"
const JOINT_SLOW = (0.55, 0.75, 1.00)
const JOINT_FAST = (0.8, 1.0)
const JOINT_KAPPA = (1.0, 2.5, 5.0)
const JOINT_BACKGROUND = (0.0, 0.03, 0.06)
pct3(x) = lpad(round(Int, 100 * x), 3, '0')
# Minimal reader for the {"reef": {"level": value, ...}, ...} JSON written by
# build_initial_share_range.py (avoids adding a JSON dependency to the sandbox runner).
function parse_share_levels(path::AbstractString)::Dict{String,Dict{String,Float64}}
    text = read(path, String)
    out = Dict{String,Dict{String,Float64}}()
    for block in eachmatch(r"\"([^\"]+)\"\s*:\s*\{([^{}]*)\}", text)
        out[String(block.captures[1])] = Dict(String(m.captures[1]) => parse(Float64, m.captures[2])
            for m in eachmatch(r"\"([^\"]+)\"\s*:\s*([-0-9.eE+]+)", block.captures[2]))
    end
    return out
end
joint_arm_name(slow, fast, kappa, mu) =
    "forced_s$(pct3(slow))_f$(pct3(fast))_k$(pct3(kappa))_b$(pct3(mu))"
joint_arm_values(name) = Tuple(parse(Int, part[2:end]) / 100 for part in split(name, "_")[2:5])
# Preregistered Mellin follow-up (2026-10-13): as the joint fit, but the background
# mortality is the Mellin et al. 2019 / COTS_Model site x year table.
const MELLIN = EXPERIMENT == "mellin_coral"
const MELLIN_BACKGROUND_CSV = joinpath(ROOT, "sandbox", "data", "derived",
    "mellin_background_mortality_1991_2024", "background_mortality.csv")
mellin_arm_name(slow, fast, kappa) = "mellin_s$(pct3(slow))_f$(pct3(fast))_k$(pct3(kappa))"
# Preregistered bleaching-response sensitivity (2026-10-13): CoralBlox growth untouched
# (s = f = 1), Mellin background with or without 3%/yr chronic loss, cyclone kappa, and
# ADRIA_DEBUG_BLEACHING_SCALAR beta; COTS forced from observed CPUE / q.
# The coral re-fit (2026-10-14) is a member of the bleaching family (forced COTS, Mellin
# background, kappa 5, beta 1.5, growth untouched): bleaching floor {1, 3, 4} x chronic
# background mu {0.03, 0.045, 0.06, 0.08} x forcing q {0.015, 0.030}; arms
# refit_f<floor>_mu<mu>_q<q>; refit_f100_mu030_q015 replays run 09 bleach_mel3_k500_b150.
# Coral re-fit v2 (2026-10-14): low-DHW bleaching scale s {1, 0.25, 0} below 4 DHW x chronic
# mu {0.06..0.14} x forcing alignment {rep: reporting-year table, srv: survey-year table};
# arms refit2_s<s>_mu<mu>_<align>; refit2_s100_mu030_rep replays run 09.
const REFIT2 = EXPERIMENT == "coral_refit_v2"
const REFIT2_SCALES = (1.0, 0.25, 0.0)
const REFIT2_MU = (0.06, 0.08, 0.10, 0.12, 0.14)
const FORCED_ADULTS_SURVEY_CSV = joinpath(ROOT, "sandbox", "data", "derived",
    "observed_cots_forcing_1991_2024_survey_year", "observed_adults_ha.csv")
const REFIT = EXPERIMENT == "coral_refit" || REFIT2
const REFIT_FLOORS = (1.0, 3.0, 4.0)
const REFIT_MU = (0.03, 0.045, 0.06, 0.08)
const REFIT_Q = (0.015, 0.030)
const REFIT_KAPPA = 5.0
const BLEACH = EXPERIMENT == "bleaching_sensitivity" || REFIT
refit_arm_name(floor, mu, q) = "refit_f$(pct3(floor))_mu$(lpad(round(Int, 1000 * mu), 3, '0'))_q$(lpad(round(Int, 1000 * q), 3, '0'))"
const BLEACH_BACKGROUND = ("mel" => 0.0, "mel3" => 0.03)
const BLEACH_KAPPA = (2.5, 5.0)
const BLEACH_BETA = (1.0, 1.5, 2.5)
bleach_arm_name(level, kappa, beta) = "bleach_$(level)_k$(pct3(kappa))_b$(pct3(beta))"
# Preregistered COTS re-baseline (2026-10-13): run-15 leading candidate (COTS dynamics
# on) on the lead's working coral configuration: growth untouched, Mellin + 3%/yr,
# kappa 5, bleaching beta {1, 1.5}, with Q re-checked; ref and single-component arms.
# The composition follow-up (2026-10-13) reuses this machinery: working coral, Q = 20,
# 1991 fast share {acrpoc, halfother, other} (run-15 levels) x beta {1, 1.5}; arms are
# named wc_b<beta>_q20_<level>.
const REBASE_COMPOSITION = EXPERIMENT == "composition_rebaseline"
# The starting-state follow-up (2026-10-13) also reuses it: 1991 CPUE surface {idw, eb
# (empirical-Bayes shrunk donors)} x inversion {lin: C/q, h10: inverse of the
# density-dependent manta model at A_half = 10} x Q {20, 30}; arms wc_b150_q<Q>_<surface>_<inv>.
const REBASE_INITIAL = EXPERIMENT == "initial_state_rebaseline"
const EB_SURFACE_DIR = joinpath(@__DIR__, "runs", "20261013T_15_1991_initial_surfaces_eb")
const INITIAL_A_HALF = 10.0
# Local-retention screen (2026-10-13): adopted eb/lin starting state; Owen annual matrices
# with source-conserving retention (connectivity_retention.jl) in the Lizard complex
# (14-116b/c/d) at rho 0.25/0.5/0.9, or on every reef at 0.5; arms wc_b150_q20_eb_lin_<ret>.
const REBASE_RETENTION = EXPERIMENT == "local_retention"
const RETENTION_LEVELS = Dict("owen" => nothing, "liz25" => ("lizard", 0.25),
    "liz50" => ("lizard", 0.5), "liz90" => ("lizard", 0.9), "all50" => ("all", 0.5))
const LIZARD_COMPLEX_IDS = ("14-116b", "14-116c", "14-116d")
# Bleaching-floor screen (2026-10-14): ADRIA_BLEACHING_DHW_FLOOR {1 (default), 4} x beta
# {1, 1.5} x starting fast-share multiplier m {1, 2, 3} on the observed acrpoc share
# (capped 0.9); eb/lin seeds; arms wc_b<beta>_q20_eb_lin_f<floor>_m<m>.
const REBASE_FLOOR = EXPERIMENT == "bleaching_floor"
const FLOOR_LEVELS = (1.0, 4.0)
const SHARE_MULTIPLIERS = (1.0, 2.0, 3.0)
const SHARE_CAP = 0.9
# Lizard hard-coral initialisation screen (2026-10-14): Lizard 1991 cover x c_L {1, 1.5, 2}
# (Pratchett 2005 reef-wide 32% vs LTMP ~13%) x observation q {0.015, 0.030} (eb seeds
# inverted with the same q) x chronic mu {0.03, 0.08}; arms
# wc_b150_q20_eb_lin_obs<q>_mu<mu>_c<c_L>.
const REBASE_LIZCORAL = EXPERIMENT == "lizard_coral"
const LIZCORAL_Q = (0.015, 0.030)
const LIZCORAL_MU = (0.03, 0.08)
const LIZCORAL_SCALES = (1.0, 1.5, 2.0)
# Starting fast-share partial calibration (2026-10-14): coral re-fit v2 configuration (cv2:
# low-DHW bleaching x0.25 below 4 DHW, Mellin + 12%/yr) x 1991 fast share level {min, q25, q50,
# q75, max} of each reef's healthy-year photo-transect range (build_initial_share_range.py);
# arms wc_b150_q20_eb_lin_ref (old coral, replays run 17), ..._cv2_acrpoc, ..._cv2_<level>,
# and ..._cv2_pick when OWEN_INITIAL_SHARE_PICKS="Lizard Island Reef=q75;..." is set.
# Timing factorial (2026-10-14, close-out B2): cv2 x predator-pit Q {10, 15, 20} x share set
# {pick: Lizard min, MacGillivray q50, North Direction q25; pnd: North Direction min} x observation
# q {0.015, 0.030} (eb seeds inverted with the arm's q); arms wc_b150_q<Q>_eb_lin_obs<q>_cv2_<set>.
# Reuses the initial_share paths; wc_b150_q20_eb_lin_obs015_cv2_pick replays run 25's cv2_pick.
const TIMING = EXPERIMENT == "timing"
const TIMING_Q = (10, 15, 20)
const TIMING_OBS_Q = (0.015, 0.030)
const TIMING_SHARE_SETS = Dict(
    "pick" => Dict("Lizard Island Reef" => "min", "MacGillivray Reef" => "q50", "North Direction Reef" => "q25"),
    "pnd" => Dict("Lizard Island Reef" => "min", "MacGillivray Reef" => "q50", "North Direction Reef" => "min"),
    "pndl50" => Dict("Lizard Island Reef" => "q50", "MacGillivray Reef" => "q50", "North Direction Reef" => "min"))
timing_obs_q(name) = parse(Int, match(r"_obs(\d+)_", name).captures[1]) / 1000
# Candidate set x replicate seeds (2026-10-14, close-out B3/B4): candidates (suffix, observation q)
# x five (COTS_CONNECTIVITY_SEED, ADRIA_RNG_SEED) pairs; pair 1 is the archived pair, so s1 arms
# replay runs 25 / 32. Arms wc_b150_q20_eb_lin_obs<q>_s<k>_<suffix>.
const CANDIDATES = EXPERIMENT == "candidates"
const CANDIDATE_ARMS = (("cv2_pnd", 0.015), ("cv2_pick", 0.015), ("cv2_pnd", 0.030),
    ("cv2_pndl50", 0.015), ("cv2_pndl50", 0.030), ("ref", 0.015))
const CANDIDATE_SEEDS = ((20260930, 90972102), (20261101, 11101), (20261102, 22202),
    (20261103, 33303), (20261104, 44404))
candidate_seed_pair(name) = CANDIDATE_SEEDS[parse(Int, match(r"_s(\d+)_", name).captures[1])]
# Observed culling (2026-10-14, close-out C2): working configuration (cv2_pnd, q 0.015, seed
# pair 1) x {no cull, observed cull with fuzzy site matching sigma 0 / 210 / 420 m}; arms
# wc_b150_q20_eb_lin_obs015_s1_<cull>_cv2_pnd. The nocull arm replays the candidate W arm.
const CULLING = EXPERIMENT == "culling"
const CULL_DIR = joinpath(ROOT, "sandbox", "data", "derived", "cots_cull_removals_lizard_v0_3_v2")
const CULL_ARMS = ("nocull", "cull000", "cull210", "cull420")
cull_table(name) = (m = match(r"_cull(\d{3})_", name); isnothing(m) ? nothing :
    joinpath(CULL_DIR, "removals_sigma$(m.captures[1]).csv"))
# Intermediate coral settings (2026-10-14, sixth round): working configuration (cv2_pnd, q 0.015,
# seed pair 1) x chronic mu {0.03, 0.06, 0.09, 0.12} x bleaching multiplier below 4 DHW
# {1, 0.5, 0.25}, plus the run-17 coral control; arms wc_b150_q20_eb_lin_obs015_s1_m<mu>_l<scale>_cv2_pnd
# and wc_b150_q20_eb_lin_obs015_s1_ref. (0.12, 0.25) replays W; ref replays C5.
const CORALMIX = EXPERIMENT == "coral_mix"
const CORALMIX_MU = (0.03, 0.06, 0.09, 0.12)
const CORALMIX_LOW = (1.0, 0.5, 0.25)
coralmix_parts(name) = (m = match(r"_m(\d{3})_l(\d{3})_", name);
    (mu=parse(Int, m.captures[1]) / 1000, low=parse(Int, m.captures[2]) / 100))
const REBASE_SHARE = EXPERIMENT in ("initial_share", "timing", "candidates", "culling", "coral_mix")
const SHARE_LEVELS_JSON = joinpath(ROOT, "sandbox", "data", "derived", "initial_fast_share_range", "share_levels.json")
const SHARE_LEVEL_NAMES = ("min", "q25", "q50", "q75", "max")
const CV2_MU = 0.12
const CV2_LOW_DHW_SCALE = 0.25
const SHARE_PICKS = get(ENV, "OWEN_INITIAL_SHARE_PICKS", "")
TIMING && !isempty(SHARE_PICKS) && error("OWEN_INITIAL_SHARE_PICKS is not used by the timing experiment")
const REBASE = EXPERIMENT == "cots_rebaseline" || REBASE_COMPOSITION || REBASE_INITIAL || REBASE_RETENTION ||
    REBASE_FLOOR || REBASE_LIZCORAL || REBASE_SHARE
lizcoral_parts(name) = (q=parse(Int, split(name, "_")[6][4:end]) / 1000,
    mu=parse(Int, split(name, "_")[7][3:end]) / 1000, c=parse(Int, split(name, "_")[8][2:end]) / 100)
const REBASE_KAPPA = 5.0
const REBASE_CHRONIC = 0.03
const REBASE_BETA = (1.0, 1.5)
const REBASE_Q = (10, 20, 30, 45, 70)
const V3RUN = INITIATION || PIT || COLLAPSE || GROWTH || COMPOSITION || GRACE || CORAL || JOINT ||
    MELLIN || BLEACH || REBASE
const V3 = joinpath(ROOT, "sandbox", "data", "Lizard_Historical_v0.3")
const DOMAIN_PATH = V3RUN ? V3 : V2
const EXTERNAL_MULTIPLIER = 5.0
const JUVENILE_PRELOAD_MULTIPLIER = 5.0
# Preregistered trough-mechanism screen (2026-10-12): predator pit (P), preferred-prey
# maturation delay (D), gradual Allee (A) on k=0.5 and geomorphic-k + cyclone bases.
const TROUGH = EXPERIMENT == "trough_mechanisms"
const GEOMORPHIC_K_CSV = joinpath(ROOT, "sandbox", "data", "Lizard_Historical_v0.2", "derived",
    "site_habitable_fraction_benthic_geomorphic2025.csv")
const JUVENILE_PREDATION_CAPACITY = 20.0
# Preregistered benthic habitable fraction x historical cyclone screen (2026-10-12).
const BENTHIC_CYCLONE = EXPERIMENT == "benthic_cyclone"
const DERIVED = joinpath(ROOT, "sandbox", "data", "Lizard_Historical_v0.2", "derived")
const BENTHIC_K_CSV = joinpath(DERIVED, "site_habitable_fraction_benthic2025.csv")
const CYCLONE_HOURS_CSV = joinpath(DERIVED, "cyclone_hours4m_rrn2025.csv")
# Survival exp(b*H) per hour of >=4 m waves: posterior median of the LTMP
# dose-response (COTS-Coral-OutcomesV2 q3_doseresponse_coefs.csv, -0.0444 / 10 h).
const CYCLONE_LOG_SURVIVAL_PER_HOUR = -0.00444
# ADRIA cannot run k = 0 sites (scenario.jl sizes habitable arrays with
# length(habitable_locs) instead of count, and constrain_recruitment! mixes
# all-site and habitable-site vectors). Benthic arms floor k at this value and
# keep zero initial cover where the benthic map gives k = 0.
const BENTHIC_K_FLOOR = 0.01
# Preregistered habitat x bleaching screen on the condition model (2026-10-12):
# RME habitable fraction (1 - sand proportion) versus the k=0.5 placeholder, and
# a labelled no-bleaching counterfactual. Cyclone/wave forcing remains zero.
const HABITAT = EXPERIMENT == "habitat_disturbance"
const RME_ID_LIST = joinpath(ROOT, "sandbox", "data", "rme_ml_2025_06_05", "data_files",
    "id", "id_list_2024_12_01.csv")
# Preregistered 2x2 on the per-capita structural base (protocol, 2026-10-12).
const SWITCHING_CONDITION = EXPERIMENT == "switching_condition"
# Preregistered prey-switching repair (protocol, 2026-10-12): diet-linked
# starvation x non-preferred feeding efficiency, all with preferred condition.
const SWITCHING_REPAIR = EXPERIMENT == "switching_repair"
const NONPREFERRED_EFFICIENCY = 0.2
const SWITCH_MIDPOINT = 0.10
const SWITCH_STEEPNESS = 70.0
const CONDITION_SATURATION = 0.10
const CONDITION_MEMORY_YEARS = 1.0
# Readout-only replay of two v2 arms with site-level second-wave diagnostics.
const SECOND_WAVE_AUDIT = EXPERIMENT == "second_wave_audit"
const Q_FIXED = 0.015
# unit_compat_k2 is a numerical diagnostic: k=2 makes every rescaling exact in
# binary floating point, isolating rounding effects from unscaled terms.
const K_MIGRATION = EXPERIMENT == "unit_compat_k2" ?
    parse(Float64, get(ENV, "OWEN_1991_COMPAT_K", "2.0")) :
    fixed_scales["Lizard Island Reef"] / Q_FIXED
const PER_CAPITA_M2_PER_ADULT_YEAR = 10.0

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

domain = ADRIA.load_domain(ADRIA.LizardDomain, DOMAIN_PATH, "historical")
baseline_k = Float64.(domain.loc_data.k)
rme_k = if HABITAT
    id_table = CSV.read(RME_ID_LIST, DataFrame; comment="#", header=false,
        types=Dict(1 => String))
    lookup = Dict(String(row[1]) => 1.0 - Float64(row[3]) for row in eachrow(id_table))
    site_rme_k = [get(lookup, String(id), NaN) for id in domain.loc_data.GBRMPA_ID]
    all(isfinite, site_rme_k) && all(0.0 .< site_rme_k .<= 1.0) ||
        error("RME habitable fraction missing or invalid for a domain site")
    site_rme_k
else
    baseline_k
end
geomorphic_k = if TROUGH
    table = CSV.read(GEOMORPHIC_K_CSV, DataFrame; types=Dict(:site_id => String))
    lookup = Dict(row.site_id => Float64(row.k_geomorphic_infill) for row in eachrow(table))
    all(haskey(lookup, String(id)) for id in domain.loc_ids) ||
        error("Geomorphic habitable fraction missing for a domain site")
    raw_geo = [lookup[String(id)] for id in domain.loc_ids]
    (k=max.(raw_geo, BENTHIC_K_FLOOR), zero=raw_geo .== 0.0)
else
    nothing
end
benthic = if BENTHIC_CYCLONE
    table = CSV.read(BENTHIC_K_CSV, DataFrame; types=Dict(:site_id => String))
    lookup = Dict(row.site_id => (Float64(row.k_benthic), Float64(row.k_benthic_mapped))
        for row in eachrow(table))
    all(haskey(lookup, String(id)) for id in domain.loc_ids) ||
        error("Benthic habitable fraction missing for a domain site")
    raw_total = [lookup[String(id)][1] for id in domain.loc_ids]
    raw_mapped = [lookup[String(id)][2] for id in domain.loc_ids]
    (total=max.(raw_total, BENTHIC_K_FLOOR), mapped=max.(raw_mapped, BENTHIC_K_FLOOR),
     total_zero=raw_total .== 0.0, mapped_zero=raw_mapped .== 0.0)
else
    nothing
end
String.(surface.site_id) == String.(domain.loc_ids) ||
    error("1991 surface order does not match model site order")
# Starting-state screen: empirical-Bayes surface (build_1991_eb_surface.py), same sites/order.
eb_surface = if REBASE_INITIAL || REBASE_RETENTION || REBASE_FLOOR || REBASE_LIZCORAL || REBASE_SHARE
    table = CSV.read(joinpath(EB_SURFACE_DIR, "surface.csv"), DataFrame)
    String.(table.site_id) == String.(domain.loc_ids) || error("EB surface order does not match model sites")
    # pandas' default float parser can differ from Julia's by one ulp on the copied column.
    all(isapprox.(Float64.(table.cots_cpue_1989_1991_idw_unshrunk),
        Float64.(surface.cots_cpue_1989_1991_idw); rtol=0.0, atol=1e-12)) ||
        error("EB surface was not built from the current IDW surface")
    eb_cpue = Float64.(table.cots_cpue_1989_1991_idw)
    all(isfinite, eb_cpue) && all(>=(0.0), eb_cpue) || error("Invalid EB surface")
    eb_cpue
else
    nothing
end
# Local-retention screen: original annual matrices (kept by reference so the owen arm
# replays exactly) and source-conserving retention variants (rows sources, columns sinks).
base_reef_connectivity = REBASE_RETENTION ? copy(domain.cots_forcing.reef_connectivity) : nothing
retention = if REBASE_RETENTION
    site_reef = domain.cots_forcing.site_to_reef
    site_gbrmpa = String.(domain.loc_data.GBRMPA_ID)
    lizard = sort(unique(site_reef[in.(site_gbrmpa, Ref(LIZARD_COMPLEX_IDS))]))
    length(lizard) == length(LIZARD_COMPLEX_IDS) || error("Lizard complex reefs not found in COTS forcing")
    n_reefs = size(first(base_reef_connectivity), 1)
    group_sets = Dict("lizard" => [lizard],
        "all" => vcat([lizard], [[r] for r in 1:n_reefs if !(r in lizard)]))
    matrices = Dict(level => [retain_within_groups(m, group_sets[spec[1]], spec[2])
                              for m in base_reef_connectivity]
        for (level, spec) in RETENTION_LEVELS if !isnothing(spec))
    row_delta = maximum(maximum(abs.(sum(new; dims=2) .- sum(old; dims=2)))
        for mats in values(matrices) for (new, old) in zip(mats, base_reef_connectivity))
    row_delta <= 1e-12 || error("Retention changed source survival: $row_delta")
    all(all(>=(0.0), m) for mats in values(matrices) for m in mats) || error("Negative connectivity")
    (matrices=matrices, lizard=lizard, row_delta=row_delta,
     owen_lizard_share=[retained_share(m, lizard) for m in base_reef_connectivity])
else
    nothing
end
retention_matrices = REBASE_RETENTION ? retention.matrices : nothing
mapping = load_reef_area_weights(domain, DOMAIN_PATH, [t.model for t in TARGETS])
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
# Package-declared cyclone mortality (v0.3), sliced to the run window.
domain_cyclones = ADRIA.DataCube(Float64.(Array(domain.cyclone_mortality_scens)[keep, :, :, :]);
    timesteps=years, locs=domain.loc_ids,
    species=ADRIA.functional_group_names(), scenarios=[1])
domain.cyclone_mortality_scens = ADRIA.ZeroDataCube(; T=Float64,
    timesteps=years, locs=domain.loc_ids,
    species=ADRIA.functional_group_names(), scenarios=[1])
zero_cyclones = domain.cyclone_mortality_scens
# Joint coral fit: cyclone strength kappa and background mortality mu per arm; the
# kappa = 1, mu = 0 arm reuses the domain cube so it replays earlier runs exactly.
joint_cyclones = if JOINT
    base_mortality = Array(domain_cyclones)
    Dict((kappa, mu) => kappa == 1.0 && mu == 0.0 ? domain_cyclones :
        ADRIA.DataCube(1.0 .- (1.0 .- base_mortality) .^ kappa .* (1.0 - mu);
            timesteps=years, locs=domain.loc_ids,
            species=ADRIA.functional_group_names(), scenarios=[1])
        for kappa in JOINT_KAPPA for mu in JOINT_BACKGROUND)
else
    nothing
end
mellin_cyclones = if MELLIN || BLEACH || REBASE
    background_table = CSV.read(MELLIN_BACKGROUND_CSV, DataFrame; types=Dict(:site_id => String))
    background_rows = Dict(String(id) => i for (i, id) in enumerate(background_table.site_id))
    all(haskey(background_rows, String(id)) for id in domain.loc_ids) ||
        error("Mellin background table does not cover every domain site")
    background_columns = String[string(Int(year)) for year in collect(years)]
    all(in(names(background_table)), background_columns) ||
        error("Mellin background table does not cover every year")
    background = permutedims(Matrix{Float64}(background_table[
        Int[background_rows[String(id)] for id in domain.loc_ids], background_columns]))
    all(0.0 .<= background .< 1.0) || error("Mellin background mortality outside [0, 1)")
    base_mortality = Array(domain_cyclones)
    Dict(kappa => ADRIA.DataCube(
            1.0 .- (1.0 .- base_mortality) .^ kappa .* (1.0 .- reshape(background,
                size(background, 1), size(background, 2), 1, 1));
            timesteps=years, locs=domain.loc_ids,
            species=ADRIA.functional_group_names(), scenarios=[1])
        for kappa in (MELLIN ? JOINT_KAPPA : REBASE ? (1.0, REBASE_KAPPA) : BLEACH_KAPPA))
else
    nothing
end
# Bleaching screen: kappa x background level (mel reuses the Mellin cube exactly).
bleach_cyclones = if BLEACH && !REFIT
    Dict((level, kappa) => mu == 0.0 ? mellin_cyclones[kappa] :
        ADRIA.DataCube(1.0 .- (1.0 .- Array(mellin_cyclones[kappa])) .* (1.0 - mu);
            timesteps=years, locs=domain.loc_ids,
            species=ADRIA.functional_group_names(), scenarios=[1])
        for (level, mu) in BLEACH_BACKGROUND for kappa in BLEACH_KAPPA)
else
    nothing
end
BLEACH && !REFIT && !all(cube -> all(0.0 .<= Array(cube) .<= 1.0), values(bleach_cyclones)) &&
    error("Bleaching-screen mortality outside [0, 1]")
# Coral re-fit: same expression as the bleaching screen's mel3 cube, so mu = 0.03 replays it.
refit_cyclones = if REFIT
    Dict(mu => ADRIA.DataCube(1.0 .- (1.0 .- Array(mellin_cyclones[REFIT_KAPPA])) .* (1.0 - mu);
            timesteps=years, locs=domain.loc_ids,
            species=ADRIA.functional_group_names(), scenarios=[1])
        for mu in (REFIT2 ? (0.03, REFIT2_MU...) : REFIT_MU))
else
    nothing
end
REFIT && !all(cube -> all(0.0 .<= Array(cube) .<= 1.0), values(refit_cyclones)) &&
    error("Re-fit mortality outside [0, 1]")
# Forcing q = 0.030 (Pratchett 2005 pairs): the q = 0.015 adults table scaled by 0.015 / q.
refit_forcing = if REFIT2
    Dict("rep" => FORCED_ADULTS_CSV, "srv" => FORCED_ADULTS_SURVEY_CSV)
elseif REFIT
    table = CSV.read(FORCED_ADULTS_CSV, DataFrame; types=Dict(:site_id => String))
    out = Dict(0.015 => FORCED_ADULTS_CSV)
    for q in REFIT_Q
        q == 0.015 && continue
        scaled = copy(table)
        for column in names(scaled)[2:end]
            scaled[!, column] = Float64.(scaled[!, column]) .* (0.015 / q)
        end
        path = joinpath(OUT, "forced_adults_q$(lpad(round(Int, 1000 * q), 3, '0')).csv")
        CSV.write(path, scaled)
        out[q] = path
    end
    out
else
    nothing
end
# Re-baseline: ref reuses the domain cube so it replays run 15; k500 is cyclones only,
# mel3 is Mellin + 3%/yr only (kappa 1), wc is the working coral configuration.
rebase_cyclones = if REBASE
    rebase_cube(mortality) = ADRIA.DataCube(mortality; timesteps=years, locs=domain.loc_ids,
        species=ADRIA.functional_group_names(), scenarios=[1])
    Dict("ref" => domain_cyclones,
        "k500" => rebase_cube(1.0 .- (1.0 .- Array(domain_cyclones)) .^ REBASE_KAPPA),
        "mel3" => rebase_cube(1.0 .- (1.0 .- Array(mellin_cyclones[1.0])) .* (1.0 - REBASE_CHRONIC)),
        "wc" => rebase_cube(1.0 .- (1.0 .- Array(mellin_cyclones[REBASE_KAPPA])) .* (1.0 - REBASE_CHRONIC)))
else
    nothing
end
REBASE && !all(cube -> all(0.0 .<= Array(cube) .<= 1.0), values(rebase_cyclones)) &&
    error("Re-baseline mortality outside [0, 1]")
lizcoral_cyclones = REBASE_LIZCORAL ? Dict(mu => mu == REBASE_CHRONIC ? rebase_cyclones["wc"] :
        ADRIA.DataCube(1.0 .- (1.0 .- Array(mellin_cyclones[REBASE_KAPPA])) .* (1.0 - mu);
            timesteps=years, locs=domain.loc_ids,
            species=ADRIA.functional_group_names(), scenarios=[1])
    for mu in LIZCORAL_MU) : nothing
REBASE_LIZCORAL && !all(cube -> all(0.0 .<= Array(cube) .<= 1.0), values(lizcoral_cyclones)) &&
    error("Lizard-coral mortality outside [0, 1]")
cv2_cyclones = REBASE_SHARE ? ADRIA.DataCube(1.0 .- (1.0 .- Array(mellin_cyclones[REBASE_KAPPA])) .* (1.0 - CV2_MU);
    timesteps=years, locs=domain.loc_ids, species=ADRIA.functional_group_names(), scenarios=[1]) : nothing
REBASE_SHARE && !all(0.0 .<= Array(cv2_cyclones) .<= 1.0) && error("cv2 mortality outside [0, 1]")
# Same expression as cv2 / run 17 ("wc"), with mu varied.
coralmix_cyclones = CORALMIX ? Dict(mu => ADRIA.DataCube(1.0 .- (1.0 .- Array(mellin_cyclones[REBASE_KAPPA])) .* (1.0 - mu);
    timesteps=years, locs=domain.loc_ids, species=ADRIA.functional_group_names(), scenarios=[1])
    for mu in CORALMIX_MU) : nothing
CORALMIX && Array(coralmix_cyclones[CV2_MU]) != Array(cv2_cyclones) && error("mu 0.12 does not reproduce cv2")
CORALMIX && Array(coralmix_cyclones[REBASE_CHRONIC]) != Array(rebase_cyclones["wc"]) && error("mu 0.03 does not reproduce run 17")
MELLIN && !all(cube -> all(0.0 .<= Array(cube) .<= 1.0), values(mellin_cyclones)) &&
    error("Mellin combined mortality outside [0, 1]")
JOINT && !all(cube -> all(0.0 .<= Array(cube) .<= 1.0), values(joint_cyclones)) &&
    error("Joint mortality outside [0, 1]")
V3RUN && maximum(Array(domain_cyclones)) > 0.0 ||
    !V3RUN || error("v0.3 domain did not provide cyclone mortality")
cyclone_cube = if BENTHIC_CYCLONE || TROUGH
    hours_table = CSV.read(CYCLONE_HOURS_CSV, DataFrame; types=Dict(:reef_id => String))
    hours = Dict(String(row.reef_id) => row for row in eachrow(hours_table))
    reef_ids = String.(domain.loc_data.GBRMPA_ID)
    all(haskey(hours, id) for id in reef_ids) || error("Cyclone hours missing for a domain reef")
    n_groups = length(ADRIA.functional_group_names())
    mortality = zeros(Float64, length(years), length(reef_ids), n_groups, 1)
    for (yi, year) in enumerate(years), (si, id) in enumerate(reef_ids)
        value = hours[id][Symbol(string(year))]
        h = ismissing(value) || !isfinite(Float64(value)) ? 0.0 : Float64(value)
        mortality[yi, si, :, 1] .= 1.0 - exp(CYCLONE_LOG_SURVIVAL_PER_HOUR * h)
    end
    all(0.0 .<= mortality .<= 1.0) || error("Cyclone mortality outside [0, 1]")
    ADRIA.DataCube(mortality; timesteps=years, locs=domain.loc_ids,
        species=ADRIA.functional_group_names(), scenarios=[1])
else
    zero_cyclones
end
md = domain.env_layer_md
domain.env_layer_md = ADRIA.EnvLayer(md.dpkg_path, md.loc_data_fn, md.loc_id_col,
    md.cluster_id_col, md.init_coral_cov_fn, md.connectivity_fn, md.DHW_fn,
    md.wave_fn, years)
baseline_cover = Float64.(Array(domain.init_coral_cover))
size(baseline_cover, 2) == nrow(surface) || error("Unexpected initial coral dimensions")
baseline_total = vec(sum(baseline_cover; dims=1))
# v0.3 holds zero cover at non-habitable sites; the manta rescaling (used only by
# the coral="manta_1989_1991" option) leaves those sites at zero.
all(>=(0.0), baseline_total) || error("Negative initial coral cover")
DOMAIN_PATH == V2 && !all(>(0.0), baseline_total) && error("Cannot rescale zero initial coral cover")
covered = baseline_total .> 0.0
interpolated_cover = Float64.(surface.coral_manta_9m_fraction_1989_1991_idw)
new_cover = baseline_cover .* reshape(
    ifelse.(covered, interpolated_cover ./ max.(baseline_total, eps()), 0.0), 1, :)
all(isapprox.(vec(sum(new_cover; dims=1))[covered], interpolated_cover[covered]; atol=1e-12)) ||
    error("Interpolated coral totals not preserved")

# Observed 1993-95 fast share of hard coral per calibration reef (LTMP photo
# transects); other reefs take the mean of the observed reefs. Total cover and the
# structure within the fast (groups 1-3) and slow (4-5) blocks are preserved.
function observed_composition_cover(site_share::Vector{Float64})
    labels = String.(collect(domain.init_coral_cover.species))
    fast_rows = [any(startswith(label, g * "_") for g in FAST_GROUP_NAMES) for label in labels]
    count(fast_rows) == 21 || error("Unexpected fast-group rows in initial cover")
    cover = copy(baseline_cover)
    unchanged = 0
    for site in axes(cover, 2)
        total = baseline_total[site]
        total > 0.0 || continue
        fast = sum(baseline_cover[fast_rows, site])
        slow = total - fast
        if fast <= 0.0 || slow <= 0.0
            unchanged += 1
            continue
        end
        cover[fast_rows, site] .*= site_share[site] * total / fast
        cover[.!fast_rows, site] .*= (1.0 - site_share[site]) * total / slow
    end
    all(isapprox.(vec(sum(cover; dims=1)), baseline_total; atol=1e-12)) ||
        error("Observed composition changed total cover")
    return cover, unchanged
end
composition = if COMPOSITION || GRACE || CORAL || JOINT || MELLIN || BLEACH || REBASE
    photo = CSV.read(PHOTO_TRANSECT_CSV, DataFrame)
    photo = photo[(photo.purpose .== "COMPOSITION") .& (photo.variable .== "HARD CORAL") .&
        (photo.depth .== 9) .& in.(photo.report_year, Ref(PHOTO_YEARS)) .&
        in.(photo.reef_name, Ref(Set(values(PHOTO_REEFS)))), :]
    levels = GRACE || REBASE_COMPOSITION ? OTHER_WEIGHTS : Dict("acrpoc" => 0.0)
    out = Dict{String,Any}()
    for (level, weight) in levels
        yearly = combine(groupby(photo, [:reef_name, :report_year]),
            [:reefpage_category, :median] => ((c, v) ->
                (sum(v[in.(c, Ref(PHOTO_FAST_TAXA))]) + weight * sum(v[c .== "Other"])) / sum(v)) => :share)
        reef_share = Dict(reef => sum(yearly.share[yearly.reef_name .== obs]) /
            count(==(obs), yearly.reef_name) for (reef, obs) in PHOTO_REEFS)
        all(0.0 < v < 1.0 for v in values(reef_share)) || error("Invalid observed fast share")
        default_share = sum(values(reef_share)) / length(reef_share)
        site_share = fill(default_share, length(domain.loc_ids))
        for (reef, share) in reef_share
            site_share[mapping[reef].indices] .= share
        end
        cover, unchanged = observed_composition_cover(site_share)
        out[level] = (cover=cover, reef_share=reef_share, default_share=default_share,
            unchanged_sites=unchanged, other_weight=weight, site_share=site_share)
    end
    (levels=out, cover=out["acrpoc"].cover, reef_share=out["acrpoc"].reef_share,
     default_share=out["acrpoc"].default_share, unchanged_sites=out["acrpoc"].unchanged_sites,
     years=collect(PHOTO_YEARS))
else
    nothing
end
# Bleaching-floor screen: starting fast share = observed acrpoc share x m (capped), totals
# and within-block structure preserved; m = 1 reproduces composition.cover exactly.
# Lizard hard-coral initialisation: acrpoc cover with every Lizard site (all sub-reefs)
# multiplied by c_L; composition share preserved; c_L = 1 is composition.cover exactly.
lizcoral_cover = REBASE_LIZCORAL ? Dict(c => let cover = copy(composition.cover)
        cover[:, mapping["Lizard Island Reef"].indices] .*= c
        maximum(vec(sum(cover; dims=1))) <= 0.9 || error("Lizard-scaled cover exceeds 0.9 of habitable area")
        cover
    end for c in LIZCORAL_SCALES) : nothing
REBASE_LIZCORAL && lizcoral_cover[1.0] != composition.cover && error("c_L = 1 does not reproduce the acrpoc composition")
# Starting fast-share levels: reef-specific values at each level (other reefs take the default).
share_levels = REBASE_SHARE ? let raw = parse_share_levels(SHARE_LEVELS_JSON)
    function level_share(reef_level)
        share = fill(raw["default"][reef_level["default"]], length(domain.loc_ids))
        for reef in ("Lizard Island Reef", "MacGillivray Reef", "North Direction Reef")
            share[mapping[reef].indices] .= raw[reef][reef_level[reef]]
        end
        return share
    end
    picks = Dict{String,String}()
    if !isempty(SHARE_PICKS)
        for item in split(SHARE_PICKS, ";")
            reef, level = split(item, "=")
            level in SHARE_LEVEL_NAMES || error("Unknown share level $level")
            picks[String(reef)] = String(level)
        end
        Set(keys(picks)) == Set(["Lizard Island Reef", "MacGillivray Reef", "North Direction Reef"]) ||
            error("OWEN_INITIAL_SHARE_PICKS must name all three calibration reefs")
    end
    uniform(level) = Dict("default" => level, "Lizard Island Reef" => level,
        "MacGillivray Reef" => level, "North Direction Reef" => level)
    covers = Dict{String,Matrix{Float64}}(level => observed_composition_cover(level_share(uniform(level)))[1]
        for level in SHARE_LEVEL_NAMES)
    isempty(picks) || (covers["pick"] = observed_composition_cover(level_share(merge(Dict("default" => "q50"), picks)))[1])
    if TIMING || CANDIDATES || CULLING || CORALMIX
        for (set, set_picks) in TIMING_SHARE_SETS
            covers[set] = observed_composition_cover(level_share(merge(Dict("default" => "q50"), set_picks)))[1]
        end
        picks = TIMING_SHARE_SETS["pick"]
    end
    (raw=raw, covers=covers, picks=picks)
end : nothing
share_scaled = REBASE_FLOOR ? Dict(m => observed_composition_cover(
        min.(composition.levels["acrpoc"].site_share .* m, SHARE_CAP))[1] for m in SHARE_MULTIPLIERS) : nothing
REBASE_FLOOR && share_scaled[1.0] != composition.cover && error("m = 1 does not reproduce the acrpoc composition")

# Plain label vector captured once. Rebuilding from `collect(cube.species)` returns a
# dimension-wrapped vector, so each rebuild nested the lookup one level deeper and the
# DataCube construction time grew geometrically across arms (0.5 s -> >60 s after ~12
# arms; diagnosed 2026-10-13, runs 20261013T_06-08).
const COVER_SPECIES = String.(Array(ADRIA.YAXArrays.DD.lookup(domain.init_coral_cover, :species)))
const COVER_LOCATIONS = String.(collect(domain.loc_ids))
function set_initial_cover!(values)
    domain.init_coral_cover = ADRIA.DataCube(values;
        species=COVER_SPECIES, locations=COVER_LOCATIONS)
end

"""Adults/ha from CPUE: linear C / alpha when a_half = 0, otherwise the inverse of the
density-dependent manta model C = alpha * A^2 / (A + a_half)."""
invert_cpue(cpue, alpha, a_half) = a_half == 0.0 ? cpue ./ alpha :
    (cpue ./ alpha .+ sqrt.((cpue ./ alpha) .^ 2 .+ 4.0 .* a_half .* cpue ./ alpha)) ./ 2.0

function make_initial_state(alpha::Float64, name::String; cohort_multiplier::Float64=0.0,
    density_multiplier::Float64=1.0, cpue::Vector{Float64}=Float64.(surface.cots_cpue_1989_1991_idw),
    a_half::Float64=0.0)::String
    # alpha is an explicitly provisional COTS/tow per adult COTS/ha conversion.
    # density_multiplier is used only by the unit-migration compatibility arm.
    adults = a_half == 0.0 ? density_multiplier .* cpue ./ alpha :
        density_multiplier .* invert_cpue(cpue, alpha, a_half)
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
] : EXPERIMENT == "ulp_noise_floor" ? [
    # Numerical diagnostic: archived control with a_F moved by one Float64 ulp.
    (name="ulp_perturbed_control", alpha=0.015, coral="inherited_2023", initial="cohort_full"),
] : CORALMIX ? [
    [(name="wc_b150_q20_eb_lin_obs015_s1_m$(lpad(round(Int, 1000 * mu), 3, '0'))_l$(pct3(low))_cv2_pnd",
        alpha=0.015, coral="observed_composition", initial="cohort_full")
     for mu in CORALMIX_MU for low in CORALMIX_LOW];
    [(name="wc_b150_q20_eb_lin_obs015_s1_ref", alpha=0.015, coral="observed_composition", initial="cohort_full")]
] : CULLING ? [
    (name="wc_b150_q20_eb_lin_obs015_s1_$(cull)_cv2_pnd", alpha=0.015,
        coral="observed_composition", initial="cohort_full")
    for cull in CULL_ARMS
] : CANDIDATES ? [
    (name="wc_b150_q20_eb_lin_obs$(lpad(round(Int, 1000 * q), 3, '0'))_s$(k)_$(suffix)", alpha=q,
        coral="observed_composition", initial="cohort_full")
    for (suffix, q) in CANDIDATE_ARMS for k in eachindex(CANDIDATE_SEEDS)
] : TIMING ? [
    (name="wc_b150_q$(Q)_eb_lin_obs$(lpad(round(Int, 1000 * q), 3, '0'))_cv2_$(set)", alpha=q,
        coral="observed_composition", initial="cohort_full")
    for q in TIMING_OBS_Q for Q in TIMING_Q for set in ("pick", "pnd")
] : REBASE_SHARE ? [
    (name="wc_b150_q20_eb_lin_ref", alpha=0.015, coral="observed_composition", initial="cohort_full");
    (name="wc_b150_q20_eb_lin_cv2_acrpoc", alpha=0.015, coral="observed_composition", initial="cohort_full");
    [(name="wc_b150_q20_eb_lin_cv2_$(level)", alpha=0.015, coral="observed_composition", initial="cohort_full")
     for level in SHARE_LEVEL_NAMES];
    [(name="wc_b150_q20_eb_lin_cv2_pick", alpha=0.015, coral="observed_composition", initial="cohort_full")
     for _ in 1:(isempty(SHARE_PICKS) ? 0 : 1)]
] : REBASE_LIZCORAL ? [
    (name="wc_b150_q20_eb_lin_obs$(lpad(round(Int, 1000 * q), 3, '0'))_mu$(lpad(round(Int, 1000 * mu), 3, '0'))_c$(pct3(c))",
        alpha=q, coral="observed_composition", initial="cohort_full")
    for q in LIZCORAL_Q for mu in LIZCORAL_MU for c in LIZCORAL_SCALES
] : REBASE_FLOOR ? [
    (name="wc_b$(pct3(beta))_q20_eb_lin_f$(pct3(floor))_m$(pct3(m))", alpha=0.015,
        coral="observed_composition", initial="cohort_full")
    for floor in FLOOR_LEVELS for beta in REBASE_BETA for m in SHARE_MULTIPLIERS
] : REBASE_RETENTION ? [
    (name="wc_b150_q20_eb_lin_$(level)", alpha=0.015, coral="observed_composition",
        initial="cohort_full")
    for level in ("owen", "liz25", "liz50", "liz90", "all50")
] : REBASE_INITIAL ? [
    (name="wc_b150_q$(q)_$(source)_$(inversion)", alpha=0.015, coral="observed_composition",
        initial="cohort_full")
    for q in (20, 30) for source in ("idw", "eb") for inversion in ("lin", "h10")
] : REBASE_COMPOSITION ? [
    (name="wc_b$(pct3(beta))_q20_$(level)", alpha=0.015, coral="observed_composition",
        initial="cohort_full")
    for beta in REBASE_BETA for level in ("acrpoc", "halfother", "other")
] : REBASE ? [
    (name="ref_q30", alpha=0.015, coral="observed_composition", initial="cohort_full");
    (name="k500_q30", alpha=0.015, coral="observed_composition", initial="cohort_full");
    (name="mel3_q30", alpha=0.015, coral="observed_composition", initial="cohort_full");
    [(name="wc_b$(pct3(beta))_q$(q)", alpha=0.015, coral="observed_composition", initial="cohort_full")
     for beta in REBASE_BETA for q in REBASE_Q]
] : REFIT2 ? [
    (name="refit2_s100_mu030_rep", alpha=0.015, coral="observed_composition", initial="no_cots");
    [(name="refit2_s$(pct3(s))_mu$(lpad(round(Int, 1000 * mu), 3, '0'))_$(align)", alpha=0.015,
        coral="observed_composition", initial="no_cots")
     for s in REFIT2_SCALES for mu in REFIT2_MU for align in ("rep", "srv")]
] : REFIT ? [
    (name=refit_arm_name(floor, mu, q), alpha=0.015, coral="observed_composition", initial="no_cots")
    for floor in REFIT_FLOORS for mu in REFIT_MU for q in REFIT_Q
] : BLEACH ? [
    (name=bleach_arm_name(level, kappa, beta), alpha=0.015, coral="observed_composition",
        initial="no_cots")
    for (level, _) in BLEACH_BACKGROUND for kappa in BLEACH_KAPPA for beta in BLEACH_BETA
] : MELLIN ? [
    (name=mellin_arm_name(slow, fast, kappa), alpha=0.015, coral="observed_composition",
        initial="no_cots")
    for slow in JOINT_SLOW for fast in JOINT_FAST for kappa in JOINT_KAPPA
] : JOINT ? [
    (name=joint_arm_name(slow, fast, kappa, mu), alpha=0.015, coral="observed_composition",
        initial="no_cots")
    for slow in JOINT_SLOW for fast in JOINT_FAST for kappa in JOINT_KAPPA for mu in JOINT_BACKGROUND
] : CORAL ? [
    (name="dynamic_ref", alpha=0.015, coral="observed_composition", initial="cohort_full");
    [(name=coral_arm_name(slow, fast), alpha=0.015, coral="observed_composition", initial="no_cots")
     for slow in SLOW_MULTIPLIERS for fast in FAST_MULTIPLIERS]
] : GRACE ? [
    (name="$(level)_g$(grace)", alpha=0.015, coral="observed_composition", initial="cohort_full")
    for level in ("acrpoc", "halfother", "other") for grace in (0, 1)
] : COMPOSITION ? [
    (name="inh_m110_q30", alpha=0.015, coral="inherited_2023", initial="cohort_full");
    [(name="obs_m110_q$(q)", alpha=0.015, coral="observed_composition", initial="cohort_full")
     for q in (20, 30, 45, 70, 100)];
    [(name="obs_m$(m)_q30", alpha=0.015, coral="observed_composition", initial="cohort_full")
     for m in (100, 120)]
] : GROWTH ? [
    (name="m$(lpad(round(Int, 100 * multiplier), 3, '0'))_$(starvation)", alpha=0.015,
        coral="inherited_2023", initial="cohort_full")
    for multiplier in GROWTH_MULTIPLIERS for starvation in ("ref", "M0", "PG")
] : COLLAPSE ? [
    (name="$(coral)_$(collapse)", alpha=0.015, coral="inherited_2023", initial="cohort_full")
    for coral in ("g1", "grme", "grmefast") for collapse in ("none", "D2", "D2J", "D4")
] : PIT ? [
    (name="ext5_A", alpha=0.015, coral="inherited_2023", initial="cohort_full");
    [(name="ext5_PA_q$(capacity)", alpha=0.015, coral="inherited_2023", initial="cohort_full")
     for capacity in PIT_CAPACITIES]
] : INITIATION ? [
    (name=string(init, "_", mech), alpha=0.015, coral="inherited_2023", initial="cohort_full")
    for init in ("none", "ext5", "juv5", "both") for mech in ("none", "A", "PA")
] : TROUGH ? [
    (name=string(base, "_", treatment), alpha=0.015, coral="inherited_2023", initial="cohort_full")
    for base in ("b05", "bgeo") for treatment in ("none", "P", "D", "A", "PDA")
] : BENTHIC_CYCLONE ? [
    (name="k05_nocyc", alpha=0.015, coral="inherited_2023", initial="cohort_full"),
    (name="k05_cyc", alpha=0.015, coral="inherited_2023", initial="cohort_full"),
    (name="ktotal_nocyc", alpha=0.015, coral="inherited_2023", initial="cohort_full"),
    (name="ktotal_cyc", alpha=0.015, coral="inherited_2023", initial="cohort_full"),
    (name="kmapped_nocyc", alpha=0.015, coral="inherited_2023", initial="cohort_full"),
    (name="kmapped_cyc", alpha=0.015, coral="inherited_2023", initial="cohort_full"),
] : HABITAT ? [
    (name="k05_bleaching", alpha=0.015, coral="inherited_2023", initial="cohort_full"),
    (name="krme_bleaching", alpha=0.015, coral="inherited_2023", initial="cohort_full"),
    (name="k05_no_bleaching", alpha=0.015, coral="inherited_2023", initial="cohort_full"),
    (name="krme_no_bleaching", alpha=0.015, coral="inherited_2023", initial="cohort_full"),
] : SWITCHING_REPAIR ? [
    (name="condition_ref", alpha=0.015, coral="inherited_2023", initial="cohort_full"),
    (name="switching_ref", alpha=0.015, coral="inherited_2023", initial="cohort_full"),
    (name="switch_diet_starvation", alpha=0.015, coral="inherited_2023", initial="cohort_full"),
    (name="switch_low_efficiency", alpha=0.015, coral="inherited_2023", initial="cohort_full"),
    (name="switch_diet_starvation_low_efficiency", alpha=0.015, coral="inherited_2023", initial="cohort_full"),
] : SWITCHING_CONDITION ? [
    (name="base", alpha=0.015, coral="inherited_2023", initial="cohort_full"),
    (name="switching", alpha=0.015, coral="inherited_2023", initial="cohort_full"),
    (name="condition", alpha=0.015, coral="inherited_2023", initial="cohort_full"),
    (name="switching_condition", alpha=0.015, coral="inherited_2023", initial="cohort_full"),
] : SECOND_WAVE_AUDIT ? [
    (name="migrated_control", alpha=0.015, coral="inherited_2023", initial="cohort_full"),
    (name="per_capita_structural", alpha=0.015, coral="inherited_2023", initial="cohort_full"),
] : EXPERIMENT == "unit_compat_k2" ? [
    (name="compat_rescaled_exact", alpha=0.015, coral="inherited_2023", initial="cohort_full"),
] : MIGRATION ? [
    (name="compat_rescaled_exact", alpha=0.015, coral="inherited_2023", initial="cohort_full"),
    (name="migrated_control", alpha=0.015, coral="inherited_2023", initial="cohort_full"),
    (name="migrated_structural", alpha=0.015, coral="inherited_2023", initial="cohort_full"),
    (name="per_capita_legacy", alpha=0.015, coral="inherited_2023", initial="cohort_full"),
    (name="per_capita_structural", alpha=0.015, coral="inherited_2023", initial="cohort_full"),
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
) : EXPERIMENT == "ulp_noise_floor" ? Dict(
    "ulp_perturbed_control" => make_initial_state(0.015, "alpha015_cohort_full"; cohort_multiplier=1.0)
) : TIMING || CANDIDATES || CULLING || CORALMIX ? let
    states = Dict(q => make_initial_state(q, "alpha$(lpad(round(Int, 1000 * q), 3, '0'))_cohort_full_eb_lin";
            cohort_multiplier=1.0, cpue=eb_surface) for q in TIMING_OBS_Q)
    Dict(t.name => states[timing_obs_q(t.name)] for t in treatments)
end : REBASE_SHARE ? let
    state = make_initial_state(0.015, "alpha015_cohort_full_eb_lin"; cohort_multiplier=1.0, cpue=eb_surface)
    Dict(t.name => state for t in treatments)
end : REBASE_LIZCORAL ? let
    # eb seeds inverted linearly with the arm's observation q; q = 0.015 reproduces run 17's eb_lin state.
    states = Dict(q => make_initial_state(q, "alpha$(lpad(round(Int, 1000 * q), 3, '0'))_cohort_full_eb_lin";
            cohort_multiplier=1.0, cpue=eb_surface) for q in LIZCORAL_Q)
    Dict(t.name => states[lizcoral_parts(t.name).q] for t in treatments)
end : REBASE_INITIAL || REBASE_RETENTION || REBASE_FLOOR ? let
    sources = Dict("idw" => Float64.(surface.cots_cpue_1989_1991_idw), "eb" => eb_surface)
    halves = Dict("lin" => 0.0, "h10" => INITIAL_A_HALF)
    states = Dict((source, inversion) => make_initial_state(0.015,
            "alpha015_cohort_full_$(source)_$(inversion)"; cohort_multiplier=1.0,
            cpue=sources[source], a_half=halves[inversion])
        for source in keys(sources) for inversion in keys(halves))
    Dict(t.name => states[(split(t.name, "_")[4], split(t.name, "_")[5])] for t in treatments)
end : PIT || COLLAPSE || GROWTH || COMPOSITION || GRACE || CORAL || JOINT || MELLIN || BLEACH || REBASE ? let
    state_path = make_initial_state(0.015, "alpha015_cohort_full"; cohort_multiplier=1.0)
    Dict(t.name => state_path for t in treatments)
end : INITIATION ? let
    state_path = make_initial_state(0.015, "alpha015_cohort_full"; cohort_multiplier=1.0)
    juvenile_path = make_initial_state(0.015, "alpha015_cohort_juvenile_x5";
        cohort_multiplier=JUVENILE_PRELOAD_MULTIPLIER)
    Dict(t.name => (startswith(t.name, "juv5") || startswith(t.name, "both") ?
        juvenile_path : state_path) for t in treatments)
end : MIGRATION ? let
    state_path = make_initial_state(0.015, "alpha015_cohort_full"; cohort_multiplier=1.0)
    compat_path = make_initial_state(0.015, "alpha015_cohort_full_times_k";
        cohort_multiplier=1.0, density_multiplier=K_MIGRATION)
    Dict(t.name => (t.name == "compat_rescaled_exact" ? compat_path : state_path)
        for t in treatments)
end : EXPERIMENT in ("crash_factorial", "lagged_hazard", "cohort_audit", "structural_crash") ? let
    state_path = make_initial_state(0.015, "alpha015_cohort_full"; cohort_multiplier=1.0)
    Dict(t.name => state_path for t in treatments)
end : Dict{String,String}()
archived_m3 = Float64(scenario[1, :m3])
archived_density_params = Dict(name => Float64(scenario[1, name])
    for name in (:a_F, :a_S, :IMM, :b_ricker))
# ADRIA otherwise seeds its coral RNG with floor(sum(parameter values))
# (utils/setup.jl), so rescaling parameter values silently changes the coral
# random stream. Migration screens pin the archived control's derived seed.
archived_rng_seed = floor(Int64, sum(Float64(scenario[1, name])
    for name in names(scenario) if name != "RCP"))
(CANDIDATES || CULLING || CORALMIX) && CANDIDATE_SEEDS[1] != (SEED, archived_rng_seed) &&
    error("Candidate seed pair 1 must be the archived pair ($(SEED), $(archived_rng_seed))")
# Coral calibration groups. The Lizard loader assigns group 1 to every site; the
# GBR-wide RME layer assigns each reef a group. ADRIA keeps only the CB-group
# factors present in the domain and indexes them by group number, so RME groups
# are relabelled to contiguous IDs that carry each RME group's default factors.
base_cb_groups = copy(domain.loc_data.CB_CALIB_GROUPS)
rme_cb = if COLLAPSE || GROWTH || COMPOSITION || GRACE || CORAL || JOINT || MELLIN || BLEACH || REBASE
    layer = ADRIA.GDF.read(RME_GBR_LAYER)
    lookup = Dict(String(row.GBRMPA_ID) => Int(row.CB_CALIB_GROUPS) for row in eachrow(layer))
    reef_ids = String.(domain.loc_data.GBRMPA_ID)
    all(haskey(lookup, id) for id in reef_ids) || error("RME layer lacks a domain reef")
    site_groups = [lookup[id] for id in reef_ids]
    rme_ids = sort(unique(site_groups))
    label_of = Dict(g => i for (i, g) in enumerate(rme_ids))
    (rme_ids=rme_ids, site_groups=site_groups, labels=[label_of[g] for g in site_groups],
     reef_groups=Dict(id => lookup[id] for id in unique(reef_ids)))
else
    nothing
end
function coral_group_scenario(base::DataFrame, fast_multiplier::Float64;
    slow_multiplier::Float64=1.0)::DataFrame
    rng_state = copy(Random.default_rng())
    domain.loc_data.CB_CALIB_GROUPS .= rme_cb.labels
    grouped = ADRIA.sample(domain, 2)[1:1, :]
    domain.loc_data.CB_CALIB_GROUPS .= base_cb_groups
    copy!(Random.default_rng(), rng_state)
    for name in names(grouped)
        occursin("cb_group_", name) && continue
        hasproperty(base, Symbol(name)) || error("Missing scenario factor $name")
        grouped[!, name] .= base[1, name]
    end
    extension = ADRIA.linear_extension_group_scale_factors()
    mortality = ADRIA.mb_rate_group_scale_factors()
    acceleration = (steepness=ADRIA.growth_acceleration_steepness(),
        midpoint=ADRIA.growth_acceleration_midpoint(),
        height=ADRIA.growth_acceleration_height())
    for (label, group) in enumerate(rme_cb.rme_ids)
        for (index, fg) in enumerate(string.(ADRIA.functional_group_names()))
            # Fast prey (COTS prey map) are functional groups 1-3.
            multiplier = index <= 3 ? fast_multiplier : slow_multiplier
            grouped[!, "linear_extension_scale_cb_group_$(label)_$(fg)"] .=
                extension[index, group] * multiplier
            grouped[!, "mb_rate_scale_cb_group_$(label)_$(fg)"] .= mortality[index, group]
        end
        for parameter in (:steepness, :midpoint, :height)
            grouped[!, "growth_acceleration_cb_group_$(label)_$(parameter)"] .=
                getfield(acceleration, parameter)[group]
        end
    end
    return grouped
end
scenario_rme = COLLAPSE ? coral_group_scenario(scenario, 1.0) : nothing
scenario_rme_fast = COLLAPSE ? coral_group_scenario(scenario, FAST_LINEAR_EXTENSION_MULTIPLIER) : nothing
growth_scenarios = GROWTH ?
    Dict(m => coral_group_scenario(scenario, m) for m in GROWTH_MULTIPLIERS) :
    COMPOSITION ? Dict(m => coral_group_scenario(scenario, m) for m in COMPOSITION_MULTIPLIERS) :
    GRACE || REBASE ? Dict(1.0 => coral_group_scenario(scenario, 1.0)) :
    BLEACH ? Dict("100_100" => coral_group_scenario(scenario, 1.0; slow_multiplier=1.0)) :
    JOINT || MELLIN ? Dict("$(pct3(slow))_$(pct3(fast))" => coral_group_scenario(scenario, fast;
        slow_multiplier=slow) for slow in JOINT_SLOW for fast in JOINT_FAST) :
    CORAL ? merge(Dict("dynamic_ref" => coral_group_scenario(scenario, 1.0)),
        Dict(coral_arm_name(slow, fast) => coral_group_scenario(scenario, fast;
            slow_multiplier=slow) for slow in SLOW_MULTIPLIERS for fast in FAST_MULTIPLIERS)) :
    nothing
raw_obs = CSV.read(OBS, DataFrame)
scores = DataFrame(treatment=String[], reef_name=String[], observation_treatment=String[],
    loss=Float64[], sim_peak_count=Int[], sim_peak_years=String[],
    sim_peak_heights_cpue=String[], matched_peaks=Int[])
trajectories = DataFrame(treatment=String[], reef_name=String[], year=Int[],
    forcing_year=Int[], adults_ha=Float64[], simulated_cpue=Float64[], coral_cover=Float64[])
flows = DataFrame(treatment=String[], reef_name=String[], year=Int[])
# Observed culling: requested and applied removals (COTS counts = density x site area) per
# calibration reef and for the whole domain, by phase and stage.
culls = DataFrame(treatment=String[], scope=String[], year=Int[],
    requested_pre_juveniles=Float64[], requested_pre_adults=Float64[],
    requested_post_juveniles=Float64[], requested_post_adults=Float64[],
    removed_pre_juveniles=Float64[], removed_pre_adults=Float64[],
    removed_post_juveniles=Float64[], removed_post_adults=Float64[])
cull_site_area_ha = if CULLING
    site_table = CSV.read(joinpath(V3, "site_to_reef.csv"), DataFrame; types=Dict(:site_id => String))
    area = Dict(r.site_id => r.area_m2 / 1e4 for r in eachrow(site_table))
    Float64[area[String(id)] for id in domain.loc_ids]
else
    Float64[]
end
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
    senescent_adults_ha=Float64[], senescent_deaths_ha=Float64[],
    preferred_diet_share=Float64[], juvenile_predation_ha=Float64[],
    maturation_fraction=Float64[], body_condition=Float64[],
    total_coral_site_fraction=Float64[], fast_coral_site_fraction=Float64[],
    outbreak_clock=Float64[], collapse_removed_ha=Float64[], low_food_years=Float64[])
site_audit = DataFrame(treatment=String[], reef_name=String[], site_index=Int[],
    reef_weight=Float64[], habitable_fraction=Float64[], year=Int[],
    recruits_ha=Float64[], subadults_ha=Float64[], adults_ha=Float64[],
    body_condition=Float64[], preferred_food_memory=Float64[], food_survival=Float64[],
    fast_coral=Float64[], total_coral=Float64[], local_production_ha=Float64[],
    settled_recruits_ha=Float64[])
initial_reef = DataFrame(treatment=String[], reef_name=String[],
    seed_adults_ha=Float64[], seed_corals_fraction=Float64[],
    idw_cpue=Float64[], idw_manta_coral_fraction=Float64[])
failures = DataFrame(treatment=String[], error=String[])
# Optional arm filter (regex on treatment name) to complete an interrupted grid in a new
# run directory; recorded in metadata. Empty runs every arm.
const ARM_FILTER = get(ENV, "OWEN_1991_ARM_FILTER", "")
run_treatments = isempty(ARM_FILTER) ? treatments :
    filter(t -> occursin(Regex(ARM_FILTER), t.name), treatments)
isempty(run_treatments) && error("OWEN_1991_ARM_FILTER matches no treatment")

for t in run_treatments
    arm_started = time()
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
    compat = MIGRATION && t.name == "compat_rescaled_exact"
    k = MIGRATION ? K_MIGRATION : 1.0
    migration_structural = MIGRATION &&
        (endswith(t.name, "_structural") || SWITCHING_CONDITION || SWITCHING_REPAIR || HABITAT ||
         BENTHIC_CYCLONE || TROUGH || V3RUN)
    ENV["COTS_PREFERRED_PREY_STARVATION"] = (structural &&
        occursin("preferred_starvation", t.name)) || migration_structural ? "true" : "false"
    ENV["COTS_STARVATION_PREFERRED_THRESHOLD"] = string(STRUCTURAL_THRESHOLD)
    ENV["COTS_STARVATION_MEMORY_YEARS"] = GROWTH && endswith(t.name, "_M0") ?
        "0.0" : string(STRUCTURAL_MEMORY_YEARS)
    ENV["COTS_STARVATION_POST_GRAZING"] = GROWTH && endswith(t.name, "_PG") ? "true" : "false"
    ENV["COTS_STARVATION_GRACE_YEARS"] = GRACE ? split(t.name, "_g")[end] : "0"
    if ((CORAL || JOINT) && startswith(t.name, "forced")) || MELLIN || BLEACH
        ENV["COTS_FORCED_ADULTS_CSV"] = REFIT2 ? refit_forcing[split(t.name, "_")[4]] : REFIT ?
            refit_forcing[parse(Int, split(t.name, "_q")[end]) / 1000] : FORCED_ADULTS_CSV
    else
        delete!(ENV, "COTS_FORCED_ADULTS_CSV")
    end
    ENV["COTS_ADULT_SENESCENCE"] = (structural &&
        occursin("senescence", t.name)) || migration_structural ? "true" : "false"
    ENV["COTS_SENESCENCE_AGE"] = string(STRUCTURAL_SENESCENCE_AGE)
    ENV["COTS_SENESCENCE_MORTALITY"] = string(STRUCTURAL_SENESCENCE_MORTALITY)
    ENV["COTS_PER_CAPITA_CONSUMPTION"] = MIGRATION &&
        (startswith(t.name, "per_capita") || SWITCHING_CONDITION || SWITCHING_REPAIR ||
         HABITAT || BENTHIC_CYCLONE || TROUGH || V3RUN) ? "true" : "false"
    ENV["COTS_PREY_SWITCHING"] = (SWITCHING_CONDITION && startswith(t.name, "switching")) ||
        (SWITCHING_REPAIR && t.name != "condition_ref") ? "true" : "false"
    ENV["COTS_DIET_LINKED_STARVATION"] = SWITCHING_REPAIR &&
        occursin("diet_starvation", t.name) ? "true" : "false"
    ENV["COTS_NONPREFERRED_FEEDING_EFFICIENCY"] = SWITCHING_REPAIR &&
        occursin("low_efficiency", t.name) ? string(NONPREFERRED_EFFICIENCY) : "1.0"
    ENV["COTS_SWITCH_MIDPOINT"] = string(SWITCH_MIDPOINT)
    ENV["COTS_SWITCH_STEEPNESS"] = string(SWITCH_STEEPNESS)
    ENV["COTS_PREFERRED_CONDITION"] = (SWITCHING_CONDITION &&
        endswith(t.name, "condition")) || SWITCHING_REPAIR || HABITAT || BENTHIC_CYCLONE ||
        TROUGH || V3RUN ? "true" : "false"
    setup_marks = [("env", time())]
    # Habitat screen: run-time habitable fraction and bleaching counterfactual.
    domain.loc_data.k .= HABITAT && startswith(t.name, "krme") ? rme_k :
        BENTHIC_CYCLONE && startswith(t.name, "ktotal") ? benthic.total :
        BENTHIC_CYCLONE && startswith(t.name, "kmapped") ? benthic.mapped :
        TROUGH && startswith(t.name, "bgeo") ? geomorphic_k.k : baseline_k
    domain.cyclone_mortality_scens = V3RUN ? domain_cyclones :
        (BENTHIC_CYCLONE && endswith(t.name, "_cyc")) ||
        (TROUGH && startswith(t.name, "bgeo")) ? cyclone_cube : zero_cyclones
    if JOINT
        slow_value, fast_value, kappa_value, mu_value = joint_arm_values(t.name)
        domain.cyclone_mortality_scens = joint_cyclones[(kappa_value, mu_value)]
    elseif MELLIN
        kappa_value = parse(Int, split(t.name, "_")[4][2:end]) / 100
        domain.cyclone_mortality_scens = mellin_cyclones[kappa_value]
    elseif REFIT
        domain.cyclone_mortality_scens = refit_cyclones[parse(Int, split(t.name, "_")[3][3:end]) / 1000]
    elseif BLEACH
        bleach_parts = split(t.name, "_")
        domain.cyclone_mortality_scens =
            bleach_cyclones[(bleach_parts[2], parse(Int, bleach_parts[3][2:end]) / 100)]
    elseif REBASE
        domain.cyclone_mortality_scens = rebase_cyclones[split(t.name, "_")[1]]
    end
    if REBASE_LIZCORAL
        domain.cyclone_mortality_scens = lizcoral_cyclones[lizcoral_parts(t.name).mu]
    end
    if REBASE_SHARE && occursin("_cv2_", t.name)
        domain.cyclone_mortality_scens = cv2_cyclones
    end
    if CORALMIX && occursin("_cv2_", t.name)
        domain.cyclone_mortality_scens = coralmix_cyclones[coralmix_parts(t.name).mu]
    end
    if (REFIT2 && split(t.name, "_")[2] != "s100") || (REBASE_SHARE && occursin("_cv2_", t.name))
        ENV["ADRIA_BLEACHING_LOW_DHW_SCALE"] = REBASE_SHARE ? string(CV2_LOW_DHW_SCALE) :
            string(parse(Int, split(t.name, "_")[2][2:end]) / 100)
        ENV["ADRIA_BLEACHING_LOW_DHW_THRESHOLD"] = "4.0"
    else
        delete!(ENV, "ADRIA_BLEACHING_LOW_DHW_SCALE")
        delete!(ENV, "ADRIA_BLEACHING_LOW_DHW_THRESHOLD")
    end
    if CORALMIX && occursin("_cv2_", t.name) && coralmix_parts(t.name).low == 1.0
        delete!(ENV, "ADRIA_BLEACHING_LOW_DHW_SCALE")
        delete!(ENV, "ADRIA_BLEACHING_LOW_DHW_THRESHOLD")
    elseif CORALMIX && occursin("_cv2_", t.name)
        ENV["ADRIA_BLEACHING_LOW_DHW_SCALE"] = string(coralmix_parts(t.name).low)
        ENV["ADRIA_BLEACHING_LOW_DHW_THRESHOLD"] = "4.0"
    end
    if REBASE_FLOOR && split(t.name, "_")[6] != "f100"
        ENV["ADRIA_BLEACHING_DHW_FLOOR"] = string(parse(Int, split(t.name, "_")[6][2:end]) / 100)
    elseif REFIT && !REFIT2 && split(t.name, "_")[2] != "f100"
        ENV["ADRIA_BLEACHING_DHW_FLOOR"] = string(parse(Int, split(t.name, "_")[2][2:end]) / 100)
    else
        delete!(ENV, "ADRIA_BLEACHING_DHW_FLOOR")
    end
    if REBASE_RETENTION
        level = split(t.name, "_")[6]
        for i in eachindex(base_reef_connectivity)
            domain.cots_forcing.reef_connectivity[i] = isnothing(RETENTION_LEVELS[level]) ?
                base_reef_connectivity[i] : retention_matrices[level][i]
        end
    end
    push!(setup_marks, ("k_cyclone", time()))
    trough_treatment = COLLAPSE || GROWTH || COMPOSITION || GRACE || CORAL || JOINT || MELLIN || BLEACH || REBASE ? "PA" : TROUGH || V3RUN ? split(t.name, "_")[2] : ""
    external_multiplier = COLLAPSE || GROWTH || COMPOSITION || GRACE || CORAL || JOINT || MELLIN || BLEACH || REBASE ||
        (V3RUN && (startswith(t.name, "ext5") || startswith(t.name, "both"))) ?
        EXTERNAL_MULTIPLIER : 1.0
    predation_capacity = COMPOSITION ? composition_arm_q(t.name) :
        REBASE ? parse(Float64, match(r"_q(\d+)", t.name).captures[1]) :
        COLLAPSE || GROWTH || GRACE || CORAL || JOINT || MELLIN || BLEACH ? COLLAPSE_Q : PIT && occursin("_q", t.name) ?
        parse(Float64, split(t.name, "_q")[end]) : JUVENILE_PREDATION_CAPACITY
    collapse_arm = COLLAPSE ? split(t.name, "_")[2] : "none"
    ENV["COTS_OUTBREAK_COLLAPSE"] = collapse_arm == "none" ? "false" : "true"
    ENV["COTS_COLLAPSE_THRESHOLD_HA"] = string(COLLAPSE_THRESHOLD_HA)
    ENV["COTS_COLLAPSE_BACKGROUND_HA"] = string(COLLAPSE_BACKGROUND_HA)
    ENV["COTS_COLLAPSE_DURATION_YEARS"] = startswith(collapse_arm, "D2") ? "2" : "4"
    ENV["COTS_COLLAPSE_JUVENILES"] = collapse_arm == "D2J" ? "true" : "false"
    domain.loc_data.CB_CALIB_GROUPS .= (COLLAPSE && startswith(t.name, "grme")) ||
        GROWTH || COMPOSITION || GRACE || CORAL || JOINT || MELLIN || BLEACH || REBASE ?
        rme_cb.labels : base_cb_groups
    push!(setup_marks, ("collapse_cbgroups", time()))
    ENV["COTS_JUVENILE_PREDATION"] = occursin("P", trough_treatment) ? "true" : "false"
    ENV["COTS_JUVENILE_PREDATION_CAPACITY"] = string(predation_capacity)
    ENV["COTS_PREFERRED_MATURATION"] = occursin("D", trough_treatment) ? "true" : "false"
    ENV["COTS_ALLEE_STRETCHED"] = occursin("A", trough_treatment) ? "true" : "false"
    ENV["COTS_ALLEE_EXPONENT"] = "0.5"
    if HABITAT
        ENV["ADRIA_DEBUG_BLEACHING_SCALAR"] = occursin("no_bleaching", t.name) ? "0.0" : "1.0"
    elseif REFIT
        ENV["ADRIA_DEBUG_BLEACHING_SCALAR"] = "1.5"
    elseif BLEACH
        ENV["ADRIA_DEBUG_BLEACHING_SCALAR"] = string(parse(Int, split(t.name, "_b")[end]) / 100)
    elseif REBASE && startswith(t.name, "wc_")
        ENV["ADRIA_DEBUG_BLEACHING_SCALAR"] = string(parse(Int, split(t.name, "_")[2][2:end]) / 100)
    else
        delete!(ENV, "ADRIA_DEBUG_BLEACHING_SCALAR")
    end
    ENV["COTS_CONDITION_SATURATION"] = string(CONDITION_SATURATION)
    ENV["COTS_CONDITION_MEMORY_YEARS"] = string(CONDITION_MEMORY_YEARS)
    ENV["COTS_CONSUMPTION_M2_PER_ADULT_YEAR"] = string(PER_CAPITA_M2_PER_ADULT_YEAR)
    # Exact density-unit reparameterisation; identity outside the migration screen.
    scenario[!, :a_F] .= EXPERIMENT == "ulp_noise_floor" ?
        nextfloat(archived_density_params[:a_F]) : archived_density_params[:a_F] / k
    scenario[!, :a_S] .= archived_density_params[:a_S] / k
    scenario[!, :b_ricker] .= archived_density_params[:b_ricker] / k
    scenario[!, :IMM] .= archived_density_params[:IMM] * k
    ENV["COTS_ALLEE_THRESHOLD"] = string(compat ? 3.0 * k : 3.0)
    if MIGRATION
        ENV["ADRIA_RNG_SEED"] = string(archived_rng_seed)
    else
        delete!(ENV, "ADRIA_RNG_SEED")
    end
    if CANDIDATES || CULLING || CORALMIX
        connectivity_seed, rng_seed = candidate_seed_pair(t.name)
        ENV["COTS_CONNECTIVITY_SEED"] = string(connectivity_seed)
        ENV["ADRIA_RNG_SEED"] = string(rng_seed)
    end
    if CULLING && !isnothing(cull_table(t.name))
        ENV["COTS_CULL_REMOVALS_CSV"] = cull_table(t.name)
    else
        delete!(ENV, "COTS_CULL_REMOVALS_CSV")
    end
    scenario[!, :m3] .= EXPERIMENT == "crash_factorial" &&
        occursin("mortality030", t.name) ? 0.30 : archived_m3
    push!(setup_marks, ("env_scenario_cols", time()))
    initial_cover = t.coral == "manta_1989_1991" ? new_cover :
        t.coral == "observed_composition" ?
        (GRACE ? composition.levels[split(t.name, "_")[1]].cover :
         REBASE_COMPOSITION ? composition.levels[split(t.name, "_")[end]].cover :
         REBASE_FLOOR ? share_scaled[parse(Int, split(t.name, "_")[7][2:end]) / 100] :
         REBASE_LIZCORAL ? lizcoral_cover[lizcoral_parts(t.name).c] :
         REBASE_SHARE && !(split(t.name, "_")[end] in ("ref", "acrpoc")) ? share_levels.covers[split(t.name, "_")[end]] :
         composition.cover) : baseline_cover
    if TROUGH && startswith(t.name, "bgeo")
        initial_cover = initial_cover .* reshape(Float64.(.!geomorphic_k.zero), 1, :)
    end
    if BENTHIC_CYCLONE && (startswith(t.name, "ktotal") || startswith(t.name, "kmapped"))
        # Sites the benthic map gives k = 0 (floored for ADRIA) start with no coral.
        zero_mask = startswith(t.name, "ktotal") ? benthic.total_zero : benthic.mapped_zero
        initial_cover = initial_cover .* reshape(Float64.(.!zero_mask), 1, :)
    end
    set_initial_cover!(initial_cover)
    push!(setup_marks, ("initial_cover", time()))
    if t.name == "inherited_grazing_off_2004_2012"
        ENV["COTS_DIAGNOSTIC_GRAZING_OFF_START_YEAR"] = "2004"
        ENV["COTS_DIAGNOSTIC_GRAZING_OFF_END_YEAR"] = "2012"
    else
        delete!(ENV, "COTS_DIAGNOSTIC_GRAZING_OFF_START_YEAR")
        delete!(ENV, "COTS_DIAGNOSTIC_GRAZING_OFF_END_YEAR")
    end
    ENV["ADRIA_COTS_ENABLED"] = t.initial == "no_cots" ? "false" : "true"
    ENV["COTS_LARVAL_FECUNDITY"] = t.initial == "no_cots" ? "0.0" : string(fecundity)
    ENV["COTS_EXTERNAL_OUTBREAK_PRODUCTION"] = t.initial == "no_cots" ? "0.0" :
        string(boundary * k * external_multiplier)
    ENV["COTS_IMMIGRATION_SCALAR"] = t.initial == "no_cots" ? "0.0" : "1.0"
    if t.initial in ("legacy", "no_cots")
        delete!(ENV, "COTS_INITIAL_STATE_CSV")
    else
        key = EXPERIMENT == "initialization" ?
            (t.alpha == 0.015 ? "idw_alpha015" : "idw_alpha150") : t.name
        ENV["COTS_INITIAL_STATE_CSV"] = initial_paths[key]
    end
    push!(setup_marks, ("env_tail", time()))
    setup_breakdown = join(["$(name) $(round(mark - (i == 1 ? arm_started : setup_marks[i - 1][2]); digits=1))"
        for (i, (name, mark)) in enumerate(setup_marks)], ", ")
    Random.seed!(SEED)
    try
        active_scenario = GROWTH ? growth_scenarios[growth_arm_multiplier(t.name)] :
            COMPOSITION ? growth_scenarios[parse(Int, split(t.name, "_")[2][2:end]) / 100] :
            GRACE || REBASE ? growth_scenarios[1.0] :
            CORAL ? growth_scenarios[t.name] :
            BLEACH ? growth_scenarios["100_100"] :
            JOINT || MELLIN ? growth_scenarios[join(split(t.name, "_")[2:3] .|> (x -> x[2:end]), "_")] :
            COLLAPSE && startswith(t.name, "grmefast") ? scenario_rme_fast :
            COLLAPSE && startswith(t.name, "grme") ? scenario_rme : scenario
        if active_scenario !== scenario
            for name in (:a_F, :a_S, :b_ricker, :IMM, :m3, :allee_threshold)
                active_scenario[!, name] .= scenario[1, name]
            end
        end
        setup_seconds = time() - arm_started
        scenario_timing = @timed ADRIA.run_scenario(domain, active_scenario[1, :])
        post_started = time()
        result = scenario_timing.value
        Int.(result.cots_calendar_year) == years || error("1991 calendar alignment failed")
        if t.initial == "no_cots"
            maximum(abs, result.cots_log) == 0.0 ||
                error("No-COTS diagnostic received nonzero COTS state")
        end
        coral = dropdims(sum(result.raw; dims=(2, 3)); dims=(2, 3))
        if SECOND_WAVE_AUDIT
            # Site-level readout for target-reef sites only (COTS ha^-1 of site
            # area; cover as fractions of habitable area).
            raw = Array(result.raw)
            habitable = Float64.(ADRIA.location_k(domain))
            for target in TARGETS
                m = mapping[target.model]
                for (weight, site) in zip(m.weights, m.indices), (idx, year) in enumerate(years)
                    push!(site_audit, (t.name, target.model, site, weight,
                        habitable[site], year,
                        result.cots_log[idx, 1, site], result.cots_log[idx, 2, site],
                        result.cots_log[idx, 3, site], result.cots_condition_log[idx, site],
                        result.cots_mechanism_log[idx, 1, site],
                        result.cots_mechanism_log[idx, 2, site],
                        sum(raw[idx, 1:3, :, site]), sum(raw[idx, :, :, site]),
                        result.cots_flow_log[idx, 1, site],
                        result.cots_flow_log[idx, 10, site]))
                end
            end
        end
        if CULLING
            removed = Array(result.cots_culled_log)
            requested = isnothing(cull_table(t.name)) ? zeros(size(removed)) :
                ADRIA.load_cots_cull_removals(cull_table(t.name), domain.loc_ids, years)
            size(requested) == size(removed) || error("Cull request and log dimensions differ")
            all(removed .<= requested .+ 1e-9) || error("Applied removals exceed requests")
            scopes = [("domain", collect(eachindex(domain.loc_ids)));
                      [(target.model, mapping[target.model].indices) for target in TARGETS]]
            for (scope, sites) in scopes, (idx, year) in enumerate(years)
                a = cull_site_area_ha[sites]
                push!(culls, (t.name, scope, year,
                    (sum(requested[idx, k, sites] .* a) for k in 1:4)...,
                    (sum(removed[idx, k, sites] .* a) for k in 1:4)...))
            end
        end
        fast_coral = dropdims(sum(Array(result.raw)[:, 1:3, :, :]; dims=(2, 3)); dims=(2, 3))
        site_k = reshape(Float64.(domain.loc_data.k), 1, :)
        coral_site = Matrix(coral) .* site_k
        fast_site = fast_coral .* site_k
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
                REBASE_INITIAL || REBASE_RETENTION || REBASE_FLOOR || REBASE_LIZCORAL || REBASE_SHARE ? Float64.(CSV.read(initial_paths[t.name], DataFrame).adult_ha) :
                Float64.(surface.cots_cpue_1989_1991_idw) ./ t.alpha .* (compat ? k : 1.0)
            push!(initial_reef, (t.name, reef,
                sum(m.weights .* state_seed[m.indices]),
                sum(m.weights .* vec(sum(Array(domain.init_coral_cover); dims=1))[m.indices]),
                sum(m.weights .* Float64.(surface.cots_cpue_1989_1991_idw[m.indices])),
                sum(m.weights .* interpolated_cover[m.indices])))
            reef_flows = [aggregate_reef_density(Matrix(result.cots_flow_log[:, channel, :]), m)
                for channel in 1:length(FLOW_NAMES)]
            reef_stages = EXPERIMENT in ("cohort_audit", "structural_crash") || MIGRATION ?
                [aggregate_reef_density(Matrix(result.cots_log[:, stage, :]), m)
                    for stage in 1:3] : Vector{Float64}[]
            reef_hazards = EXPERIMENT == "lagged_hazard" ?
                [aggregate_reef_density(Matrix(result.cots_hazard_log[:, channel, :]), m)
                    for channel in 1:3] : Vector{Float64}[]
            reef_mechanisms = EXPERIMENT == "structural_crash" || MIGRATION ?
                vcat([aggregate_reef_density(Matrix(result.cots_mechanism_log[:, channel, :]), m)
                    for channel in 1:7],
                    [aggregate_reef_density(Matrix(result.cots_condition_log), m)]) :
                Vector{Float64}[]
            reef_collapse = EXPERIMENT == "structural_crash" || MIGRATION ?
                [aggregate_reef_density(Matrix(result.cots_mechanism_log[:, channel, :]), m)
                    for channel in 8:10] : Vector{Float64}[]
            reef_fast_coral = aggregate_reef_density(fast_coral, m)
            reef_coral_site = aggregate_reef_density(coral_site, m)
            reef_fast_site = aggregate_reef_density(fast_site, m)
            scale = MIGRATION ? Q_FIXED : fixed_scales[reef]
            for (idx, year) in enumerate(years)
                if EXPERIMENT == "structural_crash" || MIGRATION
                    push!(mechanisms, (t.name, reef, year,
                        reef_stages[1][idx], reef_stages[2][idx], reef_stages[3][idx],
                        reef_flows[5][idx], reef_fast_coral[idx], coral_reef[idx],
                        (values[idx] for values in reef_mechanisms)...,
                        reef_coral_site[idx], reef_fast_site[idx],
                        (values[idx] for values in reef_collapse)...))
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
        println("Completed ", t.name, " in ", round(time() - arm_started; digits=1), " s (run_scenario ",
            round(scenario_timing.time; digits=1), " s, gc ", round(scenario_timing.gctime; digits=1),
            " s, alloc ", round(scenario_timing.bytes / 1e9; digits=2), " GB, setup ",
            round(setup_seconds; digits=1), " s [", setup_breakdown, "], post ", round(time() - post_started; digits=1), " s, live heap ",
            round(Base.gc_live_bytes() / 1e9; digits=2), " GB)")
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
    (EXPERIMENT == "structural_crash" || MIGRATION) &&
        CSV.write(joinpath(OUT, "mechanisms.csv"), mechanisms)
    SECOND_WAVE_AUDIT && CSV.write(joinpath(OUT, "site_audit.csv"), site_audit)
    CULLING && CSV.write(joinpath(OUT, "culls.csv"), culls)
end
delete!(ENV, "COTS_CULL_REMOVALS_CSV")
delete!(ENV, "COTS_INITIAL_STATE_CSV")
delete!(ENV, "ADRIA_COTS_ENABLED")
delete!(ENV, "COTS_DIAGNOSTIC_GRAZING_OFF_START_YEAR")
delete!(ENV, "COTS_DIAGNOSTIC_GRAZING_OFF_END_YEAR")
delete!(ENV, "COTS_LAGGED_ADULT_HAZARD")
delete!(ENV, "COTS_HAZARD_FOOD_COUPLED")
delete!(ENV, "COTS_HAZARD_MAX")
delete!(ENV, "ADRIA_RNG_SEED")
delete!(ENV, "ADRIA_DEBUG_BLEACHING_SCALAR")
delete!(ENV, "ADRIA_BLEACHING_DHW_FLOOR")
delete!(ENV, "ADRIA_BLEACHING_LOW_DHW_SCALE")
delete!(ENV, "ADRIA_BLEACHING_LOW_DHW_THRESHOLD")
domain.loc_data.k .= baseline_k
domain.cyclone_mortality_scens = zero_cyclones
domain.loc_data.CB_CALIB_GROUPS .= base_cb_groups
if REBASE_RETENTION
    for i in eachindex(base_reef_connectivity)
        domain.cots_forcing.reef_connectivity[i] = base_reef_connectivity[i]
    end
end
for name in ("COTS_PREFERRED_PREY_STARVATION", "COTS_STARVATION_PREFERRED_THRESHOLD",
        "COTS_STARVATION_MEMORY_YEARS", "COTS_ADULT_SENESCENCE", "COTS_SENESCENCE_AGE",
        "COTS_SENESCENCE_MORTALITY", "COTS_PER_CAPITA_CONSUMPTION",
        "COTS_CONSUMPTION_M2_PER_ADULT_YEAR", "COTS_PREY_SWITCHING",
        "COTS_SWITCH_MIDPOINT", "COTS_SWITCH_STEEPNESS", "COTS_PREFERRED_CONDITION",
        "COTS_CONDITION_SATURATION", "COTS_CONDITION_MEMORY_YEARS",
        "COTS_DIET_LINKED_STARVATION", "COTS_NONPREFERRED_FEEDING_EFFICIENCY",
        "COTS_JUVENILE_PREDATION", "COTS_JUVENILE_PREDATION_CAPACITY",
        "COTS_PREFERRED_MATURATION", "COTS_ALLEE_STRETCHED", "COTS_ALLEE_EXPONENT",
        "COTS_OUTBREAK_COLLAPSE", "COTS_COLLAPSE_THRESHOLD_HA", "COTS_COLLAPSE_BACKGROUND_HA",
        "COTS_COLLAPSE_DURATION_YEARS", "COTS_COLLAPSE_JUVENILES",
        "COTS_STARVATION_POST_GRAZING", "COTS_STARVATION_GRACE_YEARS", "COTS_FORCED_ADULTS_CSV")
    delete!(ENV, name)
end

metadata = Dict(
    "status" => "bounded_1991_$(EXPERIMENT)_screen_not_promoted",
    "experiment" => EXPERIMENT,
    "arm_filter" => ARM_FILTER, "arms_run" => [t.name for t in run_treatments],
    "grazing_off_years" => EXPERIMENT == "gap_grazing" ? [2004, 2012] : Int[],
    "year_window" => [1991, 2024], "score_start_year" => 1992,
    "seed" => SEED, "allee_threshold_adults_ha" => 3.0,
    "cots_initial_profile" => EXPERIMENT == "initialization" ?
        "adult-only for IDW treatments; inherited 6:3:1 for control" :
        EXPERIMENT == "cohort_history" ?
        "IDW adult CPUE with 0.5 or 1.0 times inherited 6:3:1 immature:adult ratio; no-COTS coral counterfactual" :
        EXPERIMENT in ("crash_factorial", "lagged_hazard", "cohort_audit", "structural_crash") || MIGRATION ?
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
    "cots_rebaseline_screen" => REBASE ? Dict(
        "domain" => DOMAIN_PATH,
        "domain_dhw_sha256" => file_sha256(joinpath(DOMAIN_PATH, "DHWs", "dhw_RCPhistorical.nc")),
        "domain_cyclone_sha256" => file_sha256(joinpath(DOMAIN_PATH, "cyclones", "cyclone_mortality.csv")),
        "domain_initial_cover_sha256" => file_sha256(joinpath(DOMAIN_PATH, "initial_cover.nc")),
        "mellin_background_csv" => MELLIN_BACKGROUND_CSV,
        "mellin_background_sha256" => file_sha256(MELLIN_BACKGROUND_CSV),
        "photo_transect_csv" => PHOTO_TRANSECT_CSV, "photo_transect_sha256" => file_sha256(PHOTO_TRANSECT_CSV),
        "base_model" => "run-15 acrpoc_g0: per-capita c=$(PER_CAPITA_M2_PER_ADULT_YEAR), H1, H2, preferred condition, no switching, v0.3, Ext x$(EXTERNAL_MULTIPLIER), P (Q per arm), gradual Allee, RME CB groups (group 5), growth untouched, observed 1993-95 composition",
        "observed_fast_share_by_reef" => composition.reef_share,
        "working_coral" => "kappa $(REBASE_KAPPA), Mellin table + $(REBASE_CHRONIC)/yr chronic, bleaching beta per arm, s = f = 1",
        "mortality" => "1 - (1 - m_cyclone)^kappa * (1 - m_mellin[site, year]) * (1 - mu); bleaching loss scaled by ADRIA_DEBUG_BLEACHING_SCALAR",
        "initial_share" => REBASE_SHARE ? Dict(
            "share_levels_json" => SHARE_LEVELS_JSON, "share_levels_sha256" => file_sha256(SHARE_LEVELS_JSON),
            "levels" => share_levels.raw, "picks" => share_levels.picks,
            "coral_v2" => "ADRIA_BLEACHING_LOW_DHW_SCALE $(CV2_LOW_DHW_SCALE) below 4 DHW; Mellin + $(CV2_MU)/yr; kappa 5; beta 1.5 (re-fit v2 selection)",
            "growth_jl_sha256" => file_sha256(joinpath(ROOT, "ADRIA", "src", "ecosystem", "corals", "growth.jl")),
        ) : "not tested",
        "coral_mix" => CORALMIX ? Dict(
            "chronic_mu" => collect(CORALMIX_MU), "low_dhw_scale" => collect(CORALMIX_LOW),
            "mortality" => "1 - (1 - m_cyclone)^5 * (1 - m_mellin) * (1 - mu); bleaching scalar 1.5; low-DHW threshold 4",
            "shares" => TIMING_SHARE_SETS["pnd"], "seed_pair" => collect(CANDIDATE_SEEDS[1]),
        ) : "not tested",
        "culling" => CULLING ? Dict(
            "switch" => "COTS_CULL_REMOVALS_CSV (unset = established path)",
            "tables" => Dict(c => Dict("path" => cull_table("_$(c)_"), "sha256" => file_sha256(cull_table("_$(c)_")))
                for c in CULL_ARMS if c != "nocull"),
            "builder_metadata_sha256" => file_sha256(joinpath(CULL_DIR, "metadata.json")),
            "timing" => "pre (Sep-Dec dives) before the COTS step; post (Jan-Aug dives) after it, before dispersal",
            "stages" => "juveniles = cohort 1 (<15 cm) from N[2]; adults = cohorts 2-4 from N[3]; capped at present",
            "cots_jl_sha256" => file_sha256(joinpath(ROOT, "ADRIA", "src", "ecosystem", "cots.jl")),
            "scenario_jl_sha256" => file_sha256(joinpath(ROOT, "ADRIA", "src", "scenario.jl")),
        ) : "not tested",
        "candidates" => CANDIDATES ? Dict(
            "arms" => [Dict("suffix" => s, "observation_q" => q) for (s, q) in CANDIDATE_ARMS],
            "seed_pairs_connectivity_adria_rng" => [collect(s) for s in CANDIDATE_SEEDS],
            "archived_pair" => [SEED, archived_rng_seed], "julia_random_seed" => SEED,
            "share_sets" => TIMING_SHARE_SETS, "other_reefs" => "q50",
            "scoring" => "runner scores at q = 0.015; obs030 arms are read via rescore_detectability.jl --q=0.030; primary score v2",
        ) : "not tested",
        "timing" => TIMING ? Dict(
            "predation_capacity_q_cots_ha_yr" => collect(TIMING_Q), "observation_q" => collect(TIMING_OBS_Q),
            "share_sets" => TIMING_SHARE_SETS, "other_reefs" => "q50",
            "seeds" => "eb surface CPUE / q (linear, arm's q), 6:3:1 preload",
            "scoring" => "runner scores at q = 0.015; obs030 arms are read via rescore_detectability.jl --q=0.030; primary score v2 (score_v2.jl)",
        ) : "not tested",
        "lizard_coral" => REBASE_LIZCORAL ? Dict(
            "lizard_scale" => collect(LIZCORAL_SCALES), "observation_q" => collect(LIZCORAL_Q), "chronic_mu" => collect(LIZCORAL_MU),
            "lizard_sites" => length(mapping["Lizard Island Reef"].indices),
            "evidence" => "Pratchett 2005 reef-wide hard coral 32.2% (Oct 1996) vs LTMP 9 m manta 12.7% / photo 10.7%",
            "scoring" => "runner scores at q = 0.015; q = 0.030 arms are read via rescore_detectability.jl --q=0.030",
            "seeds" => "eb surface CPUE / q (linear), 6:3:1 preload",
        ) : "not tested",
        "bleaching_floor" => REBASE_FLOOR ? Dict(
            "switch" => "ADRIA_BLEACHING_DHW_FLOOR (unset = HEAT_LB = 1.0, established path); floor levels $(collect(FLOOR_LEVELS))",
            "meaning" => "minimum bleaching susceptibility (DHW-weeks): no mortality at or below; tolerance distribution truncated there",
            "growth_jl_sha256" => file_sha256(joinpath(ROOT, "ADRIA", "src", "ecosystem", "corals", "growth.jl")),
            "share_multipliers" => collect(SHARE_MULTIPLIERS), "share_cap" => SHARE_CAP,
            "share_by_reef_m1" => composition.reef_share,
            "test" => "sandbox/calibration/test_bleaching_floor.jl",
        ) : "not tested",
        "local_retention" => REBASE_RETENTION ? Dict(
            "transform" => "connectivity_retention.jl retain_within_groups: rows sources, columns sinks; source s in group G keeps rho of its surviving larvae (row sum) within G, Owen within-G pattern kept, exports scaled; row sums preserved",
            "transform_sha256" => file_sha256(joinpath(@__DIR__, "connectivity_retention.jl")),
            "levels" => Dict(k => isnothing(v) ? "owen (unchanged)" : "$(v[1]) rho=$(v[2])" for (k, v) in RETENTION_LEVELS),
            "lizard_complex_gbrmpa_ids" => collect(LIZARD_COMPLEX_IDS),
            "lizard_complex_reef_indices" => retention.lizard,
            "owen_lizard_retained_share_by_year" => retention.owen_lizard_share,
            "max_abs_row_sum_change" => retention.row_delta,
            "connectivity_years" => domain.cots_forcing.years,
        ) : "not tested",
        "initial_state" => REBASE_INITIAL || REBASE_RETENTION ? Dict(
            "eb_surface_dir" => EB_SURFACE_DIR,
            "eb_surface_sha256" => file_sha256(joinpath(EB_SURFACE_DIR, "surface.csv")),
            "eb_metadata_sha256" => file_sha256(joinpath(EB_SURFACE_DIR, "metadata.json")),
            "inversion" => "lin: A = C / q; h10: A = (C/q + sqrt((C/q)^2 + 4 A_half C/q)) / 2, inverse of C = q A^2 / (A + A_half), A_half = $(INITIAL_A_HALF)",
            "preload" => "6:3:1 juvenile:subadult:adult on the inverted adults",
        ) : "idw surface, linear inversion",
        "arms" => CORALMIX ? "wc_b150_q20_eb_lin_obs015_s1_m{030,060,090,120}_l{100,050,025}_cv2_pnd + wc_b150_q20_eb_lin_obs015_s1_ref" :
            CULLING ? "wc_b150_q20_eb_lin_obs015_s1_{nocull,cull000,cull210,cull420}_cv2_pnd" :
            CANDIDATES ? "wc_b150_q20_eb_lin_obs{015,030}_s{1..5}_{cv2_pnd,cv2_pick,cv2_pndl50,ref} (see candidates)" :
            TIMING ? "wc_b150_q{10,15,20}_eb_lin_obs{015,030}_cv2_{pick,pnd}" :
            REBASE_SHARE ? "wc_b150_q20_eb_lin_{ref, cv2_acrpoc, cv2_min..cv2_max, cv2_pick}" :
            REBASE_LIZCORAL ?
            "wc_b150_q20_eb_lin_obs{015,030}_mu{030,080}_c{100,150,200}" :
            REBASE_FLOOR ?
            "wc_b{100,150}_q20_eb_lin_f{100,400}_m{100,200,300}" :
            REBASE_RETENTION ?
            "wc_b150_q20_eb_lin_{owen,liz25,liz50,liz90,all50}" :
            REBASE_INITIAL ?
            "wc_b150_q{20,30}_{idw,eb}_{lin,h10}" :
            REBASE_COMPOSITION ?
            "wc_b{100,150}_q20_{acrpoc,halfother,other}: 1991 fast share (Acropora + Pocilloporidae + w * Other) / hard coral, w 0 / 0.5 / 1" :
            "ref_q30 (kappa 1, no background, beta unset); k500_q30 (kappa 5 only); mel3_q30 (Mellin + 3% only); wc_b{100,150}_q{$(join(REBASE_Q, ","))}",
        "fast_share_by_level" => Dict(level => Dict("reefs" => v.reef_share, "default" => v.default_share,
            "unchanged_sites" => v.unchanged_sites) for (level, v) in composition.levels),
        "bleaching_scalar" => REBASE_INITIAL || REBASE_RETENTION || REBASE_LIZCORAL ? [1.5] : collect(REBASE_BETA),
        "predation_capacity_q_cots_ha_yr" => CANDIDATES || CULLING || CORALMIX ? [20] : TIMING ? collect(TIMING_Q) : REBASE_COMPOSITION || REBASE_RETENTION || REBASE_FLOOR || REBASE_LIZCORAL ? [20] : REBASE_INITIAL ? [20, 30] : collect(REBASE_Q),
        "selection" => REBASE_FLOOR ?
            "m: min summed loss among arms with G2' Lizard/MacGillivray trough no worse than ref and manta total RMSE <= 1.1 x ref; floor 4 by evidence if manta coral J not worse; beta by coral J" :
            REBASE_RETENTION ?
            "none: diagnostic; gate G2' (evaluate_gate_g2prime.py) on the re-score; discrimination vs all50 control" :
            REBASE_INITIAL ?
            "none: diagnostic of first-wave ranking and heights; h10 arms read primarily under the matching A_half = 10 observation re-score" :
            REBASE_COMPOSITION ?
            "adoption only if Lizard and MacGillivray first peaks rise, mean trough min <= 4 adults/ha, summed loss worsens <= 0.5, fast RMSE worsens <= 0.01; lead decides" :
            "Q per beta: min summed smoothed_3y loss over 4 reefs, first peak on every reef, ties to smaller Q; beta by lower frozen coral J at its selected Q",
        "reference" => CORALMIX ?
            "m120_l025 must replay run 37 wc_b150_q20_eb_lin_obs015_s1_cv2_pnd; ref must replay run 37 wc_b150_q20_eb_lin_obs015_s1_ref" :
            CULLING ?
            "wc_b150_q20_eb_lin_obs015_s1_nocull_cv2_pnd must replay run 32 wc_b150_q20_eb_lin_obs015_cv2_pnd" :
            CANDIDATES ?
            "s1 arms replay: obs015_s1_cv2_pick -> run 25 cv2_pick; obs015_s1_cv2_pnd / obs030_s1_cv2_pnd -> run 32 obs015 / obs030 cv2_pnd; obs015_s1_ref -> run 25 ref" :
            TIMING ?
            "wc_b150_q20_eb_lin_obs015_cv2_pick must replay wc_b150_q20_eb_lin_cv2_pick of 20261014T_25_initial_share_1991_v2" :
            REBASE_SHARE ?
            "wc_b150_q20_eb_lin_ref must replay wc_b150_q20_eb_lin of 20261013T_17_initial_state_rebaseline_1991_v1" :
            REBASE_LIZCORAL ?
            "wc_b150_q20_eb_lin_obs015_mu030_c100 must replay wc_b150_q20_eb_lin of 20261013T_17_initial_state_rebaseline_1991_v1" :
            REBASE_FLOOR ?
            "wc_b150_q20_eb_lin_f100_m100 must replay wc_b150_q20_eb_lin of 20261013T_17_initial_state_rebaseline_1991_v1" :
            REBASE_RETENTION ?
            "wc_b150_q20_eb_lin_owen must replay wc_b150_q20_eb_lin of 20261013T_17_initial_state_rebaseline_1991_v1" :
            REBASE_INITIAL ?
            "wc_b150_q20_idw_lin must replay wc_b150_q20 of 20261013T_11_cots_rebaseline_1991_v1" :
            REBASE_COMPOSITION ?
            "wc_b{100,150}_q20_acrpoc must replay wc_b{100,150}_q20 of 20261013T_11_cots_rebaseline_1991_v1" :
            "ref_q30 must replay acrpoc_g0 of 20261012T_15_food_grace_1991_v1",
    ) : "not tested",
    "coral_refit_v2_screen" => REFIT2 ? Dict(
        "forced_tables" => Dict(k => Dict("path" => v, "sha256" => file_sha256(v)) for (k, v) in refit_forcing),
        "survey_table_metadata_sha256" => file_sha256(joinpath(dirname(FORCED_ADULTS_SURVEY_CSV), "metadata.json")),
        "low_dhw" => "ADRIA_BLEACHING_LOW_DHW_SCALE s in $(collect(REFIT2_SCALES)) below 4 DHW (unset for s = 1)",
        "chronic_mu" => collect(REFIT2_MU), "kappa" => REFIT_KAPPA, "bleaching_scalar" => 1.5, "forcing_q" => 0.015,
        "growth_jl_sha256" => file_sha256(joinpath(ROOT, "ADRIA", "src", "ecosystem", "corals", "growth.jl")),
        "mellin_background_sha256" => file_sha256(MELLIN_BACKGROUND_CSV),
        "reference" => "refit2_s100_mu030_rep must replay bleach_mel3_k500_b150 of 20261013T_09_bleaching_sensitivity_1991_v1",
    ) : "not tested",
    "coral_refit_screen" => REFIT && !REFIT2 ? Dict(
        "domain" => DOMAIN_PATH,
        "forced_adults_csv" => FORCED_ADULTS_CSV, "forced_adults_sha256" => file_sha256(FORCED_ADULTS_CSV),
        "forced_tables" => Dict(string(q) => Dict("path" => path, "sha256" => file_sha256(path)) for (q, path) in refit_forcing),
        "forcing_q_note" => "q = 0.030 is the Poisson fit to the Pratchett (2005) Lizard pairs; adults = 0.015-table x 0.015 / q",
        "mellin_background_csv" => MELLIN_BACKGROUND_CSV, "mellin_background_sha256" => file_sha256(MELLIN_BACKGROUND_CSV),
        "growth_jl_sha256" => file_sha256(joinpath(ROOT, "ADRIA", "src", "ecosystem", "corals", "growth.jl")),
        "floors" => collect(REFIT_FLOORS), "chronic_mu" => collect(REFIT_MU), "forcing_q" => collect(REFIT_Q),
        "kappa" => REFIT_KAPPA, "bleaching_scalar" => 1.5,
        "mortality" => "1 - (1 - m_cyclone)^kappa * (1 - m_mellin) * (1 - mu); bleaching floor via ADRIA_BLEACHING_DHW_FLOOR (unset for floor 1)",
        "reference" => "refit_f100_mu030_q015 must replay bleach_mel3_k500_b150 of 20261013T_09_bleaching_sensitivity_1991_v1",
    ) : "not tested",
    "bleaching_sensitivity_screen" => BLEACH && !REFIT ? Dict(
        "domain" => DOMAIN_PATH,
        "domain_dhw_sha256" => file_sha256(joinpath(DOMAIN_PATH, "DHWs", "dhw_RCPhistorical.nc")),
        "domain_cyclone_sha256" => file_sha256(joinpath(DOMAIN_PATH, "cyclones", "cyclone_mortality.csv")),
        "forced_adults_csv" => FORCED_ADULTS_CSV, "forced_adults_sha256" => file_sha256(FORCED_ADULTS_CSV),
        "mellin_background_csv" => MELLIN_BACKGROUND_CSV,
        "mellin_background_sha256" => file_sha256(MELLIN_BACKGROUND_CSV),
        "growth" => "CoralBlox calibrated group factors untouched (s = f = 1), RME groups",
        "background_levels" => Dict(BLEACH_BACKGROUND), "cyclone_kappa" => collect(BLEACH_KAPPA),
        "bleaching_scalar" => collect(BLEACH_BETA),
        "mortality" => "1 - (1 - m_cyclone)^kappa * (1 - m_mellin) * (1 - mu); bleaching loss scaled by ADRIA_DEBUG_BLEACHING_SCALAR",
        "reference" => "bleach_mel_k250_b100 must replay mellin_s100_f100_k250 of 20261013T_03_mellin_coral_1991_v1",
    ) : "not tested",
    "mellin_coral_screen" => MELLIN ? Dict(
        "domain" => DOMAIN_PATH,
        "domain_cyclone_sha256" => file_sha256(joinpath(DOMAIN_PATH, "cyclones", "cyclone_mortality.csv")),
        "forced_adults_csv" => FORCED_ADULTS_CSV, "forced_adults_sha256" => file_sha256(FORCED_ADULTS_CSV),
        "mellin_background_csv" => MELLIN_BACKGROUND_CSV,
        "mellin_background_sha256" => file_sha256(MELLIN_BACKGROUND_CSV),
        "mellin_metadata_sha256" => file_sha256(joinpath(dirname(MELLIN_BACKGROUND_CSV), "metadata.toml")),
        "slow_multipliers" => collect(JOINT_SLOW), "fast_multipliers" => collect(JOINT_FAST),
        "cyclone_kappa" => collect(JOINT_KAPPA),
        "mortality" => "1 - (1 - m_cyclone)^kappa * (1 - m_mellin[site, year]), all groups and sizes",
        "paired_controls" => "mu = 0 arms of 20261013T_02_joint_coral_1991_v1",
    ) : "not tested",
    "joint_coral_screen" => JOINT ? Dict(
        "domain" => DOMAIN_PATH,
        "domain_initial_cover_sha256" => file_sha256(joinpath(DOMAIN_PATH, "initial_cover.nc")),
        "domain_cyclone_sha256" => file_sha256(joinpath(DOMAIN_PATH, "cyclones", "cyclone_mortality.csv")),
        "forced_adults_csv" => FORCED_ADULTS_CSV,
        "forced_adults_sha256" => file_sha256(FORCED_ADULTS_CSV),
        "composition" => "observed 1993-95 photo-transect Acropora+Pocilloporidae share, inherited 2023 totals",
        "slow_multipliers" => collect(JOINT_SLOW), "fast_multipliers" => collect(JOINT_FAST),
        "cyclone_kappa" => collect(JOINT_KAPPA), "background_mortality" => collect(JOINT_BACKGROUND),
        "mortality" => "m' = 1 - (1 - m_cyclone)^kappa * (1 - mu) per site-year, all groups and sizes, via the cyclone mortality step before COTS grazing",
        "rme_groups_present" => rme_cb.rme_ids,
        "reference" => "forced_s100_f100_k100_b000 must replay forced_s100_f100 of 20261013T_01_coral_calibration_1991_v2",
    ) : "not tested",
    "coral_calibration_screen" => CORAL ? Dict(
        "domain" => DOMAIN_PATH,
        "domain_initial_cover_sha256" => file_sha256(joinpath(DOMAIN_PATH, "initial_cover.nc")),
        "forced_adults_csv" => FORCED_ADULTS_CSV,
        "forced_adults_sha256" => file_sha256(FORCED_ADULTS_CSV),
        "forced_metadata_sha256" => file_sha256(joinpath(dirname(FORCED_ADULTS_CSV), "metadata.json")),
        "forcing" => "forced arms: ADRIA_COTS_ENABLED=false; COTS_FORCED_ADULTS_CSV per-capita preferred-first grazing c=$(PER_CAPITA_M2_PER_ADULT_YEAR) m2/adult/yr divided by site k",
        "composition" => "observed 1993-95 photo-transect Acropora+Pocilloporidae share, inherited 2023 totals",
        "observed_fast_share_by_reef" => composition.reef_share,
        "slow_multipliers" => collect(SLOW_MULTIPLIERS), "fast_multipliers" => collect(FAST_MULTIPLIERS),
        "multiplier_scope" => "linear_extension_scale of small_massives+large_massives (slow) and groups 1-3 (fast) for every RME CB group present",
        "rme_groups_present" => rme_cb.rme_ids,
        "reference" => "dynamic_ref (COTS dynamics, leading candidate) must replay acrpoc_g0 of 20261012T_15_food_grace_1991_v1",
    ) : "not tested",
    "food_grace_screen" => GRACE ? Dict(
        "domain" => DOMAIN_PATH,
        "domain_initial_cover_sha256" => file_sha256(joinpath(DOMAIN_PATH, "initial_cover.nc")),
        "base_model" => "run-14 obs_m100_q30: working model, v0.3, Ext x$(EXTERNAL_MULTIPLIER), P Q=$(COLLAPSE_Q), gradual Allee, RME CB groups (group 5), fast extension x1.0",
        "photo_transect_csv" => PHOTO_TRANSECT_CSV, "photo_transect_sha256" => file_sha256(PHOTO_TRANSECT_CSV),
        "composition" => "1991 fast share = (Acropora + Pocilloporidae + w * Other) / all hard coral, 1993-95 mean, 9 m; w: acrpoc 0, halfother 0.5, other 1; totals and within-block structure kept",
        "fast_share_by_level" => Dict(level => Dict("reefs" => v.reef_share, "default" => v.default_share,
            "unchanged_sites" => v.unchanged_sites) for (level, v) in composition.levels),
        "grace" => "COTS_STARVATION_GRACE_YEARS: H1 starvation only after more than G consecutive years with remembered preferred food <= threshold; G = 0 or 1",
        "reference" => "acrpoc_g0 must replay obs_m100_q30 of 20261012T_14_initial_composition_1991_v1",
    ) : "not tested",
    "initial_composition_screen" => COMPOSITION ? Dict(
        "domain" => DOMAIN_PATH,
        "domain_gpkg_sha256" => file_sha256(joinpath(DOMAIN_PATH, "spatial", "lizard_cluster.gpkg")),
        "domain_initial_cover_sha256" => file_sha256(joinpath(DOMAIN_PATH, "initial_cover.nc")),
        "base_model" => "per-capita c=$(PER_CAPITA_M2_PER_ADULT_YEAR), H1, H2, preferred condition, no switching, bleaching on, v0.3, Ext x$(EXTERNAL_MULTIPLIER), gradual Allee, per-reef RME CB groups (group 5 for calibration reefs)",
        "photo_transect_csv" => PHOTO_TRANSECT_CSV, "photo_transect_sha256" => file_sha256(PHOTO_TRANSECT_CSV),
        "photo_selection" => "purpose=COMPOSITION, variable=HARD CORAL, depth=9, years $(first(PHOTO_YEARS))-$(last(PHOTO_YEARS)); fast = Acropora + Pocilloporidae over sum of all hard-coral categories; median",
        "photo_reef_mapping" => PHOTO_REEFS,
        "observed_fast_share_by_reef" => composition.reef_share,
        "default_fast_share_other_reefs" => composition.default_share,
        "sites_left_unchanged_zero_block" => composition.unchanged_sites,
        "arms" => "inh_m110_q30 (inherited 2023 composition; replays m110_ref of run 13), obs_m110_q{20,30,45,70,100}, obs_m{100,120}_q30",
    ) : "not tested",
    "growth_starvation_screen" => GROWTH ? Dict(
        "domain" => DOMAIN_PATH,
        "domain_datapackage_sha256" => file_sha256(joinpath(DOMAIN_PATH, "datapackage.json")),
        "domain_cyclone_sha256" => file_sha256(joinpath(DOMAIN_PATH, "cyclones", "cyclone_mortality.csv")),
        "domain_gpkg_sha256" => file_sha256(joinpath(DOMAIN_PATH, "spatial", "lizard_cluster.gpkg")),
        "base_model" => "per-capita c=$(PER_CAPITA_M2_PER_ADULT_YEAR), H1, H2, preferred condition, no switching, bleaching on, v0.3, Ext x$(EXTERNAL_MULTIPLIER), P Q=$(COLLAPSE_Q) (tuned), gradual Allee, per-reef RME CB groups (project lead: group 5)",
        "rme_layer" => RME_GBR_LAYER, "rme_layer_sha256" => file_sha256(RME_GBR_LAYER),
        "rme_groups_present" => rme_cb.rme_ids,
        "rme_group_of_calibration_reefs" => Dict(t.model => rme_cb.reef_groups[String(domain.loc_data.GBRMPA_ID[first(mapping[t.model].indices)])] for t in TARGETS),
        "fast_linear_extension_multipliers" => collect(GROWTH_MULTIPLIERS),
        "fast_groups" => "tabular_Acropora, corymbose_Acropora, corymbose_non_Acropora (COTS prey map 1-3)",
        "starvation" => "ref: memory $(STRUCTURAL_MEMORY_YEARS) yr on pre-grazing fast cover; M0: memory 0; PG: memory $(STRUCTURAL_MEMORY_YEARS) yr on post-grazing fast cover max(F - N3*c/(1e4 k), 0)",
        "reference" => "m100_ref must replay grme_none of 20261012T_12_outbreak_collapse_1991_v1",
    ) : "not tested",
    "outbreak_collapse_screen" => COLLAPSE ? Dict(
        "domain" => DOMAIN_PATH,
        "domain_datapackage_sha256" => file_sha256(joinpath(DOMAIN_PATH, "datapackage.json")),
        "domain_cyclone_sha256" => file_sha256(joinpath(DOMAIN_PATH, "cyclones", "cyclone_mortality.csv")),
        "domain_gpkg_sha256" => file_sha256(joinpath(DOMAIN_PATH, "spatial", "lizard_cluster.gpkg")),
        "base_model" => "per-capita c=$(PER_CAPITA_M2_PER_ADULT_YEAR), H1, H2, preferred condition, no switching, bleaching on, v0.3, Ext x$(EXTERNAL_MULTIPLIER), P Q=$(COLLAPSE_Q) (tuned), gradual Allee",
        "collapse" => "after D consecutive years with post-transition adults >= $(COLLAPSE_THRESHOLD_HA) COTS/ha, adults (D2J: N2+adults) rescaled to $(COLLAPSE_BACKGROUND_HA) COTS/ha keeping proportions; D2: D=2, D2J: D=2 with N2, D4: D=4",
        "rme_layer" => RME_GBR_LAYER, "rme_layer_sha256" => file_sha256(RME_GBR_LAYER),
        "rme_groups_present" => rme_cb.rme_ids,
        "rme_group_site_counts" => Dict(string(g) => count(==(g), rme_cb.site_groups) for g in rme_cb.rme_ids),
        "rme_group_of_calibration_reefs" => Dict(t.model => rme_cb.reef_groups[String(domain.loc_data.GBRMPA_ID[first(mapping[t.model].indices)])] for t in TARGETS),
        "coral_arms" => "g1: CB group 1 everywhere (Lizard loader default); grme: per-reef RME groups relabelled contiguously with each RME group's default factors; grmefast: grme with fast-group (1-3) linear-extension scale x$(FAST_LINEAR_EXTENSION_MULTIPLIER)",
        "reference" => "g1_none must replay ext5_PA_q30 of 20261012T_11_predation_capacity_1991_v1",
    ) : "not tested",
    "predation_capacity_grid" => PIT ? Dict(
        "domain" => DOMAIN_PATH,
        "domain_datapackage_sha256" => file_sha256(joinpath(DOMAIN_PATH, "datapackage.json")),
        "domain_provenance_sha256" => file_sha256(joinpath(DOMAIN_PATH, "provenance.json")),
        "domain_cyclone_sha256" => file_sha256(joinpath(DOMAIN_PATH, "cyclones", "cyclone_mortality.csv")),
        "domain_gpkg_sha256" => file_sha256(joinpath(DOMAIN_PATH, "spatial", "lizard_cluster.gpkg")),
        "base_model" => "per-capita c=$(PER_CAPITA_M2_PER_ADULT_YEAR), H1, H2, preferred condition, no switching, bleaching on, v0.3 geomorphic k + RRN cyclones, gradual Allee exp(-(3/A)^0.5)",
        "external_multiplier" => "x$(EXTERNAL_MULTIPLIER) on migrated external outbreak production, all arms (working assumption, project lead 2026-10-12)",
        "capacities_cots_ha_yr" => collect(PIT_CAPACITIES),
        "Q_status" => "tuned (calibrated) parameter; not evidence-bounded on the model's effective age-1 juvenile scale",
        "selection" => "min summed fixed_q0015 smoothed_3y loss over 4 reefs, first peak required on every reef; ties to smaller Q",
        "references" => "ext5_A and ext5_PA_q20 must replay ext5_A and ext5_PA of 20261012T_10_initiation_1991_v3",
    ) : "not tested",
    "initiation_factorial" => INITIATION ? Dict(
        "domain" => DOMAIN_PATH,
        "domain_datapackage_sha256" => file_sha256(joinpath(DOMAIN_PATH, "datapackage.json")),
        "domain_provenance_sha256" => file_sha256(joinpath(DOMAIN_PATH, "provenance.json")),
        "domain_cyclone_sha256" => file_sha256(joinpath(DOMAIN_PATH, "cyclones", "cyclone_mortality.csv")),
        "domain_gpkg_sha256" => file_sha256(joinpath(DOMAIN_PATH, "spatial", "lizard_cluster.gpkg")),
        "base_model" => "per-capita c=$(PER_CAPITA_M2_PER_ADULT_YEAR), H1, H2, preferred condition, no switching, bleaching on, v0.3 geomorphic k + RRN cyclones",
        "external_multiplier" => "x$(EXTERNAL_MULTIPLIER) on migrated external outbreak production in *ext5*/*both* arms (all years)",
        "juvenile_preload" => "x$(JUVENILE_PRELOAD_MULTIPLIER) on 6:3:1 juvenile/subadult preload in *juv5*/*both* arms; adults unchanged",
        "P" => "juvenile predation Q=$(JUVENILE_PREDATION_CAPACITY)", "A" => "fertilisation exp(-(3/A)^0.5)",
        "reference" => "none_none must replay bgeo_none of 20261012T_07_trough_mechanisms_1991_v1",
    ) : "not tested",
    "trough_mechanisms_factorial" => TROUGH ? Dict(
        "base_model" => "per-capita c=$(PER_CAPITA_M2_PER_ADULT_YEAR), H1, H2, preferred condition, no switching, bleaching on",
        "b05" => "k=0.5, no cyclones (b05_none replays benthic v2 k05_nocyc)",
        "bgeo" => "geomorphic-infilled k (floor $(BENTHIC_K_FLOOR); zero initial cover where k=0) + RRN cyclones",
        "geomorphic_k_csv" => GEOMORPHIC_K_CSV, "geomorphic_k_sha256" => file_sha256(GEOMORPHIC_K_CSV),
        "geomorphic_k_provenance_sha256" => file_sha256(replace(GEOMORPHIC_K_CSV, ".csv" => ".provenance.json")),
        "geomorphic_k_mean" => sum(geomorphic_k.k) / length(geomorphic_k.k),
        "geomorphic_k_zero_sites" => count(geomorphic_k.zero),
        "cyclone_hours_csv" => CYCLONE_HOURS_CSV, "cyclone_hours_sha256" => file_sha256(CYCLONE_HOURS_CSV),
        "cyclone_mortality" => "1 - exp($(CYCLONE_LOG_SURVIVAL_PER_HOUR) * hours)",
        "P" => "juvenile predation Q*J/(J+Q), Q=$(JUVENILE_PREDATION_CAPACITY) COTS/ha/yr",
        "D" => "maturation fraction sqrt(remembered fast cover); retained juveniles herbivorous (no coral starvation)",
        "A" => "fertilisation exp(-(3/A)^0.5) (CoCoNet form, scale 3 COTS/ha)",
    ) : "not tested",
    "benthic_cyclone_factorial" => BENTHIC_CYCLONE ? Dict(
        "base" => "per-capita c=$(PER_CAPITA_M2_PER_ADULT_YEAR), H1, H2, preferred condition, no switching, bleaching on",
        "benthic_k_csv" => BENTHIC_K_CSV, "benthic_k_sha256" => file_sha256(BENTHIC_K_CSV),
        "benthic_k_provenance_sha256" => file_sha256(replace(BENTHIC_K_CSV, ".csv" => ".provenance.json")),
        "k_total_mean" => sum(benthic.total) / length(benthic.total),
        "k_mapped_mean" => sum(benthic.mapped) / length(benthic.mapped),
        "k_zero_rule" => "benthic k floored at $(BENTHIC_K_FLOOR) (ADRIA cannot run k = 0 sites); initial cover zero where benthic k = 0",
        "k_zero_sites_floored_total" => count(benthic.total_zero),
        "k_zero_sites_floored_mapped" => count(benthic.mapped_zero),
        "cyclone_hours_csv" => CYCLONE_HOURS_CSV, "cyclone_hours_sha256" => file_sha256(CYCLONE_HOURS_CSV),
        "cyclone_hours_provenance_sha256" => file_sha256(replace(CYCLONE_HOURS_CSV, ".csv" => ".provenance.json")),
        "cyclone_mortality" => "1 - exp($(CYCLONE_LOG_SURVIVAL_PER_HOUR) * hours of >=4 m waves), all groups and sizes",
        "derived_inputs_access" => "internal/licensed sources; outputs in git-ignored sandbox/data",
    ) : "not tested",
    "habitat_disturbance_factorial" => HABITAT ? Dict(
        "base" => "per-capita c=$(PER_CAPITA_M2_PER_ADULT_YEAR), H1, H2, preferred condition, no switching",
        "k_rme_source" => "1 - sand_proportion from RME id_list_2024_12_01.csv via GBRMPA_ID",
        "k_rme_sha256" => file_sha256(RME_ID_LIST),
        "k_rme_mean" => sum(rme_k) / length(rme_k),
        "k_rme_min" => minimum(rme_k), "k_rme_max" => maximum(rme_k),
        "bleaching_counterfactual" => "ADRIA_DEBUG_BLEACHING_SCALAR=0 in *no_bleaching* arms; diagnostic only",
        "cyclones_waves" => "zero; no historical cyclone hindcast available locally",
    ) : "not tested",
    "switching_repair_factorial" => SWITCHING_REPAIR ?
        "all arms: per-capita c=$(PER_CAPITA_M2_PER_ADULT_YEAR), H1 memory, H2 senescence, preferred condition saturation=$(CONDITION_SATURATION) memory=$(CONDITION_MEMORY_YEARS) yr; switching midpoint=$(SWITCH_MIDPOINT) steepness=$(SWITCH_STEEPNESS) except condition_ref; diet-linked starvation in *diet_starvation*; non-preferred efficiency=$(NONPREFERRED_EFFICIENCY) in *low_efficiency* (else 1)" :
        "not tested",
    "switching_condition_factorial" => SWITCHING_CONDITION ?
        "base = per_capita_structural (c=$(PER_CAPITA_M2_PER_ADULT_YEAR) m2/adult/yr; H1 threshold=$(STRUCTURAL_THRESHOLD) memory=$(STRUCTURAL_MEMORY_YEARS) yr; H2 age=$(STRUCTURAL_SENESCENCE_AGE) mortality=$(STRUCTURAL_SENESCENCE_MORTALITY)); prey switching midpoint=$(SWITCH_MIDPOINT) habitable fraction steepness=$(SWITCH_STEEPNESS); preferred condition saturation=$(CONDITION_SATURATION) memory=$(CONDITION_MEMORY_YEARS) yr" :
        "not tested",
    "adria_rng_seed" => MIGRATION ? archived_rng_seed :
        "legacy: floor(sum(parameter values)) per scenario (archived control value $(archived_rng_seed))",
    "observation_model" => MIGRATION ? "fixed_q0015: simulated COTS/tow = 0.015 * adult COTS/ha on every reef" :
        "per-reef frozen scales from 20261001T113401_lizard_domain_gate",
    "unit_migration" => MIGRATION ? Dict(
        "q_cots_tow_per_adult_ha" => Q_FIXED,
        "k_model_unit_to_cots_ha" => K_MIGRATION,
        "k_basis" => "Lizard frozen scale / q; exact for Lizard only",
        "archived" => Dict(string(key) => value for (key, value) in archived_density_params),
        "migrated" => Dict("a_F" => archived_density_params[:a_F] / K_MIGRATION,
            "a_S" => archived_density_params[:a_S] / K_MIGRATION,
            "b_ricker" => archived_density_params[:b_ricker] / K_MIGRATION,
            "IMM" => archived_density_params[:IMM] * K_MIGRATION,
            "external_production_pre_pelagic_recruits_ha_year" => boundary * K_MIGRATION),
        "compat_arm" => "Allee 3k and initial state x k; must replay archived control x k",
        "per_capita_m2_per_adult_year" => PER_CAPITA_M2_PER_ADULT_YEAR,
        "structural_arms" => "preferred starvation threshold=$(STRUCTURAL_THRESHOLD) memory=$(STRUCTURAL_MEMORY_YEARS) yr; senescence age=$(STRUCTURAL_SENESCENCE_AGE) mortality=$(STRUCTURAL_SENESCENCE_MORTALITY)",
    ) : "not applicable",
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
        "archived_cohort_trajectories" => EXPERIMENT in ("crash_factorial", "lagged_hazard", "cohort_audit", "structural_crash") || MIGRATION ?
            file_sha256(joinpath(@__DIR__, "runs", "20261007T199111_owen_1991_cohort_history_final", "trajectories.csv")) : "not_applicable",
        "archived_cohort_initial_state" => EXPERIMENT in ("crash_factorial", "lagged_hazard", "cohort_audit", "structural_crash") || MIGRATION ?
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
              SECOND_WAVE_AUDIT ? ["mechanisms.csv", "site_audit.csv"] :
              EXPERIMENT == "structural_crash" || MIGRATION ? ["mechanisms.csv"] : String[],
              basename.(collect(values(initial_paths))))),
)
open(joinpath(OUT, "metadata.toml"), "w") do io
    TOML.print(io, metadata)
end
isempty(failures.treatment) || error("1991 screen has failed treatments; see failures.csv")
