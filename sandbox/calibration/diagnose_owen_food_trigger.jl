"""Reconstruct site-level food stress from the archived joint-fourfold Owen control.

The transition is unchanged. Adult state plus logged maturation identifies the
food-survival factor; successive body-condition states identify the pre-feeding
food signal when unsaturated. This diagnostic informs, but does not fit, a new
opt-in low-food mortality hypothesis.
"""

using Pkg
const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
Pkg.activate(joinpath(ROOT, "sandbox"))
cd(ROOT)

using ADRIA, CSV, DataFrames, Dates, Random, TOML, Statistics
include(joinpath(@__DIR__, "reef_observation_mapping.jl"))
include(joinpath(@__DIR__, "calibration_run.jl"))

const SEED = 20260930
const RUN_ID = get(ENV, "OWEN_FOOD_DIAGNOSTIC_RUN_ID",
    Dates.format(now(), "yyyymmddTHHMMSS") * "_owen_food_trigger_diagnostic")
occursin(r"^[A-Za-z0-9_.-]+$", RUN_ID) || error("Invalid run ID")
const OUT = joinpath(@__DIR__, "runs", RUN_ID)
ispath(OUT) && error("Refusing to overwrite $OUT")
mkpath(OUT)
const V2 = joinpath(ROOT, "sandbox", "data", "Lizard_Historical_v0.2")
const SCREEN = joinpath(@__DIR__, "runs", "20261001T162000_owen_embedded_mortality_screen")
const REPLICATES = joinpath(@__DIR__, "runs", "20261001T163000_owen_embedded_replicates")
const BEST = joinpath(@__DIR__, "runs", "expanded_cotsconn_pilot_seed20260930",
    "best_summary.csv")
const TARGETS = ["Lizard Island Reef", "MacGillivray Reef",
    "North Direction Reef", "Eyrie Reef"]

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
delete!(ENV, "COTS_INTERNAL_SETTLEMENT_BLACKOUT_START_YEAR")
delete!(ENV, "COTS_INTERNAL_SETTLEMENT_BLACKOUT_END_YEAR")

domain = ADRIA.load_domain(ADRIA.LizardDomain, V2, "historical")
mapping = load_reef_area_weights(domain, V2, TARGETS)
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

# Strict baseline replay before deriving any food diagnostics.
archived = CSV.read(joinpath(REPLICATES, "trajectories.csv"), DataFrame)
archived = archived[(archived.seed .== SEED) .&
    (archived.treatment .== "embedded_joint_4x"), :]
cover = dropdims(sum(result.raw; dims=(2, 3)); dims=(2, 3))
for reef in TARGETS
    m = mapping[reef]
    reference = sort(archived[archived.reef_name .== reef, :], :year)
    adults = aggregate_reef_density(Matrix(result.cots_log[:, 3, :]), m)
    coral = aggregate_reef_density(cover, m)
    nrow(reference) == length(years) || error("Incomplete archived reef")
    years == Int.(reference.year) || error("Baseline years changed")
    Int.(result.cots_forcing_year) == Int.(reference.forcing_year) ||
        error("Baseline hydrodynamic years changed")
    all(isapprox.(adults, Float64.(reference.adults_ha); atol=1e-12, rtol=0)) ||
        error("Baseline adult state changed")
    all(isapprox.(coral, Float64.(reference.coral_cover); atol=1e-12, rtol=0)) ||
        error("Baseline coral state changed")
end

site_rows = DataFrame(reef_name=String[], site_id=String[], year=Int[],
    forcing_year=Int[], area_weight=Float64[], prey_cover_postfeeding=Float64[],
    adult_before_ha=Float64[], adult_after_ha=Float64[],
    juvenile_before_ha=Float64[], matured_juveniles_ha=Float64[],
    condition_before=Float64[], condition_after=Float64[],
    food_signal=Float64[], prey_cover_before_feeding_est=Union{Missing,Float64}[],
    food_survival=Union{Missing,Float64}[], adult_food_deaths_ha=Union{Missing,Float64}[],
    juvenile_food_deaths_ha=Union{Missing,Float64}[])
alpha = 1 / Float64(candidate.tau_condition)
threshold = 0.15 * Float64(candidate.C_max)
max_replay_error = 0.0
for reef in TARGETS
    m = mapping[reef]
    for (site, weight) in zip(m.indices, m.weights), t in 2:length(years)
        adult_before = Float64(result.cots_log[t - 1, 3, site])
        adult_after = Float64(result.cots_log[t, 3, site])
        juvenile_before = Float64(result.cots_log[t - 1, 2, site])
        matured = Float64(result.cots_flow_log[t, 5, site])
        bc_before = Float64(result.cots_condition_log[t - 1, site])
        bc_after = Float64(result.cots_condition_log[t, site])
        signal = clamp((bc_after - (1 - alpha) * bc_before) / alpha, 0.0, 1.0)
        pre_cover = signal < 1 - 1e-9 ? signal * Float64(candidate.C_max) * 0.5 : missing
        food = adult_before > 1e-10 ?
            clamp((adult_after - matured) / ((1 - Float64(candidate.m3)) * adult_before),
                0.0, 1.0) : missing
        if !ismissing(pre_cover) && !ismissing(food)
            expected = pre_cover > threshold ? 1.0 :
                (1 - Float64(candidate.p_tilde)) + Float64(candidate.p_tilde) *
                (pre_cover / threshold)^3
            global max_replay_error = max(max_replay_error, abs(food - expected))
        end
        adult_food_deaths = ismissing(food) ? missing :
            adult_before * (1 - Float64(candidate.m3)) * (1 - food)
        juvenile_food_deaths = ismissing(food) ? missing :
            juvenile_before * (1 - Float64(candidate.m2)) * (1 - food)
        post_cover = sum(result.raw[t, g, s, site] for g in 1:5
            for s in axes(result.raw, 3))
        push!(site_rows, (reef, String(domain.loc_ids[site]), years[t],
            Int(result.cots_forcing_year[t]), weight, post_cover,
            adult_before, adult_after, juvenile_before, matured,
            bc_before, bc_after, signal, pre_cover, food,
            adult_food_deaths, juvenile_food_deaths))
    end
end
max_replay_error <= 1e-8 || error("Food-survival reconstruction failed: $max_replay_error")
CSV.write(joinpath(OUT, "site_food_diagnostics.csv"), site_rows)

summary = combine(groupby(site_rows, [:reef_name, :year]),
    :forcing_year => first => :forcing_year,
    [:area_weight, :adult_before_ha] => ((w, a) -> sum(w .* a)) => :adults_before_ha,
    [:area_weight, :adult_after_ha] => ((w, a) -> sum(w .* a)) => :adults_after_ha,
    [:area_weight, :prey_cover_postfeeding] => ((w, c) -> sum(w .* c)) => :prey_cover_postfeeding,
    [:area_weight, :condition_after] => ((w, c) -> sum(w .* c)) => :condition_area_mean,
    [:area_weight, :adult_before_ha, :condition_after] =>
        ((w, a, c) -> sum(w .* a .* c) / max(sum(w .* a), eps())) =>
        :condition_adult_weighted,
    [:area_weight, :adult_food_deaths_ha] =>
        ((w, d) -> sum(w .* coalesce.(d, 0.0))) => :adult_food_deaths_ha,
    [:area_weight, :juvenile_food_deaths_ha] =>
        ((w, d) -> sum(w .* coalesce.(d, 0.0))) => :juvenile_food_deaths_ha,
)
CSV.write(joinpath(OUT, "reef_food_diagnostics.csv"), summary)

metadata = Dict(
    "status" => "diagnostic_only_no_transition_change",
    "seed" => SEED,
    "state_units" => "COTS ha^-1; cover and condition dimensionless",
    "site_rows" => nrow(site_rows),
    "reef_rows" => nrow(summary),
    "condition_inversion" => "pre-feeding prey cover = food_signal * C_max * 0.5 only when food_signal < 1; saturated signals have missing cover estimate",
    "food_survival_inversion" => "(adult_after - matured_juveniles) / ((1-m3) * adult_before), undefined at zero adult_before",
    "prey_cover_postfeeding_note" => "result.raw includes model coral state after current-year predation; do not use as the survival input",
    "maximum_food_survival_replay_error" => max_replay_error,
    "control_replay" => "adult/coral 1e-12, forcing year exact",
    "julia_version" => string(VERSION),
    "adria_revision" => recorded_revision(ROOT, "OWEN_ADRIA_REVISION"),
    "cotsmod_revision" => recorded_revision(joinpath(ROOT, "..", "COTSMod.jl"),
        "OWEN_COTSMOD_REVISION"),
    "inputs" => Dict(
        "screen_metadata" => file_sha256(joinpath(SCREEN, "metadata.toml")),
        "reference_trajectories" => file_sha256(joinpath(REPLICATES, "trajectories.csv")),
        "candidate" => file_sha256(BEST),
        "script" => file_sha256(@__FILE__),
        "adapter" => file_sha256(joinpath(ROOT, "ADRIA", "src", "ecosystem", "cots.jl")),
        "core" => file_sha256(joinpath(ROOT, "..", "COTSMod.jl", "src", "COTSMod.jl")),
        "owen_forcing" => file_sha256(joinpath(V2, "cots_connectivity", "datasets",
            "owen_global_2018_2023", "provenance.json")),
    ),
    "outputs" => Dict(file => file_sha256(joinpath(OUT, file)) for file in
        ("site_food_diagnostics.csv", "reef_food_diagnostics.csv")),
)
open(joinpath(OUT, "metadata.toml"), "w") do io
    TOML.print(io, metadata; sorted=true)
end
println("Saved ", OUT)
