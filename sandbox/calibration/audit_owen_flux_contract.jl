"""Audit the frozen V2/Owen production-to-settlement chain without fitting ecology."""

using Pkg
const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
Pkg.activate(joinpath(ROOT, "sandbox"))
cd(ROOT)

using ADRIA, CSV, DataFrames, Dates, Random, SHA, Statistics, TOML
include(joinpath(@__DIR__, "reef_observation_mapping.jl"))

const RUN_ID = get(ENV, "OWEN_FLUX_AUDIT_RUN_ID",
    Dates.format(now(), "yyyymmddTHHMMSS") * "_owen_flux_audit")
occursin(r"^[A-Za-z0-9_.-]+$", RUN_ID) || error("Invalid run ID")
const OUT = joinpath(@__DIR__, "runs", RUN_ID)
ispath(OUT) && error("Refusing to overwrite $OUT")
mkpath(OUT)
const SEED = 20260930
const V2 = joinpath(ROOT, "sandbox", "data", "Lizard_Historical_v0.2")
const FORCING = joinpath(V2, "cots_connectivity", "datasets",
    "owen_global_2018_2023", "annual_forcing")
const RME = joinpath(ROOT, "sandbox", "data", "rme_ml_2025_06_05", "data_files",
    "water_csv", "COTS_LARVAL_REDUCTION_q3baseline.csv")
const CANDIDATE = joinpath(@__DIR__, "runs", "expanded_cotsconn_pilot_seed20260930",
    "best_summary.csv")
const COTSMOD = normpath(joinpath(ROOT, "..", "COTSMod.jl"))
const TARGETS = ["Lizard Island Reef", "MacGillivray Reef",
    "North Direction Reef", "Eyrie Reef"]
sha(path) = bytes2hex(open(sha256, path))

ENV["ADRIA_COTS_CONNECTIVITY_MODE"] = "cots"
ENV["ADRIA_COTS_CONNECTIVITY_DATASET"] = "owen_global_2018_2023"
ENV["COTS_CONNECTIVITY_TEMPORAL_MODE"] = "sample"
ENV["COTS_CALENDAR_START_YEAR"] = "1985"
ENV["COTS_ALLEE_THRESHOLD"] = "3.0"
ENV["COTS_EXTERNAL_PULSE"] = "false"
ENV["COTS_EXTERNAL_SOURCE_DENSITY"] = "0.0"
ENV["COTS_JUVENILE_STORAGE"] = "false"
ENV["COTS_HABITAT_MEDIATION"] = "false"
ENV["COTS_SEPARATED_RECRUITMENT"] = "true"
ENV["COTS_APPLY_LARVAL_SURVIVAL"] = "true"
ENV["COTS_SETTLEMENT_PROBABILITY"] = "1.0"
ENV["COTS_SEED_FIRST_N"] = "10"
ENV["ADRIA_DEBUG_SEED_FIRST_N"] = "10"
best = sort(CSV.read(CANDIDATE, DataFrame), :loss)[1, :]
ENV["COTS_INITIAL_MULTIPLIER"] = string(best.seed_mult)

domain = ADRIA.load_domain(ADRIA.LizardDomain, V2, "historical")
forcing = domain.cots_forcing
mapping = load_reef_area_weights(domain, V2, TARGETS)
legacy_dir = joinpath(V2, "cots_connectivity", "annual_forcing")
legacy = ADRIA.load_cots_connectivity_forcing(legacy_dir, String.(domain.loc_ids))
reference_survival = mean(legacy.source_survival)
fecundity = Float64(best.a_ricker) / reference_survival
ENV["COTS_LARVAL_FECUNDITY"] = string(fecundity)
production(adults) = fecundity * adults * exp(-Float64(best.b_ricker) * adults) *
    adults^2 / (3.0^2 + adults^2)
boundary_ceiling = maximum(production.(range(0.0, 100.0; length=10001)))
ENV["COTS_EXTERNAL_OUTBREAK_PRODUCTION"] = string(boundary_ceiling)

# Check that every applied local source factor equals the RME 2017 column.
site_map = CSV.read(joinpath(FORCING, "site_reef_map.csv"), DataFrame)
survival_table = CSV.read(joinpath(FORCING, "source_larval_survival.csv"), DataFrame)
length(unique(string.(survival_table.reef_id))) == nrow(survival_table) ||
    error("Duplicate local survival reef IDs")
rme = CSV.read(RME, DataFrame)
rme_lookup = Dict(string(row.id) => row for row in eachrow(rme))
reef_ids = [only(unique(string.(site_map.reef_id[site_map.reef_index .== i])))
    for i in 1:length(forcing.reef_habitable_area_ha)]
length(unique(reef_ids)) == length(reef_ids) || error("Duplicate local reef IDs")
source_checks = DataFrame(reef_index=Int[], reef_id=String[], rme_2017=Float64[],
    applied_2018=Float64[], rme_2010_2017_mean=Float64[],
    ratio_2017_to_historical_mean=Float64[])
for (i, id) in enumerate(reef_ids)
    haskey(rme_lookup, id) || error("Missing RME source $id")
    row = rme_lookup[id]
    historical = Float64[row[Symbol(string(y))] for y in 2010:2017]
    applied = forcing.source_survival[i, :]
    supplied = survival_table[survival_table.reef_id .== id, :]
    nrow(supplied) == 1 || error("Missing supplied survival row for $id")
    all(isapprox.(applied, Float64.([supplied[1, Symbol(string(y))]
        for y in forcing.years]); rtol=1e-12, atol=1e-14)) ||
        error("Loaded survival differs from supplied table at $id")
    all(isapprox.(applied, Float64(row[Symbol("2017")]); rtol=1e-10, atol=1e-14)) ||
        error("Applied survival differs from RME 2017 at $id")
    push!(source_checks, (i, id, Float64(row[Symbol("2017")]), applied[1],
        mean(historical), applied[1] / mean(historical)))
end
CSV.write(joinpath(OUT, "source_survival_checks.csv"), source_checks)

Random.seed!(SEED)
scenario = ADRIA.sample(domain, 2)[1:1, :]
defaults = ADRIA.param_table(domain)
for name in names(scenario)
    scenario[!, name] .= defaults[1, name]
    if name != "allee_threshold" && hasproperty(best, Symbol(name))
        scenario[!, name] .= best[Symbol(name)]
    end
end
scenario[!, :allee_threshold] .= 3.0
ENV["COTS_CONNECTIVITY_SEED"] = string(SEED)
Random.seed!(SEED)
result = ADRIA.run_scenario(domain, scenario[1, :])

flows = result.cots_flow_log
years = Int.(result.cots_calendar_year)
years == collect(1985:2024) || error("Unexpected calendar years")
reef_of = forcing.site_to_reef
n_reefs = length(forcing.reef_habitable_area_ha)
sites = [findall(==(reef), reef_of) for reef in 1:n_reefs]
all(!isempty, sites) || error("A forcing reef has no sites")
area = forcing.reef_habitable_area_ha
site_area = Float64.(domain.loc_data.area)
checks = DataFrame(year=Int[], forcing_year=Int[], reef_name=String[],
    reef_id=String[], local_fecundity_ha=Float64[],
    pelagic_expected_ha=Float64[], pelagic_logged_ha=Float64[],
    local_retention_expected_ha=Float64[], local_retention_logged_ha=Float64[],
    internal_expected_ha=Float64[], internal_logged_ha=Float64[],
    external_potential_expected_ha=Float64[], external_potential_logged_ha=Float64[],
    external_pelagic_expected_ha=Float64[], external_pelagic_logged_ha=Float64[],
    external_settled_ha=Float64[], total_settled_ha=Float64[],
    background_immigration_ha=Float64[], recruit_stage_ha=Float64[],
    juvenile_stage_ha=Float64[], adult_stage_ha=Float64[],
    coral_cover=Float64[], source_area_weighted_to_equal_site_ratio=Float64[])
max_residual = Dict(key => 0.0 for key in
    ("pelagic", "local", "internal", "external_potential", "external_pelagic",
     "external_settlement", "settled_balance", "recruit_balance"))

for t in eachindex(years)
    hydro = Int(result.cots_forcing_year[t])
    hydro == 0 && continue
    hydro_index = findfirst(==(hydro), forcing.years)
    isnothing(hydro_index) && error("Unknown selected hydrodynamic year $hydro")
    conn = forcing.reef_connectivity[hydro_index]
    local_production = Float64.(flows[t, 1, :])
    equal_source = [mean(local_production[sites[r]]) for r in 1:n_reefs]
    weighted_source = [sum(local_production[sites[r]] .* site_area[sites[r]]) /
        sum(site_area[sites[r]]) for r in 1:n_reefs]
    source_density = equal_source .* forcing.source_survival[:, hydro_index]
    source_ratio = [equal_source[r] > 0 ? weighted_source[r] / equal_source[r] : NaN
        for r in 1:n_reefs]
    expected_local = zeros(n_reefs)
    expected_internal = zeros(n_reefs)
    for sink in 1:n_reefs, source in 1:n_reefs
        value = source_density[source] * area[source] / area[sink] * conn[source, sink]
        if source == sink
            expected_local[sink] += value
        else
            expected_internal[sink] += value
        end
    end
    boundary_index = findfirst(==(years[t]), forcing.outbreak_years)
    for reef_name in TARGETS
        m = mapping[reef_name]
        mapped_reefs = reef_of[m.indices]
        mapped_ids = join(sort(unique(reef_ids[mapped_reefs])), ";")
        avg(channel) = sum(m.weights .* Float64.(flows[t, channel, m.indices]))
        site_gates = Float64.(flows[t, 7, m.indices])
        # Gate is site-varying; compare the exact site-weighted expected flux.
        local_expected = sum(m.weights .* expected_local[mapped_reefs] .* site_gates)
        internal_expected = sum(m.weights .* expected_internal[mapped_reefs] .* site_gates)
        potential_expected = isnothing(boundary_index) ? 0.0 :
            boundary_ceiling * sum(m.weights .* forcing.external_outbreak_potential[
                mapped_reefs, boundary_index, hydro_index])
        pelagic_external_expected = isnothing(boundary_index) ? 0.0 :
            boundary_ceiling * sum(m.weights .* forcing.external_outbreak_pelagic[
                mapped_reefs, boundary_index, hydro_index])
        external_settled_expected = isnothing(boundary_index) ? 0.0 :
            boundary_ceiling * sum(m.weights .* site_gates .* forcing.external_outbreak_pelagic[
                mapped_reefs, boundary_index, hydro_index])
        pelagic_expected = sum(m.weights .* (expected_local[mapped_reefs] .+
            expected_internal[mapped_reefs]))
        residuals = Dict(
            "pelagic" => abs(avg(9) - pelagic_expected),
            "local" => abs(avg(8) - local_expected),
            "internal" => abs(avg(3) - internal_expected),
            "external_potential" => abs(avg(11) - potential_expected),
            "external_pelagic" => abs(avg(12) - pelagic_external_expected),
            "external_settlement" => abs(avg(4) - external_settled_expected),
            "settled_balance" => abs(avg(10) - avg(8) - avg(3) - avg(4)),
            "recruit_balance" => abs(sum(m.weights .* Float64.(result.cots_log[t, 1, m.indices])) -
                avg(2) - avg(8) - avg(3) - avg(4)))
        for (key, value) in residuals
            max_residual[key] = max(max_residual[key], value)
        end
        coral = dropdims(sum(result.raw[t, :, :, m.indices]; dims=(1, 2)); dims=(1, 2))
        push!(checks, (years[t], hydro, reef_name, mapped_ids, avg(1),
            pelagic_expected, avg(9), local_expected, avg(8), internal_expected,
            avg(3), potential_expected, avg(11), pelagic_external_expected,
            avg(12), avg(4), avg(10), avg(2),
            sum(m.weights .* Float64.(result.cots_log[t, 1, m.indices])),
            sum(m.weights .* Float64.(result.cots_log[t, 2, m.indices])),
            sum(m.weights .* Float64.(result.cots_log[t, 3, m.indices])),
            sum(m.weights .* coral),
            sum(m.weights .* source_ratio[mapped_reefs])))
    end
end
CSV.write(joinpath(OUT, "target_reef_fluxes.csv"), checks)
threshold = 1e-7
all(value < threshold for value in values(max_residual)) ||
    error("Flux contract residual exceeded $threshold: $max_residual")

meta = Dict(
    "run_id" => RUN_ID, "status" => "pass", "seed" => SEED,
    "treatment" => "v2_owen_boundary", "hydrodynamic_mode" => "sample",
    "fecundity_pre_pelagic_ha_year" => fecundity,
    "boundary_ceiling_pre_pelagic_ha_year" => boundary_ceiling,
    "source_survival_policy" => "RME 2017 reused for Owen 2018-2022",
    "source_survival_2017_mean" => mean(source_checks.rme_2017),
    "source_survival_2010_2017_mean" => mean(source_checks.rme_2010_2017_mean),
    "source_area_weighting" => "equal site mean in current COTSMod; area-weighted ratio diagnostic only",
    "adria_revision" => readchomp(`git -C $ROOT rev-parse HEAD`),
    "cotsmod_revision" => readchomp(`git -C $COTSMOD rev-parse HEAD`),
    "code_sha256" => Dict(path => sha(path) for path in
        (@__FILE__, joinpath(ROOT, "ADRIA", "src", "scenario.jl"),
         joinpath(ROOT, "ADRIA", "src", "ecosystem", "cots.jl"),
         joinpath(COTSMOD, "src", "COTSMod.jl"))),
    "max_absolute_residual_ha" => max_residual,
    "source_sha256" => Dict(path => sha(path) for path in
        (RME, CANDIDATE, joinpath(FORCING, "source_larval_survival.csv"),
         joinpath(FORCING, "site_reef_map.csv"),
         joinpath(FORCING, "external_outbreak_boundary.csv"),
         joinpath(FORCING, "provenance.toml"),
         (joinpath(FORCING, "reef_connectivity_$year.csv") for year in forcing.years)...)),
    "output_sha256" => Dict(name => sha(joinpath(OUT, name)) for name in
        ("source_survival_checks.csv", "target_reef_fluxes.csv"))
)
open(joinpath(OUT, "metadata.toml"), "w") do io
    TOML.print(io, meta)
end
println("Audit passed: ", OUT)
println("Maximum absolute residuals COTS ha^-1: ", max_residual)
