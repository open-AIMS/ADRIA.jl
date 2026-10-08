using Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "..")))

using ADRIA
using Random

const DATA_ROOT = normpath(joinpath(@__DIR__, "..", "data"))
const DOMAIN_ROOT = joinpath(DATA_ROOT, "Lizard_Historical_v0.2")
const SMOKE_SEED = 20260930

ENV["ADRIA_COTS_CONNECTIVITY_MODE"] = "cots"
ENV["ADRIA_COTS_CONNECTIVITY_DATASET"] = get(
    ENV, "ADRIA_COTS_CONNECTIVITY_DATASET", "owen_global_2018_2023"
)
ENV["COTS_CONNECTIVITY_TEMPORAL_MODE"] = "cycle"
ENV["COTS_CONNECTIVITY_SEED"] = string(SMOKE_SEED)
ENV["COTS_ALLEE_THRESHOLD"] = "3.0"
ENV["COTS_EXTERNAL_PULSE"] = "false"
ENV["COTS_EXTERNAL_SOURCE_DENSITY"] = "0.0"
ENV["COTS_EXTERNAL_OUTBREAK_PRODUCTION"] = "0.0"

Random.seed!(SMOKE_SEED)
domain = ADRIA.load_domain(ADRIA.LizardDomain, DOMAIN_ROOT, "historical")
scenario = ADRIA.sample(domain, 2)[1:1, :]
defaults = ADRIA.param_table(domain)
for name in names(scenario)
    scenario[!, name] .= defaults[1, name]
end

result = ADRIA.run_scenario(domain, scenario[1, :])
@assert size(result.cots_log, 1) == 40
@assert size(result.cots_log, 3) == 2895
@assert all(isfinite, result.cots_log)
@assert all(result.cots_log .>= 0.0)
@assert all(isfinite, result.raw)
@assert all(result.raw .>= 0.0)
@assert result.cots_forcing_year[1] == 0
@assert all(result.cots_forcing_year[2:end] .∈ Ref([2018, 2019, 2020, 2021, 2022]))

println(
    "V2 scenario smoke test passed: seed=", SMOKE_SEED,
    ", COTS dataset=", ENV["ADRIA_COTS_CONNECTIVITY_DATASET"],
    ", forcing years=", join(unique(result.cots_forcing_year), ','),
    ", max adults=", maximum(result.cots_log[:, 3, :])
)
