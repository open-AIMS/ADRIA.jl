"""
Diagnose the per-arm runtime growth seen in long `run_1991_initialization_test.jl`
processes (2026-10-13): run the same forced-COTS coral scenario repeatedly in one
process and record wall time, allocations, GC time, live heap and ADRIA output-directory
size per repetition.

    julia --project=sandbox sandbox/calibration/diagnose_runtime_growth.jl [n_reps] [gc_between]
"""

using ADRIA, DataFrames, CSV, Random, Printf

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const N = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 12
const GC_BETWEEN = length(ARGS) >= 2 && ARGS[2] == "gc"

dom = ADRIA.load_domain(ADRIA.LizardDomain,
    joinpath(ROOT, "sandbox", "data", "Lizard_Historical_v0.3"), "historical")
Random.seed!(20260930)
scenario = ADRIA.sample(dom, 2)[1:1, :]
defaults = ADRIA.param_table(dom)
for name in names(scenario)
    scenario[!, name] .= defaults[1, name]
end
ENV["ADRIA_COTS_ENABLED"] = "false"
ENV["ADRIA_RNG_SEED"] = "90972102"
delete!(ENV, "COTS_FORCED_ADULTS_CSV")  # full 1985-2024 domain; forcing starts 1991
# (was) joinpath(ROOT, "sandbox", "data", "derived",
#     "observed_cots_forcing_1991_2024", "observed_adults_ha.csv")

outputs_size() = isdir("Outputs") ?
    sum(filesize(joinpath(r, f)) for (r, _, fs) in walkdir("Outputs") for f in fs; init=0) : 0

rows = DataFrame(rep=Int[], seconds=Float64[], alloc_gb=Float64[], gc_seconds=Float64[],
    live_heap_gb=Float64[], maxrss_gb=Float64[], outputs_mb=Float64[], env_vars=Int[])
for rep in 1:N
    GC_BETWEEN && GC.gc()
    stats = @timed ADRIA.run_scenario(dom, scenario[1, :])
    push!(rows, (rep, stats.time, stats.bytes / 1e9, stats.gctime, Base.gc_live_bytes() / 1e9,
        Sys.maxrss() / 1e9, outputs_size() / 1e6, length(ENV)))
    @printf("rep %2d  %7.1f s  alloc %6.2f GB  gc %6.1f s  live %5.2f GB  rss %5.2f GB\n",
        rep, stats.time, stats.bytes / 1e9, stats.gctime, Base.gc_live_bytes() / 1e9, Sys.maxrss() / 1e9)
end
out = joinpath(ROOT, "sandbox", "calibration", "runs",
    "diagnose_runtime_growth_$(GC_BETWEEN ? "gc" : "nogc").csv")
CSV.write(out, rows)
println("written ", out)
