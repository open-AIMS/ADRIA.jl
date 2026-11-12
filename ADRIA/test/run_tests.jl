#!/usr/bin/env julia
# Test launcher: runs the ADRIA test suite as supervised subprocesses with hard timeouts.
#
#   julia ADRIA/test/run_tests.jl                 # all shards, concurrently
#   julia ADRIA/test/run_tests.jl fast            # one shard (what CI matrix jobs use)
#   julia ADRIA/test/run_tests.jl fast spec       # a subset, concurrently
#   julia ADRIA/test/run_tests.jl --no-prepare .. # skip env preparation (already done)
#
# Shards (see ADRIA_TEST_SHARD in runtests.jl):
#   fast  tiers 1-2 + precompile budget test
#   spec  tier 3 (spec/sampling/decisions)
#   sim   tiers 4-5 (32-scenario simulation + dependent tests; must stay together)
#
# Each shard runs in its own `julia` process in its own process group. If it exceeds its
# limit the whole group is SIGKILLed (a plain SIGTERM is ignored by Julia while it is busy
# compiling), the shard is reported as TIMEOUT, and the launcher exits non-zero.
# Limits can be overridden with ADRIA_TEST_LIMIT_<SHARD>_S, e.g. ADRIA_TEST_LIMIT_SIM_S=1800.
#
# Output goes to ADRIA/test/logs/<shard>.log (and is echoed if a shard is the only one run).

using Pkg

const TEST_DIR = @__DIR__
const ADRIA_DIR = dirname(TEST_DIR)
const ENV_DIR = joinpath(TEST_DIR, ".testenv")   # shared, prepared once, git-ignored
const LOG_DIR = joinpath(TEST_DIR, "logs")

# Default wall-clock limits in seconds. `fast` and `spec` are 5 min as a "something is
# wrong" tripwire; `sim` runs 32 full simulations and needs longer.
const DEFAULT_LIMITS = Dict("fast" => 300, "spec" => 300, "sim" => 1200)

limit_for(shard) = parse(
    Float64, get(ENV, "ADRIA_TEST_LIMIT_$(uppercase(shard))_S", string(DEFAULT_LIMITS[shard]))
)

"""Create/refresh the shared test environment once, so concurrent shards don't each
instantiate and precompile their own copy."""
function prepare_env()
    mkpath(ENV_DIR)
    cp(joinpath(TEST_DIR, "Project.toml"), joinpath(ENV_DIR, "Project.toml"); force=true)
    @info "Preparing shared test environment" ENV_DIR
    Pkg.activate(ENV_DIR) do
        Pkg.develop(; path=ADRIA_DIR)
        Pkg.instantiate()
        Pkg.precompile()
    end
end

function run_with_timeout(cmd::Cmd, limit_s::Real, logfile::AbstractString)
    open(logfile, "w") do io
        proc = run(pipeline(detach(cmd); stdout=io, stderr=io); wait=false)
        t0 = time()
        while process_running(proc)
            if time() - t0 > limit_s
                println(stderr, "TIMEOUT after $(round(Int, limit_s)) s: killing process group")
                if Sys.isunix()
                    ccall(:kill, Cint, (Cint, Cint), -getpid(proc), 9)   # whole group
                end
                kill(proc, Base.SIGKILL)
                wait(proc)
                return (:timeout, time() - t0)
            end
            sleep(0.5)
        end
        return (success(proc) ? :pass : :fail, time() - t0)
    end
end

function run_shard(shard::String)
    julia = joinpath(Sys.BINDIR, Base.julia_exename())
    cov = get(ENV, "ADRIA_TEST_COVERAGE", "none")
    cmd = setenv(
        `$julia --project=$ENV_DIR --code-coverage=$cov $(joinpath(TEST_DIR, "runtests.jl"))`,
        merge(Dict(ENV), Dict("ADRIA_TEST_SHARD" => shard))
    )
    return run_with_timeout(cmd, limit_for(shard), joinpath(LOG_DIR, "$shard.log"))
end

function main(args)
    prepare = !("--no-prepare" in args)
    shards = filter(a -> !startswith(a, "--"), args)
    isempty(shards) && (shards = ["fast", "spec", "sim"])
    for s in shards
        haskey(DEFAULT_LIMITS, s) || error("Unknown shard $(repr(s)); choose from $(collect(keys(DEFAULT_LIMITS)))")
    end

    mkpath(LOG_DIR)
    prepare && prepare_env()

    tasks = Dict(s => Threads.@spawn(run_shard(s)) for s in shards)
    results = Dict(s => fetch(t) for (s, t) in tasks)

    println("\n", "="^60, "\nADRIA test shards")
    for s in shards
        status, secs = results[s]
        println(rpad(s, 6), rpad(uppercase(String(status)), 9), round(secs; digits=0), " s  (limit ",
            round(Int, limit_for(s)), " s)  log: ", joinpath(LOG_DIR, "$s.log"))
    end
    println("="^60)

    bad = [s for s in shards if results[s][1] != :pass]
    for s in bad  # show the tail of failing logs so CI output is useful
        lines = readlines(joinpath(LOG_DIR, "$s.log"))
        println("\n--- last lines of $s.log ---")
        foreach(println, lines[max(1, end - 40):end])
    end
    exit(isempty(bad) ? 0 : 1)
end

main(ARGS)
