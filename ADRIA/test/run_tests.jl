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
# Output goes to ADRIA/test/logs/<shard>.log. (It is not echoed live; the tail of the log
# of any failing shard is printed at the end.) If the launcher is interrupted or crashes,
# all running shard process groups are killed.
#
# Shards run like `Pkg.test` does: cwd = ADRIA/test (so test/config.toml is found) and the
# same julia flags (--check-bounds=yes etc). Coverage is controlled by
# ADRIA_TEST_COVERAGE: "none" (default), or "user"/"adria" to track ADRIA's own sources
# only (`--code-coverage=@<ADRIA dir>`; plain `user` would discard precompiled code of
# all dependencies), or any other value passed to `--code-coverage` as-is.

using Pkg

const TEST_DIR = abspath(@__DIR__)
const ADRIA_DIR = dirname(TEST_DIR)
const ENV_DIR = joinpath(TEST_DIR, ".testenv")   # shared, prepared once, git-ignored
const LOG_DIR = joinpath(TEST_DIR, "logs")

# Default wall-clock limits in seconds, measured per shard process (the launcher's own
# environment preparation is not counted). `fast` and `spec` are 6 min as a "something is
# wrong" tripwire; `sim` runs 32 full simulations and needs longer. NB: in CI the job-level
# `timeout-minutes` (UnitTests.yml) also covers checkout/cache/precompile, so for `fast`
# and `spec` the job timeout, not this limit, is the effective bound.
const DEFAULT_LIMITS = Dict("fast" => 360, "spec" => 360, "sim" => 1200)

limit_for(shard) = parse(
    Float64, get(ENV, "ADRIA_TEST_LIMIT_$(uppercase(shard))_S", string(DEFAULT_LIMITS[shard]))
)

"""Create/refresh the shared test environment once, so concurrent shards don't each
instantiate and precompile their own copy.

Note: dependency versions are resolved fresh, NOT taken from ADRIA's committed
`Manifest-v1.12.toml`. Copying that manifest into the test environment was tried: it is
the manifest of a different project (no test-only deps), and `Pkg.develop` re-resolves
and upgrades many packages anyway, so it does not pin versions. The julia-actions cache
keyed on the manifest still speeds up CI."""
function prepare_env()
    mkpath(ENV_DIR)
    cp(joinpath(TEST_DIR, "Project.toml"), joinpath(ENV_DIR, "Project.toml"); force=true)
    @info "Preparing shared test environment" ENV_DIR
    Pkg.activate(ENV_DIR) do
        Pkg.develop(; path=ADRIA_DIR)
        Pkg.instantiate()
    end
    # Precompile in a subprocess with the same flags the shards use: `--check-bounds` is
    # part of the precompile cache key, so caches built by an in-process
    # `Pkg.precompile()` (check-bounds=auto) would be rejected by shards started with
    # `--check-bounds=yes`, and every shard would re-precompile ADRIA and its deps itself
    # (inside its time limit, and racing the other shards).
    julia = joinpath(Sys.BINDIR, Base.julia_exename())
    run(`$julia --project=$ENV_DIR $JULIA_FLAGS -e "using Pkg; Pkg.precompile()"`)
end

# Same flags as `Pkg.test`. Only `--check-bounds` affects the precompile cache key.
const JULIA_FLAGS = [
    "--check-bounds=yes", "--warn-overwrite=yes", "--depwarn=yes", "--startup-file=no"
]

const RUNNING = Set{Base.Process}()
const RUNNING_LOCK = ReentrantLock()

function kill_group(proc)
    if Sys.isunix()
        ccall(:kill, Cint, (Cint, Cint), -getpid(proc), 9)   # whole group
    elseif Sys.iswindows()
        # No process groups: `taskkill /T` kills the whole process tree.
        run(ignorestatus(`taskkill /F /T /PID $(getpid(proc))`))
    end
    kill(proc, Base.SIGKILL)
end

function kill_all_running()
    lock(RUNNING_LOCK) do
        for p in RUNNING
            process_running(p) && try
                kill_group(p)
            catch
            end
        end
    end
end

function run_with_timeout(cmd::Cmd, limit_s::Real, logfile::AbstractString)
    open(logfile, "w") do io
        proc = run(pipeline(detach(cmd); stdout=io, stderr=io); wait=false)
        lock(() -> push!(RUNNING, proc), RUNNING_LOCK)
        t0 = time()
        try
            while process_running(proc)
                if time() - t0 > limit_s
                    println(stderr, "TIMEOUT after $(round(Int, limit_s)) s: killing process group")
                    kill_group(proc)
                    wait(proc)
                    return (:timeout, time() - t0)
                end
                sleep(0.5)
            end
            return (success(proc) ? :pass : :fail, time() - t0)
        finally
            process_running(proc) && try
                kill_group(proc)
            catch
            end
            lock(() -> delete!(RUNNING, proc), RUNNING_LOCK)
        end
    end
end

function coverage_flag(cov::AbstractString)
    cov in ("user", "adria") && return "--code-coverage=@$(ADRIA_DIR)"
    return "--code-coverage=$cov"
end

function run_shard(shard::String)
    julia = joinpath(Sys.BINDIR, Base.julia_exename())
    cov = coverage_flag(get(ENV, "ADRIA_TEST_COVERAGE", "none"))
    runtests = joinpath(TEST_DIR, "runtests.jl")
    cmd = setenv(
        `$julia --project=$ENV_DIR $cov $JULIA_FLAGS $runtests`,
        merge(Dict(ENV), Dict("ADRIA_TEST_SHARD" => shard));
        dir=TEST_DIR
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

    # Tasks (not threads) are enough: they spend their time in `sleep`. Cleanup on
    # Ctrl-C / error / normal `exit` kills any still-running shard process groups.
    atexit(kill_all_running)
    results = try
        tasks = Dict(s => Threads.@spawn(run_shard(s)) for s in shards)
        Dict(s => fetch(t) for (s, t) in tasks)
    catch
        kill_all_running()
        rethrow()
    end

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
