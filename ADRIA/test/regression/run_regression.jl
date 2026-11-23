# Opt-in baseline regression runner (see README.md). Not part of the test suite.
#
#   julia run_regression.jl ref [REV]            build + cache the reference for REV (default HEAD)
#   julia run_regression.jl compare [CHECKOUT] [REV] [--tol X] [--driver-variant]
#                                                (options may come first)
#                                                run the driver on CHECKOUT (default: this repo)
#                                                and compare against the cached reference of REV
#                                                (default: the one recorded by the last `ref`)
#   julia run_regression.jl selftest             reference-vs-itself passes; perturbations fail
#   julia run_regression.jl compare-dirs REF_DIR CAND_DIR [--tol X] [--driver-variant]
#                                                compare two existing run folders (the self-test
#                                                uses this; no driver is run)
#
# Exit codes: 0 pass, 1 differences, 2 refusal (usage error, mismatched conditions, missing
# or empty META values, no reference), 3 crash (the driver, compare.jl or the runner failed).
#
# Env: ADRIA_REGRESSION_REF (cache dir, default ~/.cache/adria-regression),
#      ADRIA_REGRESSION_DRIVER (driver script, default driver.jl next to this file),
#      ADRIA_REGRESSION_THREADS (default 4), TMPDIR/TEMP (scratch).
#
# The runner runs under the Julia that launched it, and starts the driver and compare.jl
# with that same Julia (`Base.julia_cmd()`). It uses the standard library only; the
# comparison runs in a child process in the checkout's ADRIA environment (Zarr, JSON).

include(joinpath(@__DIR__, "meta.jl"))

"""The folder holding this script, the driver and compare.jl."""
const HERE = @__DIR__

"""The ADRIA project of the checkout this script lives in (`<repo>/ADRIA`)."""
const ADRIA_ENV = normpath(joinpath(HERE, "..", ".."))

"""Exit code: identical (or within `--tol`)."""
const PASS = 0

"""Exit code: differences found."""
const DIFF = 1

"""Exit code: refusal (usage error, mismatched conditions, missing or empty META values)."""
const REFUSED = 2

"""Exit code: crash (the driver, compare.jl or the runner failed)."""
const CRASH = 3

"""A condition that makes the runner refuse to go ahead (exit code 2), with its message."""
struct Refusal <: Exception
    msg::String
end

"""A driver run, or the comparison, that failed (exit code 3), with its message."""
struct Crash <: Exception
    msg::String
end

"""
    die(msg)

Refuse to go ahead with `msg` (exit code 2).
"""
die(msg) = throw(Refusal(msg))

"""Settings from the environment: the reference cache, the driver script, the threads."""
struct Config
    cache::String
    driver::String
    threads::Int
end

"""
    load_config(env=ENV)::Config

The settings of `ADRIA_REGRESSION_REF`, `ADRIA_REGRESSION_DRIVER` and
`ADRIA_REGRESSION_THREADS`, with their defaults.
"""
function load_config(env=ENV)::Config
    cache = get(
        env, "ADRIA_REGRESSION_REF", joinpath(homedir(), ".cache", "adria-regression")
    )
    driver = get(env, "ADRIA_REGRESSION_DRIVER", joinpath(HERE, "driver.jl"))
    threads = tryparse(Int, get(env, "ADRIA_REGRESSION_THREADS", "4"))
    (isnothing(threads) || threads < 1) &&
        die("ADRIA_REGRESSION_THREADS must be a positive integer")
    return Config(abspath(cache), abspath(driver), threads)
end

"""
    expected_meta(cfg)::Dict{String,String}

What a reference must record to be reusable here: this Julia version, the thread count, the
driver checksum, and this OS and architecture.
"""
function expected_meta(cfg::Config)::Dict{String,String}
    return merge(
        Dict(
            "julia_version" => string(VERSION),
            "threads" => string(cfg.threads),
            "driver_sha256" => driver_sha256(cfg.driver)
        ),
        host_meta()
    )
end

"""
    repo_root()::String

The git work tree this script lives in.
"""
function repo_root()::String
    root = git_read(HERE, "rev-parse", "--show-toplevel")
    isnothing(root) && die("$HERE is not inside a git work tree")
    return root
end

"""
    git_read(dir, args...)::Union{Nothing,String}

The trimmed output of `git -C dir args...`, or `nothing` if git fails.
"""
function git_read(dir, args...)::Union{Nothing,String}
    cmd = pipeline(`git -C $dir $(collect(args))`; stderr=devnull)
    return try
        strip(read(cmd, String))
    catch
        nothing
    end
end

"""
    resolve_rev(repo, rev)::String

The full commit sha of `rev` (a sha, tag or branch), resolved in `repo`.
"""
function resolve_rev(repo, rev)::String
    sha = git_read(repo, "rev-parse", "--verify", "--quiet", "$rev^{commit}")
    (isnothing(sha) || isempty(sha)) && die("cannot resolve '$rev' to a commit in $repo")
    return sha
end

"""
    current_rev(cfg)::String

The rev recorded by the last `ref` (the cache's `CURRENT` file).
"""
function current_rev(cfg::Config)::String
    f = joinpath(cfg.cache, "CURRENT")
    isfile(f) ||
        die(
            "no reference recorded in $(cfg.cache) (no CURRENT file); run: run_regression.jl ref <rev>"
        )
    rev = strip(read(f, String))
    isempty(rev) && die("$f is empty; run: run_regression.jl ref <rev>")
    return rev
end

"""
    require_ref(cfg, rev)::String

The folder of the cached reference of `rev`; refuses if there is none.
"""
function require_ref(cfg::Config, rev)::String
    dir = joinpath(cfg.cache, rev)
    isdir(dir) ||
        die("no cached reference for $rev in $(cfg.cache); run: run_regression.jl ref $rev")
    return dir
end

"""
    remove_tree(path)

Remove a folder, retrying a few times (on Windows an antivirus scanner or a lingering
handle can hold files briefly). Never throws.
"""
function remove_tree(path)
    for attempt in 1:4
        try
            rm(path; recursive=true, force=true)
        catch
        end
        isdir(path) || return nothing
        sleep(0.25 * attempt)
    end
    return nothing
end

"""
    remove_worktree(repo, wt)

Remove the git worktree `wt` of `repo`. Every child process has exited by now (`run` waits
for them), so the folder is not in use; if `git worktree remove` fails the folder is
deleted directly and the stale worktree record pruned.
"""
function remove_worktree(repo, wt)
    ok = success(pipeline(`git -C $repo worktree remove --force $wt`; stderr=devnull))
    if !ok || isdir(wt)
        remove_tree(wt)
        success(pipeline(`git -C $repo worktree prune`; stderr=devnull))
    end
    return nothing
end

"""
    meta_get(dir, key)

The value of `key` in `<dir>/META` (trimmed), or the empty string.
"""
function meta_get(dir, key)
    m = read_meta(dir)
    return isnothing(m) ? "" : get(m, key, "")
end

"""
    write_meta(cfg, outdir, source)

Complete the META the driver wrote (its own Julia version and thread count) with the
driver, its checksum, the source, and this OS and architecture. Refuses if the driver wrote
no META, if any value is empty, or if the driver ran with another thread count.
"""
function write_meta(cfg::Config, outdir, source)
    isfile(joinpath(outdir, "META")) || die("the driver wrote no $outdir/META")
    open(joinpath(outdir, "META"), "a") do io
        println(io, "driver=", basename(cfg.driver))
        println(io, "driver_sha256=", driver_sha256(cfg.driver))
        println(io, "source=", source)
        for (k, v) in host_meta()
            println(io, k, "=", v)
        end
    end
    lacking = meta_lacking(read_meta(outdir))
    isempty(lacking) || die("empty or missing $(join(lacking, ", ")) in $outdir/META")
    got = meta_get(outdir, "threads")
    got == string(cfg.threads) ||
        die("the driver ran with $got threads, not $(cfg.threads)")
    return nothing
end

"""
    run_driver(cfg, checkout, outdir, source, scratch)

Run the driver on the ADRIA project of `checkout`, writing the results and META to
`outdir`. Julia is started from an empty folder in `scratch`, so no config.toml or Outputs
is picked up from (or written to) the current directory. Throws `Crash` if it fails.
"""
function run_driver(cfg::Config, checkout, outdir, source, scratch)
    env = joinpath(checkout, "ADRIA")
    cwd = joinpath(scratch, "cwd")
    mkpath(cwd)
    cmd = Cmd(
        `$(Base.julia_cmd()) --project=$env -t $(cfg.threads) $(cfg.driver) $env $outdir`;
        dir=cwd
    )
    proc = run(ignorestatus(cmd))
    success(proc) ||
        throw(Crash("the driver failed on $checkout (exit code $(proc.exitcode))"))
    write_meta(cfg, outdir, source)
    return nothing
end

"""
    compare_dirs(ref, cand, extra=String[]; env=ADRIA_ENV)::Int

Run compare.jl on two run folders in the ADRIA environment `env` (Zarr and JSON come from
there) and return its exit code (0 pass, 1 differences, 2 refusal, 3 crash). A compare.jl
killed by a signal counts as a crash.
"""
function compare_dirs(ref, cand, extra=String[]; env=ADRIA_ENV)::Int
    cmd = `$(Base.julia_cmd()) --project=$env $(joinpath(HERE, "compare.jl")) $ref $cand $extra`
    proc = run(ignorestatus(cmd))
    proc.termsignal == 0 || return CRASH
    return Int(proc.exitcode)
end

"""
    parse_compare_args(args)::NamedTuple

Split the arguments of `compare` / `compare-dirs` into the positionals and the options for
compare.jl (`--tol X`, `--driver-variant`), which may come before, between or after them.
Refuses on an unknown option, a missing or non-numeric `--tol` value, or too many
positionals.
"""
function parse_compare_args(args, max_positional::Int=2)
    pos = String[]
    extra = String[]
    i = 1
    while i <= length(args)
        a = args[i]
        if a == "--tol"
            i < length(args) || die("--tol needs a value")
            isnothing(tryparse(Float64, args[i + 1])) && die("--tol needs a number")
            push!(extra, a, args[i + 1])
            i += 2
        elseif a == "--driver-variant"
            push!(extra, a)
            i += 1
        elseif startswith(a, "--")
            die("unknown option $a")
        else
            length(pos) < max_positional || die("too many arguments: $a")
            push!(pos, a)
            i += 1
        end
    end
    return (; pos, extra)
end

"""
    cmd_ref(cfg, args)::Int

`ref [REV]`: build and cache the reference for `REV` (default HEAD) from a throwaway git
worktree, or check that an already cached one matches the current Julia version, thread
count, driver and OS/architecture. Records `REV` as the current reference.
"""
function cmd_ref(cfg::Config, args)::Int
    length(args) <= 1 || die("too many arguments: $(args[2])")
    repo = repo_root()
    rev = resolve_rev(repo, isempty(args) ? "HEAD" : args[1])
    dest = joinpath(cfg.cache, rev)
    if isdir(dest)
        println("reference $rev already cached at $dest")
        check_cached(cfg, dest)
    else
        build_reference(cfg, repo, rev, dest)
    end
    mkpath(cfg.cache)
    write(joinpath(cfg.cache, "CURRENT"), rev * "\n")
    return PASS
end

"""
    check_cached(cfg, dest)

Refuse to reuse the cached reference `dest` if it has no META, a META with a missing or
empty value, or one that does not match the current Julia version, thread count, driver
checksum, OS and architecture.
"""
function check_cached(cfg::Config, dest)
    fix = "Remove $dest to rebuild it, or use another ADRIA_REGRESSION_REF."
    m = read_meta(dest)
    isnothing(m) && die(
        "the cached reference has no META (built before META was recorded), so its " *
        "Julia version, threads, driver and OS cannot be checked. $fix"
    )
    lacking = meta_lacking(m)
    isempty(lacking) ||
        die("the cached reference's META lacks a value for $(join(lacking, ", ")). $fix")
    want = expected_meta(cfg)
    differing = meta_differing(m, want)
    isempty(differing) && return nothing
    have = join(["$k=$(m[k])" for k in differing], ", ")
    need = join(["$k=$(want[k])" for k in differing], ", ")
    die("the cached reference was built with $have; now $need. $fix")
end

"""
    build_reference(cfg, repo, rev, dest)

Check out `rev` in a throwaway worktree, run the driver on it, and move the results to
`dest`. The worktree and the scratch folder are removed afterwards, also on failure.
"""
function build_reference(cfg::Config, repo, rev, dest)
    scratch = mktempdir(; prefix="adr-")
    wt = joinpath(scratch, "w")
    added = false
    try
        run(`git -C $repo -c core.longpaths=true worktree add --detach $wt $rev`)
        added = true
        mkpath(cfg.cache)
        run_driver(cfg, wt, joinpath(scratch, "out"), rev, scratch)
        mv(joinpath(scratch, "out"), dest)
        println("cached reference $rev at $dest")
    finally
        added && remove_worktree(repo, wt)
        remove_tree(scratch)
    end
    return nothing
end

"""
    cmd_compare(cfg, args)::Int

`compare [CHECKOUT] [REV] [--tol X] [--driver-variant]`: run the driver on `CHECKOUT`
(default: this repository) and compare it against the cached reference of `REV` (default:
the one recorded by the last `ref`). Returns compare.jl's exit code.
"""
function cmd_compare(cfg::Config, args)::Int
    (; pos, extra) = parse_compare_args(args)
    repo = repo_root()
    cand = isempty(pos) ? repo : pos[1]
    rev = length(pos) >= 2 ? resolve_rev(repo, pos[2]) : current_rev(cfg)
    ref = require_ref(cfg, rev)
    isdir(joinpath(cand, "ADRIA")) ||
        die("$cand is not an ADRIA checkout (no ADRIA/ folder)")
    cand = abspath(cand)
    source = "$cand at $(candidate_rev(cand))"

    scratch = mktempdir(; prefix="adr-")
    try
        out = joinpath(scratch, "cand")
        run_driver(cfg, cand, out, source, scratch)
        return compare_dirs(ref, out, extra; env=joinpath(cand, "ADRIA"))
    finally
        remove_tree(scratch)
    end
end

"""
    candidate_rev(checkout)::String

The HEAD sha of `checkout`, with `+dirty` appended if tracked files are modified, or
`unknown` if it is not a git work tree.
"""
function candidate_rev(checkout)::String
    sha = git_read(checkout, "rev-parse", "HEAD")
    (isnothing(sha) || isempty(sha)) && return "unknown"
    status = git_read(checkout, "status", "--porcelain", "--untracked-files=no")
    return isnothing(status) || isempty(status) ? sha : sha * "+dirty"
end

"""
    cmd_compare_dirs(args)::Int

`compare-dirs REF_DIR CAND_DIR [--tol X] [--driver-variant]`: compare two existing run
folders without running the driver. Returns compare.jl's exit code.
"""
function cmd_compare_dirs(args)::Int
    (; pos, extra) = parse_compare_args(args)
    length(pos) == 2 || die(
        "usage: run_regression.jl compare-dirs REF_DIR CAND_DIR [--tol X] [--driver-variant]"
    )
    return compare_dirs(pos[1], pos[2], extra)
end

"""
    cmd_selftest(cfg, args)::Int

`selftest`: the cached reference must pass against itself, then selftest.jl checks that
perturbed copies of it fail with the right exit codes.
"""
function cmd_selftest(cfg::Config, args)::Int
    isempty(args) || die("selftest takes no arguments")
    ref = require_ref(cfg, current_rev(cfg))
    println("== reference vs itself (must pass)")
    compare_dirs(ref, ref) == PASS || return DIFF
    cmd = `$(Base.julia_cmd()) --project=$ADRIA_ENV $(joinpath(HERE, "selftest.jl")) $ref $(@__FILE__)`
    proc = run(ignorestatus(cmd))
    return proc.termsignal == 0 ? Int(proc.exitcode) : CRASH
end

"""The usage text printed on a missing or unknown command."""
const USAGE = """
usage: julia run_regression.jl ref [REV]
       julia run_regression.jl compare [CHECKOUT] [REV] [--tol X] [--driver-variant]
       julia run_regression.jl selftest
       julia run_regression.jl compare-dirs REF_DIR CAND_DIR [--tol X] [--driver-variant]
exit codes: 0 pass, 1 differences, 2 refusal, 3 crash (see README.md)"""

"""
    main(args)::Int

Dispatch on the subcommand and return the exit code. Refusals (`die`) print their message
and return 2; a `Crash` or any other exception returns 3.
"""
function main(args)::Int
    return try
        isempty(args) && die(USAGE)
        cmd, rest = args[1], args[2:end]
        if cmd == "ref"
            cmd_ref(load_config(), rest)
        elseif cmd == "compare"
            cmd_compare(load_config(), rest)
        elseif cmd == "selftest"
            cmd_selftest(load_config(), rest)
        elseif cmd == "compare-dirs"
            cmd_compare_dirs(rest)
        else
            die("unknown command '$cmd'\n" * USAGE)
        end
    catch err
        if err isa Refusal
            println(stderr, "ERROR: ", err.msg)
            REFUSED
        elseif err isa Crash
            println(stderr, "ERROR: ", err.msg)
            CRASH
        else
            println(stderr, "ERROR: run_regression.jl failed:")
            showerror(stderr, err, catch_backtrace())
            println(stderr)
            CRASH
        end
    end
end

abspath(PROGRAM_FILE) == (@__FILE__) && exit(main(ARGS))
