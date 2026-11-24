# Logic of the opt-in baseline regression harness (ADRIA/test/regression/): the comparison
# and its exit codes, the META rules and the runner's argument handling, on small synthetic
# Zarr stores. No git, no driver and no reference build (the harness itself is run by hand,
# see regression/README.md).

"""The comparison script, included as a module (its entry point only runs as a script)."""
module RegressionCompare
using ADRIA: Zarr, JSON
using .Zarr, .JSON
include(joinpath(@__DIR__, "regression", "compare.jl"))
end

"""The runner, included as a module (its entry point only runs as a script)."""
module RegressionRunner
include(joinpath(@__DIR__, "regression", "run_regression.jl"))
end

@testset "regression harness" begin
    RC = RegressionCompare
    RR = RegressionRunner

    META = Dict(
        "julia_version" => "1.12.7", "threads" => "4", "driver" => "driver.jl",
        "driver_sha256" => "abc", "source" => "rev", "os" => "Linux", "arch" => "x86_64"
    )

    """Write a synthetic run folder: one case with one store and a META."""
    function write_run(
        dir; cover=[1.0 2.0; 3.0 NaN], ranks=UInt16[1, 2, 3], meta=META, case="c1"
    )
        mkpath(dir)
        store = RC.Zarr.zgroup(joinpath(dir, case, "store.zarr"))
        results = RC.Zarr.zgroup(store, "results")
        z = RC.Zarr.zcreate(Float64, results, "cover", 2, 2; chunks=(2, 2))
        z[:, :] = cover
        logs = RC.Zarr.zgroup(store, "logs")
        r = RC.Zarr.zcreate(UInt16, logs, "rankings", 3; chunks=(3,))
        r[:] = ranks
        if !isnothing(meta)
            write(joinpath(dir, "META"), join(["$k=$v" for (k, v) in meta], "\n") * "\n")
        end
        return dir
    end

    """Call `f()` with stdout and stderr discarded."""
    function quiet(f)
        return redirect_stdout(devnull) do
            redirect_stderr(devnull) do
                f()
            end
        end
    end

    """The exit code of comparing `cand` with `ref`, output suppressed."""
    code(ref, cand, extra=String[]) = quiet(() -> RC.run_main([ref, cand, extra...]))

    mktempdir() do tmp
        ref = write_run(joinpath(tmp, "ref"))

        @testset "comparison and exit codes" begin
            @test code(ref, write_run(joinpath(tmp, "same"))) == 0

            # NaN-aware: NaN == NaN in the same position
            @test code(ref, write_run(joinpath(tmp, "nan"); cover=[1.0 2.0; 3.0 NaN])) == 0

            d = write_run(joinpath(tmp, "value"); cover=[1.0 2.0; 3.5 NaN])
            @test code(ref, d) == 1
            @test code(ref, d, ["--tol", "1.0"]) == 0
            @test code(ref, d, ["--tol", "0.1"]) == 1

            # NaN moved: differs even under --tol
            d = write_run(joinpath(tmp, "nanmoved"); cover=[1.0 NaN; 3.0 2.0])
            @test code(ref, d, ["--tol", "1e3"]) == 1

            # Ranks are always exact
            d = write_run(joinpath(tmp, "ranks"); ranks=UInt16[1, 2, 4])
            @test code(ref, d) == 1
            @test code(ref, d, ["--tol", "1e3"]) == 1

            # A different case set
            d = write_run(joinpath(tmp, "othercase"); case="c2")
            @test code(ref, d) == 1

            # An unexpected layout is a crash (3), not a difference
            d = write_run(joinpath(tmp, "twostores"))
            cp(joinpath(d, "c1", "store.zarr"), joinpath(d, "c1", "other.zarr"))
            @test code(ref, d) == 3

            # Usage errors
            @test code(ref, ref, ["--tol"]) == 2
            @test code(ref, ref, ["--tol", "x"]) == 2
            @test quiet(() -> RC.main([ref])) == 2
        end

        @testset "flat and grouped submodel log layouts match by contents" begin
            """A run folder with a COTS log, flat (`logs/cots_pop_log`) or grouped."""
            function with_log(name; grouped, vals=Float32[1, 2])
                d = write_run(joinpath(tmp, name))
                logs = RC.Zarr.zopen(joinpath(d, "c1", "store.zarr", "logs"), "w")
                grp = grouped ? RC.Zarr.zgroup(logs, "cots") : logs
                z = RC.Zarr.zcreate(
                    Float32, grp, grouped ? "pop" : "cots_pop_log", 2; chunks=(2,)
                )
                z[:] = vals
                return d
            end

            flat = with_log("flat"; grouped=false)
            @test code(flat, with_log("grouped"; grouped=true)) == 0
            other = with_log("grouped_diff"; grouped=true, vals=Float32[1, 3])
            @test code(flat, other) == 1

            # A log missing from the candidate is still a difference
            @test code(flat, write_run(joinpath(tmp, "no_log"))) == 1

            # A store holding a log both flat and grouped is a difference, not a pick
            both = with_log("both"; grouped=true)
            logs = RC.Zarr.zopen(joinpath(both, "c1", "store.zarr", "logs"), "w")
            z = RC.Zarr.zcreate(Float32, logs, "cots_pop_log", 2; chunks=(2,))
            z[:] = Float32[1, 2]
            @test code(flat, both) == 1
            clashes = String[]
            RC.find_arrays(joinpath(both, "c1", "store.zarr"), clashes)
            @test startswith(only(clashes), "two arrays map to logs/cots/pop: ")
            @test occursin("logs/cots_pop_log", only(clashes))

            @test RC.canonical_path("logs/cots_condition_log") == "logs/cots/condition"
            @test RC.canonical_path("logs/rankings") == "logs/rankings"
            @test RC.is_submodel_group("logs/cots") && !RC.is_submodel_group("logs")
        end

        @testset "the failure names the differing array" begin
            d = write_run(joinpath(tmp, "named"); cover=[9.0 2.0; 3.0 NaN])
            out = joinpath(tmp, "out.txt")
            open(out, "w") do io
                redirect_stdout(io) do
                    RC.main([ref, d])
                end
            end
            @test occursin("c1/results/cover", read(out, String))
        end

        @testset "META rules" begin
            """A run folder whose META has the given `key => value` pairs changed."""
            function with(name, kv::Pair...)
                return write_run(joinpath(tmp, name); meta=merge(META, Dict(kv...)))
            end

            for (key, other) in
                (("julia_version", "0.0.1"), ("threads", "123"), ("driver_sha256", "0"))
                d = with("m_$key", key => other)
                @test code(ref, d) == 2
                @test code(ref, d, ["--tol", "1e-4"]) == 0
            end

            # A driver change is accepted in exact mode only with --driver-variant
            d = with("m_variant", "driver_sha256" => "0")
            @test code(ref, d, ["--driver-variant"]) == 0

            # A different OS or architecture is refused even under --tol
            for (key, other) in (("os", "Windows"), ("arch", "aarch64"))
                d = with("m_$key", key => other)
                @test code(ref, d) == 2
                @test code(ref, d, ["--tol", "1e-4"]) == 2
                @test code(ref, d, ["--driver-variant"]) == 2
            end

            # Missing or empty values never match, not even two empty ones
            d = write_run(joinpath(tmp, "m_none"); meta=nothing)
            @test code(ref, d) == 2
            @test code(ref, d, ["--tol", "1e-4"]) == 0
            @test code(d, ref) == 2
            d = with("m_empty", "julia_version" => "")
            @test code(d, d) == 2
            d = write_run(
                joinpath(tmp, "m_noarch"); meta=filter(p -> first(p) != "arch", META)
            )
            @test code(ref, d) == 2

            # Padded keys and values and CRLF line ends are trimmed
            d = write_run(joinpath(tmp, "m_pad"))
            write(
                joinpath(d, "META"),
                join([" $k = $v \r" for (k, v) in META], "\n") * "\n"
            )
            @test code(ref, d) == 0
        end
    end

    @testset "driver checksum" begin
        mktempdir() do tmp
            lf = joinpath(tmp, "lf.jl")
            crlf = joinpath(tmp, "crlf.jl")
            other = joinpath(tmp, "other.jl")
            write(lf, "a = 1\nb = 2\n")
            write(crlf, "a = 1\r\nb = 2\r\n")
            write(other, "a = 1\nb = 3\n")
            @test RR.driver_sha256(lf) == RR.driver_sha256(crlf)
            @test RR.driver_sha256(lf) != RR.driver_sha256(other)
            @test length(RR.driver_sha256(lf)) == 64

            # A lone CR is not a line end
            write(other, "a = 1\rb = 2\n")
            @test RR.driver_sha256(lf) != RR.driver_sha256(other)
        end
    end

    @testset "runner arguments and refusals" begin
        parse = RR.parse_compare_args
        @test parse(String[]) == (; pos=String[], extra=String[])
        @test parse(["--tol", "1e-4", "co", "rev"]) ==
            (; pos=["co", "rev"], extra=["--tol", "1e-4"])
        @test parse(["co", "--driver-variant", "rev", "--tol", "2"]) ==
            (; pos=["co", "rev"], extra=["--driver-variant", "--tol", "2"])
        @test_throws RR.Refusal parse(["--tol"])
        @test_throws RR.Refusal parse(["--tol", "x"])
        @test_throws RR.Refusal parse(["--nope"])
        @test_throws RR.Refusal parse(["a", "b", "c"])

        # Every refusal is exit code 2
        @test quiet(() -> RR.main(String[])) == 2
        @test quiet(() -> RR.main(["bogus"])) == 2
        @test quiet(() -> RR.main(["compare", "--tol"])) == 2
        @test quiet(() -> RR.main(["compare-dirs", "only-one"])) == 2
        @test quiet(() -> RR.main(["ref", "a", "b"])) == 2
        @test quiet(() -> RR.main(["selftest", "x"])) == 2

        cfg = RR.load_config(
            Dict("ADRIA_REGRESSION_REF" => "r", "ADRIA_REGRESSION_THREADS" => "2")
        )
        @test cfg.threads == 2
        @test cfg.cache == abspath("r")
        @test endswith(RR.load_config(Dict()).driver, "driver.jl")
        @test RR.load_config(Dict()).threads == 4
        @test RR.load_config(Dict()).cache ==
            joinpath(homedir(), ".cache", "adria-regression")
        @test_throws RR.Refusal RR.load_config(Dict("ADRIA_REGRESSION_THREADS" => "0"))
        @test_throws RR.Refusal RR.load_config(Dict("ADRIA_REGRESSION_THREADS" => "x"))
    end

    @testset "runner META handling" begin
        mktempdir() do tmp
            driver = joinpath(tmp, "driver.jl")
            write(driver, "# driver\n")
            cfg = RR.Config(joinpath(tmp, "cache"), driver, 4)
            want = RR.expected_meta(cfg)
            @test want["julia_version"] == string(VERSION)
            @test want["threads"] == "4"
            @test want["os"] == string(Sys.KERNEL)
            @test want["arch"] == string(Sys.ARCH)

            """A fresh run folder as the driver writes it."""
            function driver_out(name; threads=4)
                d = joinpath(tmp, name)
                mkpath(d)
                write(joinpath(d, "META"), "julia_version=$VERSION\nthreads=$threads\n")
                return d
            end

            # The runner completes the META; the result is a valid reference
            d = driver_out("out")
            RR.write_meta(cfg, d, "somerev")
            m = RR.read_meta(d)
            @test isempty(RR.meta_lacking(m))
            @test m["driver"] == "driver.jl"
            @test m["source"] == "somerev"
            @test isempty(RR.meta_differing(m, want))
            @test RR.check_cached(cfg, d) === nothing

            # Refusals (exit code 2): wrong thread count, no META, an empty source
            @test_throws RR.Refusal RR.write_meta(cfg, driver_out("t"; threads=3), "rev")
            mkpath(joinpath(tmp, "nometa"))
            @test_throws RR.Refusal RR.write_meta(cfg, joinpath(tmp, "nometa"), "rev")
            @test_throws RR.Refusal RR.write_meta(cfg, driver_out("e"), "")

            # A cached reference is reused only if every condition matches
            @test_throws RR.Refusal RR.check_cached(cfg, joinpath(tmp, "nometa"))
            for (k, v) in
                (("julia_version", "0.0.1"), ("threads", "9"), ("driver_sha256", "0"),
                ("os", "Plan9"), ("arch", "riscv64"), ("source", ""))
                bad = deepcopy(m)
                bad[k] = v
                dir = joinpath(tmp, "bad_$k")
                mkpath(dir)
                write(joinpath(dir, "META"), join(["$a=$b" for (a, b) in bad], "\n") * "\n")
                @test_throws RR.Refusal RR.check_cached(cfg, dir)
            end

            # A reference without os/arch (built before they were recorded) is refused
            old = deepcopy(m)
            delete!(old, "os")
            delete!(old, "arch")
            dir = joinpath(tmp, "old")
            mkpath(dir)
            write(joinpath(dir, "META"), join(["$a=$b" for (a, b) in old], "\n") * "\n")
            @test_throws RR.Refusal RR.check_cached(cfg, dir)
        end
    end
end
