# Self-test for compare.jl: perturbed copies of a reference must FAIL, the unperturbed copy
# (and changes compare.jl is meant to accept) must PASS.
#   julia --project=<ADRIA env> selftest.jl <ref_dir> <compare.jl>
using Zarr
using JSON

src_ref, cmp = ARGS
julia = Base.julia_cmd()
project = dirname(Base.active_project())
# The reference the candidates are compared against: a copy of `src_ref` with a META file
const REF = Ref("")

passes(cand, extra=``) =
    run(
        pipeline(
            ignorestatus(`$julia --project=$project $cmp $(REF[]) $cand $extra`);
            stdout=devnull, stderr=devnull
        )
    ).exitcode == 0

store(root, c) = only(filter(isdir, readdir(joinpath(root, c); join=true)))
results = Bool[]
function expect(name, cand, want, extra=``)
    ok = passes(cand, extra) == want
    println(ok ? "ok   " : "FAIL ", name, " (expected ", want ? "pass" : "fail", ")")
    push!(results, ok)
end

mktempdir() do tmp
    REF[] = joinpath(tmp, "ref")
    cp(src_ref, REF[])
    if !isfile(joinpath(REF[], "META"))
        # A reference built before META was recorded: give the copy one, so that the META
        # cases below test the check rather than the "no META" warning
        println("note: reference has no META; using a synthetic one for the META cases")
        write(
            joinpath(REF[], "META"),
            "julia_version=$(VERSION)\nthreads=4\ndriver=driver.jl\ndriver_sha256=x\n"
        )
    end
    fresh(n) = (d=joinpath(tmp, n); cp(REF[], d); d)
    expect("identical copy", fresh("same"), true)

    # First outcome array of a store
    first_outcome(s) = first(
        p for
        p in readdir(joinpath(s, "results"); join=true) if isfile(joinpath(p, ".zarray"))
    )

    # 1. one value in the first float outcome array
    d = fresh("value")
    z = zopen(first_outcome(store(d, "mixed")), "w")
    z[ntuple(_ -> 1, ndims(z))...] += 1
    expect("one outcome value perturbed", d, false)

    # 2. attr change (non-allow-listed)
    d = fresh("attr")
    f = joinpath(store(d, "mixed"), "inputs", ".zattrs")
    a = JSON.parsefile(f)
    a["RCP"] = "perturbed"
    write(f, JSON.json(a))
    expect("inputs attr perturbed", d, false)

    # 3. allow-listed attr change must still pass
    d = fresh("allowed")
    f2 = joinpath(store(d, "mixed"), "inputs", ".zattrs")
    a = JSON.parsefile(f2)
    a["invoke_time"] = "other"
    write(f2, JSON.json(a))
    expect("allow-listed attr (invoke_time) changed", d, true)

    # 4. model_spec.csv
    d = fresh("spec")
    f3 = joinpath(store(d, "mixed"), "model_spec", "model_spec.csv")
    write(f3, read(f3, String) * "\n")
    expect("model_spec.csv perturbed", d, false)

    # 5. --tol must still require exact inputs table and intervention logs
    for p in ("inputs", joinpath("logs", "shading_log"), joinpath("logs", "rankings"))
        d = fresh("tol_" * replace(p, "/" => "_", "\\" => "_"))
        arr = joinpath(store(d, "mixed"), p)
        isfile(joinpath(arr, ".zarray")) ||
            (push!(results, false); println("FAIL missing ", p); continue)
        z = zopen(arr, "w")
        isfloat = eltype(z) <: AbstractFloat
        isfloat || (println("skip non-float ", p); continue)
        idx = ntuple(_ -> 1, ndims(z))
        v = z[idx...]
        z[idx...] = isnan(v) ? zero(v) : nextfloat(v)   # far below any tolerance
        expect("$p perturbed by < tol under --tol", d, false, `--tol 1e-4`)
    end

    # 6. run metadata (META, written by run_regression.sh)
    function with_meta(name; kw...)
        d = fresh(name)
        meta = isfile(joinpath(d, "META")) ? readlines(joinpath(d, "META")) : String[]
        filter!(l -> !any(k -> startswith(l, string(k) * "="), keys(kw)), meta)
        append!(meta, ["$k=$v" for (k, v) in kw])
        write(joinpath(d, "META"), join(meta, "\n") * "\n")
        return d
    end
    d = with_meta("meta_version"; julia_version="0.0.1")
    expect("other Julia version, exact", d, false)
    expect("other Julia version, under --tol 1e-4", d, true, `--tol 1e-4`)
    d = with_meta("meta_threads"; threads="123")
    expect("other thread count (warning only)", d, true)
    d = with_meta("meta_driver"; driver_sha256="0")
    expect("other driver (warning only)", d, true)
end
all(results) || (println("SELFTEST FAILED"); exit(1))
println("SELFTEST PASSED")
