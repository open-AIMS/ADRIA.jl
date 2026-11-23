# Self-test for compare.jl: perturbed copies of a reference must fail with exit code 1
# (differences), environment mismatches with 2 (refusal), and the unperturbed copy (and
# changes compare.jl is meant to accept) must pass (0). Exit codes are checked exactly.
#   julia --project=<ADRIA env> selftest.jl <ref_dir> <compare.jl>
using Zarr
using JSON

length(ARGS) == 2 || (println(stderr, "usage: selftest.jl ref_dir compare.jl"); exit(2))
src_ref, cmp = ARGS
isdir(src_ref) || (println(stderr, "ERROR: reference folder $src_ref not found"); exit(2))
isfile(joinpath(src_ref, "META")) || (
    println(stderr, "ERROR: reference $src_ref has no META (built before META was ",
        "recorded); rebuild it with run_regression.sh ref"); exit(2)
)
isfile(cmp) || (println(stderr, "ERROR: $cmp not found"); exit(2))
julia = Base.julia_cmd()
project = dirname(Base.active_project())
# The reference the candidates are compared against: a copy of `src_ref`
const REF = Ref("")

# Exit codes of compare.jl
const PASS, DIFF, REFUSED, CRASH = 0, 1, 2, 3
const CODE_NAMES = Dict(PASS => "pass", DIFF => "differences", REFUSED => "refusal",
    CRASH => "crash")

exitcode(cand, extra=``; ref=REF[]) =
    run(
        pipeline(
            ignorestatus(`$julia --project=$project $cmp $ref $cand $extra`);
            stdout=devnull, stderr=devnull
        )
    ).exitcode

store(root, c) = only(filter(isdir, readdir(joinpath(root, c); join=true)))
results = Bool[]
# The exact exit code is checked, so a crash (3) cannot pass as expected differences (1)
function expect(name, cand, want::Int, extra=``; ref=REF[])
    got = exitcode(cand, extra; ref=ref)
    ok = got == want
    println(ok ? "ok   " : "FAIL ", name, " (expected exit $want = ",
        CODE_NAMES[want], ok ? "" : ", got $got", ")")
    push!(results, ok)
end

mktempdir() do tmp
    REF[] = joinpath(tmp, "ref")
    cp(src_ref, REF[])
    fresh(n) = (d=joinpath(tmp, n); cp(REF[], d); d)
    expect("identical copy", fresh("same"), PASS)

    # First outcome array of a store
    first_outcome(s) = first(
        p for
        p in readdir(joinpath(s, "results"); join=true) if isfile(joinpath(p, ".zarray"))
    )

    # 1. one value in the first float outcome array
    d = fresh("value")
    z = zopen(first_outcome(store(d, "mixed")), "w")
    z[ntuple(_ -> 1, ndims(z))...] += 1
    expect("one outcome value perturbed", d, DIFF)

    # 1b. a perturbation within tolerance: fails exact, passes under --tol
    d = fresh("within_tol")
    z = zopen(first_outcome(store(d, "mixed")), "w")
    idx = ntuple(_ -> 1, ndims(z))
    z[idx...] = z[idx...] + eltype(z)(1e-6)
    expect("outcome perturbed by 1e-6, exact", d, DIFF)
    expect("outcome perturbed by 1e-6, under --tol 1e-4", d, PASS, `--tol 1e-4`)

    # 1c. an array missing from, or extra in, the candidate
    d = fresh("missing")
    rm(first_outcome(store(d, "mixed")); recursive=true)
    expect("outcome array missing", d, DIFF)
    d = fresh("extra")
    src = first_outcome(store(d, "mixed"))
    cp(src, joinpath(dirname(src), "extra_" * basename(src)))
    expect("extra outcome array", d, DIFF)

    # 2. attr change (non-allow-listed)
    d = fresh("attr")
    f = joinpath(store(d, "mixed"), "inputs", ".zattrs")
    a = JSON.parsefile(f)
    a["RCP"] = "perturbed"
    write(f, JSON.json(a))
    expect("inputs attr perturbed", d, DIFF)

    # 3. allow-listed attr change must still pass
    d = fresh("allowed")
    f2 = joinpath(store(d, "mixed"), "inputs", ".zattrs")
    a = JSON.parsefile(f2)
    a["invoke_time"] = "other"
    write(f2, JSON.json(a))
    expect("allow-listed attr (invoke_time) changed", d, PASS)

    # 4. model_spec.csv
    d = fresh("spec")
    f3 = joinpath(store(d, "mixed"), "model_spec", "model_spec.csv")
    write(f3, read(f3, String) * "\n")
    expect("model_spec.csv perturbed", d, DIFF)

    # 4b. the attached-submodel case is compared too
    d = fresh("attached_value")
    z = zopen(first_outcome(store(d, "attached")), "w")
    z[ntuple(_ -> 1, ndims(z))...] += 1
    expect("attached case outcome value perturbed", d, DIFF)

    # 4c. case sets differ: a case missing from, or extra in, the candidate
    d = fresh("case_missing")
    rm(joinpath(d, "cots_off"); recursive=true)
    expect("case missing in candidate", d, DIFF)
    d = fresh("case_extra")
    cp(joinpath(d, "cots_off"), joinpath(d, "extra_case"))
    expect("extra case in candidate", d, DIFF)

    # 4d. a layout compare.jl cannot handle (two stores in a case) is a crash, not a diff
    d = fresh("two_stores")
    s1 = store(d, "cots_off")
    cp(s1, s1 * "_copy")
    expect("two stores in a case (compare.jl error)", d, CRASH)

    # 5. --tol must still require exact inputs table and intervention logs
    for p in ("inputs", joinpath("logs", "shading_log"), joinpath("logs", "rankings"))
        d = fresh("tol_" * replace(p, "/" => "_", "\\" => "_"))
        arr = joinpath(store(d, "mixed"), p)
        isfile(joinpath(arr, ".zarray")) ||
            (push!(results, false); println("FAIL missing ", p); continue)
        z = zopen(arr, "w")
        T = eltype(z)
        idx = ntuple(_ -> 1, ndims(z))
        v = z[idx...]
        if T <: AbstractFloat
            z[idx...] = isnan(v) ? zero(v) : nextfloat(v)   # far below any tolerance
        elseif T <: Integer
            z[idx...] = v == typemax(T) ? v - one(T) : v + one(T)  # e.g. ranks (UInt16)
        else
            push!(results, false)
            println("FAIL cannot perturb $p of eltype $T")
            continue
        end
        expect("$p ($T) perturbed under --tol", d, DIFF, `--tol 1e-4`)
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
    expect("other Julia version, exact", d, REFUSED)
    expect("other Julia version, under --tol 1e-4", d, PASS, `--tol 1e-4`)
    d = with_meta("meta_threads"; threads="123")
    expect("other thread count (warning only)", d, PASS)
    d = with_meta("meta_driver"; driver_sha256="0")
    expect("other driver (warning only)", d, PASS)
    d = fresh("meta_missing")
    rm(joinpath(d, "META"))
    expect("candidate without META, exact", d, REFUSED)
    expect("candidate without META, under --tol 1e-4", d, PASS, `--tol 1e-4`)
    expect("reference without META, exact", REF[], REFUSED; ref=d)
    # Empty values on both sides must not count as a match
    d = with_meta("meta_empty"; julia_version="")
    expect("empty julia_version on both sides, exact", d, REFUSED; ref=d)
    d = with_meta("meta_no_source"; source="")
    expect("empty source in candidate, exact", d, REFUSED)
end
all(results) || (println("SELFTEST FAILED"); exit(1))
println("SELFTEST PASSED")
