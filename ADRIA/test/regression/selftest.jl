# Self-test for compare.jl: perturbed copies of a reference must FAIL, the unperturbed copy
# must PASS.  julia --project=<ADRIA env> selftest.jl <ref_dir> <compare.jl>
using Zarr
using JSON

ref, cmp = ARGS
julia = Base.julia_cmd()
project = dirname(Base.active_project())

passes(cand) =
    run(
        pipeline(ignorestatus(`$julia --project=$project $cmp $ref $cand`);
            stdout=devnull)
    ).exitcode == 0

store(root, c) = only(filter(isdir, readdir(joinpath(root, c); join=true)))
results = Bool[]
function expect(name, cand, want)
    ok = passes(cand) == want
    println(ok ? "ok   " : "FAIL ", name, " (expected ", want ? "pass" : "fail", ")")
    push!(results, ok)
end

mktempdir() do tmp
    fresh(n) = (d=joinpath(tmp, n); cp(ref, d); d)
    expect("identical copy", fresh("same"), true)

    # 1. one value in the first float outcome array
    d = fresh("value")
    s = store(d, "mixed")
    arr = first(
        p for
        p in readdir(joinpath(s, "results"); join=true) if isfile(joinpath(p, ".zarray"))
    )
    z = zopen(arr, "w")
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
end
all(results) || (println("SELFTEST FAILED"); exit(1))
println("SELFTEST PASSED")
