# Compare two result directories written by driver.jl by decoded contents.
#
#   julia --project=<checkout>/ADRIA compare.jl <ref_dir> <cand_dir> [--tol 1e-4]
#
# Exit code 0 = identical (or within --tol for outcomes), 1 = differences, 2 = usage error.
# Uses only Zarr, JSON (direct ADRIA dependencies), so any ADRIA project environment works.
#
# Compared, per case directory and per Zarr array: relative path set, dtype, shape, chunks,
# fill_value and every value (NaN-aware exact equality, `isequal`), and attrs against an
# allow-list. `model_spec/model_spec.csv` is compared as text. The geopackage is ignored.
# Ranks, intervention logs (seed/fog/shade) and the inputs table are ALWAYS exact; `--tol`
# only relaxes the other arrays (outcomes, other logs) to a max abs difference.
using Zarr
using JSON

# Attr keys ignored entirely (run-specific); `*_file` attrs compare by basename
const IGNORED_ATTRS = Set(["invoke_time", "ADRIA_VERSION", "name"])
# Always-exact arrays, by path prefix inside a store
const ALWAYS_EXACT = ["inputs", "ranks", "logs/seed", "logs/fog", "logs/shade"]
const IGNORED_DIRS = Set(["spatial_data"])

mutable struct Report
    diffs::Vector{String}
    n_arrays::Int
end

function find_arrays(root::String)
    found = String[]
    for (dir, dirs, files) in walkdir(root)
        filter!(d -> d ∉ IGNORED_DIRS, dirs)
        ".zarray" in files && push!(found, relpath(dir, root))
    end
    return sort!(found)
end

function normalize_attr(key::String, val)
    endswith(key, "_file") && val isa AbstractString && return basename(val)
    return val
end

function compare_attrs!(rep, id, ref_dir, cand_dir)
    read_attrs(d) =
        isfile(joinpath(d, ".zattrs")) ? JSON.parsefile(joinpath(d, ".zattrs")) : Dict()
    a, b = read_attrs(ref_dir), read_attrs(cand_dir)
    for k in sort!(collect(union(keys(a), keys(b))))
        k in IGNORED_ATTRS && continue
        if !haskey(a, k) || !haskey(b, k)
            push!(rep.diffs, "$id: attr '$k' present in only one side")
        elseif !isequal(normalize_attr(k, a[k]), normalize_attr(k, b[k]))
            push!(rep.diffs, "$id: attr '$k' differs")
        end
    end
end

function compare_array!(rep, id, ref_dir, cand_dir, tol)
    rep.n_arrays += 1
    ma = JSON.parsefile(joinpath(ref_dir, ".zarray"))
    mb = JSON.parsefile(joinpath(cand_dir, ".zarray"))
    for f in ("dtype", "shape", "chunks", "fill_value", "order")
        isequal(ma[f], mb[f]) ||
            push!(rep.diffs, "$id: $f differs ($(ma[f]) vs $(mb[f]))")
    end
    compare_attrs!(rep, id, ref_dir, cand_dir)
    (isequal(ma["shape"], mb["shape"]) && isequal(ma["dtype"], mb["dtype"])) ||
        return nothing

    a = Array(zopen(ref_dir))
    b = Array(zopen(cand_dir))
    exact = any(p -> id == p || startswith(id, p * "/"), ALWAYS_EXACT)
    if exact || isnothing(tol) || !(eltype(a) <: AbstractFloat)
        if !isequal(a, b)
            n = count(.!isequal.(a, b))
            push!(rep.diffs, "$id: $n of $(length(a)) values differ")
        end
    else
        nan_a, nan_b = isnan.(a), isnan.(b)
        if nan_a != nan_b
            push!(rep.diffs, "$id: NaN positions differ")
        else
            keep = .!nan_a
            d = maximum(abs.(Float64.(a[keep]) .- Float64.(b[keep])); init=0.0)
            d <= tol || push!(rep.diffs, "$id: max abs diff $d exceeds tol $tol")
        end
    end
    return nothing
end

# Lines of model_spec.csv without the `# Generated with ... on <time>` header comment
spec_lines(f) = filter(l -> !startswith(l, "#"), readlines(f))

function compare_store!(rep, name, ref, cand, tol)
    ra, rb = find_arrays(ref), find_arrays(cand)
    for p in setdiff(ra, rb)
        push!(rep.diffs, "$name/$p: missing in candidate")
    end
    for p in setdiff(rb, ra)
        push!(rep.diffs, "$name/$p: extra in candidate")
    end
    for p in intersect(ra, rb)
        compare_array!(rep, "$name/$p", joinpath(ref, p), joinpath(cand, p), tol)
    end
    ms = joinpath("model_spec", "model_spec.csv")
    fa, fb = joinpath(ref, ms), joinpath(cand, ms)
    if isfile(fa) != isfile(fb)
        push!(rep.diffs, "$name/$ms: present in only one side")
    elseif isfile(fa) && spec_lines(fa) != spec_lines(fb)
        push!(rep.diffs, "$name/$ms: differs")
    end
end

# Each case dir holds exactly one store directory (named with the invoke time)
function the_store(case_dir)
    stores = filter(isdir, readdir(case_dir; join=true))
    length(stores) == 1 || error("expected one store in $case_dir, found $(length(stores))")
    return only(stores)
end

function main(args)
    tol = nothing
    pos = String[]
    i = 1
    while i <= length(args)
        if args[i] == "--tol"
            tol = parse(Float64, args[i + 1])
            i += 2
        else
            push!(pos, args[i])
            i += 1
        end
    end
    length(pos) == 2 ||
        (println(stderr, "usage: compare.jl ref_dir cand_dir [--tol X]"); return 2)
    ref, cand = pos
    rep = Report(String[], 0)
    ca = sort!(filter(d -> isdir(joinpath(ref, d)), readdir(ref)))
    cb = sort!(filter(d -> isdir(joinpath(cand, d)), readdir(cand)))
    ca == cb || push!(rep.diffs, "case sets differ: $ca vs $cb")
    for c in intersect(ca, cb)
        compare_store!(
            rep, c, the_store(joinpath(ref, c)), the_store(joinpath(cand, c)), tol
        )
    end
    mode = isnothing(tol) ? "exact" : "tol=$tol (outcomes only)"
    if isempty(rep.diffs)
        println("PASS ($mode): $(rep.n_arrays) arrays compared")
        return 0
    end
    println("FAIL ($mode): $(length(rep.diffs)) difference(s)")
    foreach(d -> println("  ", d), first(rep.diffs, 50))
    return 1
end

exit(main(ARGS))
