# Compare two result directories written by driver.jl by decoded contents.
#
#   julia --project=<checkout>/ADRIA compare.jl <ref_dir> <cand_dir> [--tol 1e-4]
#
# Exit code 0 = identical (or within --tol for outcomes), 1 = differences, 2 = usage error
# or an environment mismatch that makes an exact comparison meaningless (see META below).
# Uses only Zarr, JSON (direct ADRIA dependencies), so any ADRIA project environment works.
#
# Compared, per case directory and per Zarr array: relative path set, dtype, shape, chunks,
# fill_value and every value (NaN-aware exact equality, `isequal`), and attrs against an
# allow-list. `model_spec/model_spec.csv` is compared as text. The geopackage is ignored.
# Ranks, intervention logs (rankings, coral_aquaculture, larval_methods, shading_log) and the inputs table are ALWAYS exact; `--tol`
# only relaxes the other arrays (outcomes, other logs) to a max abs difference.
#
# META: run_regression.sh writes `<dir>/META` (key=value lines: julia_version, threads,
# driver, driver_sha256, source) next to the case dirs. A different Julia version fails an
# exact comparison (exit 2: results can differ in the last digits across versions); under
# `--tol` it is a warning. A different thread count or driver is a warning (a per-stage
# driver variant is expected to differ). A missing META is a warning.
using Zarr
using JSON

# Attr keys ignored entirely (run-specific); `*_file` attrs compare by basename
const IGNORED_ATTRS = Set(["invoke_time", "ADRIA_VERSION"])
# Always-exact arrays, by store-relative path (rankings, intervention logs, inputs table)
const ALWAYS_EXACT = [
    "inputs", "logs/rankings", "logs/coral_aquaculture", "logs/larval_methods",
    "logs/shading_log"
]
const IGNORED_DIRS = Set(["spatial", "spatial_data"])

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

function compare_array!(rep, id, path, ref_dir, cand_dir, tol)
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
    exact = any(p -> path == p || startswith(path, p * "/"), ALWAYS_EXACT)
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
        compare_array!(rep, "$name/$p", p, joinpath(ref, p), joinpath(cand, p), tol)
    end
    # Group-level attrs (store root, logs/, ...)
    groups(root) = sort!(
        [
        relpath(dir, root) for (dir, dirs, files) in walkdir(root) if
        ".zgroup" in files && !any(d -> startswith(relpath(dir, root), d), IGNORED_DIRS)
    ]
    )
    ga, gb = groups(ref), groups(cand)
    ga == gb || push!(rep.diffs, "$name: group sets differ: $ga vs $gb")
    for g in intersect(ga, gb)
        compare_attrs!(rep, "$name/$g (group)", joinpath(ref, g), joinpath(cand, g))
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

function read_meta(dir)::Union{Nothing,Dict{String,String}}
    f = joinpath(dir, "META")
    isfile(f) || return nothing
    kv = [split(l, "="; limit=2) for l in readlines(f) if occursin("=", l)]
    return Dict(String(k) => String(strip(v)) for (k, v) in kv)
end

# Print the run metadata of both sides; return false when the comparison must not go ahead
function check_meta(ref, cand, tol)::Bool
    ma, mb = read_meta(ref), read_meta(cand)
    for (side, m) in (("reference", ma), ("candidate", mb))
        if isnothing(m)
            println(
                "WARNING: no META in the $side ($side built before META was recorded?): ",
                "Julia version, threads and driver cannot be checked")
        else
            println("$side: Julia ", get(m, "julia_version", "?"), ", ",
                get(m, "threads", "?"), " threads, ", get(m, "driver", "?"), " (sha256 ",
                first(get(m, "driver_sha256", "?"), 12), "), ", get(m, "source", "?"))
        end
    end
    (isnothing(ma) || isnothing(mb)) && return true
    get_(m, k) = get(m, k, "?")
    if get_(ma, "julia_version") != get_(mb, "julia_version")
        msg =
            "reference built with Julia $(get_(ma, "julia_version")), candidate run " *
            "with Julia $(get_(mb, "julia_version")): results can differ in the last " *
            "digits across Julia versions"
        if isnothing(tol)
            println(stderr, "ERROR: $msg, so an exact comparison is not meaningful. ",
                "Rebuild the reference with the same Julia (JULIA=...), or pass --tol X ",
                "to compare across versions.")
            return false
        end
        println("WARNING: $msg; comparing under --tol $tol")
    end
    get_(ma, "threads") == get_(mb, "threads") || println(
        "WARNING: thread counts differ ($(get_(ma, "threads")) vs $(get_(mb, "threads")));",
        " use the same ADRIA_REGRESSION_THREADS for reference and candidate"
    )
    get_(ma, "driver_sha256") == get_(mb, "driver_sha256") || println(
        "WARNING: drivers differ ($(get_(ma, "driver")) vs $(get_(mb, "driver"))): ",
        "expected only for a per-stage driver variant with the same ensemble"
    )
    return true
end

function main(args)
    tol = nothing
    pos = String[]
    i = 1
    while i <= length(args)
        if args[i] == "--tol"
            i < length(args) || (println(stderr, "--tol needs a value"); return 2)
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
    check_meta(ref, cand, tol) || return 2
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
