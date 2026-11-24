# Run metadata (META) helpers shared by run_regression.jl and compare.jl.
# Standard library only, so the runner works under any Julia environment.

using SHA

"""Keys every META must hold, each with a non-empty value."""
const META_KEYS = (
    "julia_version", "threads", "driver", "driver_sha256", "source", "os", "arch"
)

"""
Keys whose value must match between the reference and the candidate for the comparison to
mean anything. A mismatch of `os` or `arch` is always a refusal, even under `--tol`.
"""
const ENV_KEYS = ("julia_version", "threads", "driver_sha256", "os", "arch")

"""
    driver_sha256(path::AbstractString)::String

The sha256 of the driver script with CRLF line ends normalised to LF, so that
`core.autocrlf` (Windows checkouts) does not change the recorded checksum.
"""
function driver_sha256(path::AbstractString)::String
    bytes = read(path)
    out = UInt8[]
    sizehint!(out, length(bytes))
    for (i, b) in enumerate(bytes)
        (b == 0x0d && i < length(bytes) && bytes[i + 1] == 0x0a) && continue
        push!(out, b)
    end
    return bytes2hex(sha256(out))
end

"""
    host_meta()::Dict{String,String}

The operating system and CPU architecture of this process.
"""
function host_meta()::Dict{String,String}
    return Dict("os" => string(Sys.KERNEL), "arch" => string(Sys.ARCH))
end

"""
    read_meta(dir)::Union{Nothing,Dict{String,String}}

The `key=value` lines of `<dir>/META` with keys and values whitespace-trimmed (also a CR),
or `nothing` if there is no META.
"""
function read_meta(dir)::Union{Nothing,Dict{String,String}}
    f = joinpath(dir, "META")
    isfile(f) || return nothing
    kv = [split(l, "="; limit=2) for l in readlines(f) if occursin("=", l)]
    return Dict(String(strip(k)) => String(strip(v)) for (k, v) in kv)
end

"""
    meta_lacking(m)::Vector{String}

The keys of `META_KEYS` that `m` lacks a non-empty value for.
"""
meta_lacking(m) = String[k for k in META_KEYS if isempty(get(m, k, ""))]

"""
    meta_differing(m, expected)::Vector{String}

The keys of `expected` whose value in `m` differs (a missing key counts as different).
"""
function meta_differing(m, expected)::Vector{String}
    return String[k for k in sort!(collect(keys(expected))) if get(m, k, "") != expected[k]]
end
