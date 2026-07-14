# build/sim/make.jl
#
# Build adria_sim.so: sysimage for ADRIA (simulation only — no plotting).
#
# Run from repo root inside the Docker builder stage:
#   julia --project=build/sim build/sim/make.jl
#
# Relevant env vars:
#   SYSIMAGE_OUT      output path (default: build/sim/adria_sim.so)
#   JULIA_CPU_TARGET  passed verbatim to PackageCompiler

using PackageCompiler

const _EXT = Sys.iswindows() ? ".dll" : ".so"
const SYSIMAGE_OUT = abspath(
    get(ENV, "SYSIMAGE_OUT",
        joinpath(@__DIR__, "adria_sim" * _EXT))
)

const CPU_TARGET = get(ENV, "JULIA_CPU_TARGET",
    "generic;haswell;skylake;skylake-avx512;tigerlake")

@info "Building adria_sim.so" sysimage_path = SYSIMAGE_OUT cpu_target = CPU_TARGET

mkpath(dirname(SYSIMAGE_OUT))

create_sysimage(
    [
        # ── Core simulation ─────────────────────────────────────────────────
        :ADRIA,
        # Heavy deps whose specialisations are worth baking in
        :ArchGDAL, :YAXArrays, :Zarr, :DimensionalData,
        :CoralBlox, :JMcDM,
        :Distributions, :Interpolations, :OnlineStats,
        # Data
        :DataFrames, :CSV, :JSON,
        :GeoDataFrames, :GeoInterface, :GeoFormatTypes,
        :Graphs, :SimpleWeightedGraphs
    ];
    sysimage_path=SYSIMAGE_OUT,
    precompile_execution_file=joinpath(@__DIR__, "precompile_sim.jl"),
    cpu_target=CPU_TARGET,
    incremental=true,   # build on top of Julia's default sys.so
    filter_stdlibs=false  # keep all stdlibs for compatibility
)

let sz = round(stat(SYSIMAGE_OUT).size / 1024^2; digits=1)
    @info "adria_sim.so ready" size_MB = sz path = SYSIMAGE_OUT
end

# Record the exact Julia build used to compile this sysimage, so consumers
# (see Dockerfile adria-sim-runtime) can verify the runtime `julia` matches
# before loading it with `-J`, rather than failing at container start.
write(SYSIMAGE_OUT * ".version", string(VERSION))
