# build/analysis/make.jl
#
# Build adria_analysis.so: sysimage for ADRIA + ADRIAviz + ADRIAanalysis.
# MLJ and SIRUS are intentionally excluded — they must live in a weak-dep
# extension (separate refactor) before being candidates for the sysimage.
#
# Run from repo root inside the Docker builder stage:
#   julia --project=build/analysis build/analysis/make.jl
#
# Relevant env vars:
#   SYSIMAGE_OUT      output path (default: build/analysis/adria_analysis.so)
#   JULIA_CPU_TARGET  passed verbatim to PackageCompiler

using PackageCompiler

const _EXT = Sys.iswindows() ? ".dll" : ".so"
const SYSIMAGE_OUT = get(ENV, "SYSIMAGE_OUT",
    joinpath(@__DIR__, "adria_analysis" * _EXT))

const CPU_TARGET = get(ENV, "JULIA_CPU_TARGET",
    "x86_64;haswell;skylake;skylake-avx512;tigerlake")

@info "Building adria_analysis.so" sysimage_path = SYSIMAGE_OUT cpu_target = CPU_TARGET

mkpath(dirname(SYSIMAGE_OUT))

create_sysimage(
    [
        # ── Core simulation (superset of adria_sim.so) ──────────────────────
        :ADRIA,
        :ArchGDAL, :YAXArrays, :Zarr, :DimensionalData,
        :CoralBlox, :JMcDM,
        :Distributions, :Interpolations, :BlackBoxOptim, :OnlineStats,
        :DataFrames, :CSV, :JSON,
        :GeoDataFrames, :GeoInterface, :GeoFormatTypes,
        :Graphs, :SimpleWeightedGraphs,
        # ── Visualisation ───────────────────────────────────────────────────
        :ADRIAviz, :CairoMakie,
        # ── Analysis (MLJ / SIRUS deliberately absent) ───────────────────────
        :ADRIAanalysis,
        :Bootstrap, :Clustering, :DataEnvelopmentAnalysis, :HypothesisTests,
        :Distances, :StatsBase, :Statistics
    ];
    sysimage_path=SYSIMAGE_OUT,
    precompile_execution_file=joinpath(@__DIR__, "precompile_analysis.jl"),
    cpu_target=CPU_TARGET,
    incremental=true,
    filter_stdlibs=false
)

let sz = round(stat(SYSIMAGE_OUT).size / 1024^2; digits=1)
    @info "adria_analysis.so ready" size_MB = sz path = SYSIMAGE_OUT
end
