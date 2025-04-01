module ADRIA

# Utility packages
using
    PrecompileTools,
    RelocatableFolders,
    TOML,
    CpuId,
    PkgVersion,
    ProgressMeter,
    Dates

# IO and functionality packages
import GeoDataFrames as GDF
import GeoFormatTypes as GFT
import ArchGDAL as AG
using
    CSV,
    DataFrames,
    Distributed,
    FLoops,
    FileIO,
    GeoInterface,
    Graphs,
    ImageIO,
    MAT,
    SimpleWeightedGraphs

# Modelling packages
using
    Combinatorics,
    DataStructures,
    DimensionalData,
    Distances,
    Distributions,
    CoralBlox,
    LinearAlgebra,
    OrderedCollections,
    Setfield,
    Statistics,
    StaticArrays,
    SparseArrays,
    Random,
    YAXArrays

include("utils/text_display.jl")  # need better name for this file
include("utils/setup.jl")
include("utils/scale.jl")
include("factors/Factors.jl")
include("factors/const_params.jl")

include("ecosystem/corals/growth.jl")
include("ecosystem/corals/CoralGrowth.jl")
include("ecosystem/Ecosystem.jl")
include("ecosystem/corals/Corals.jl")
include("ecosystem/connectivity.jl")

include("Domain.jl")
include("io/datacubes.jl")
include("io/inputs.jl")  # Need to define input types before MCDA to make types available
include("io/initial_coral_cover.jl")

include("decision/dMCDA.jl")
include("interventions/Interventions.jl")
include("interventions/seeding.jl")
include("interventions/fogging.jl")

include("io/ResultSet.jl")
include("spatial/spatial.jl")
include("io/result_io.jl")
include("io/rme_result_io.jl")
include("io/result_post_processing.jl")
include("io/sampling.jl")
include("metrics/metrics.jl")
include("metrics/performance.jl")

include("scenario.jl")
include("analysis/analysis.jl")
include("analysis/sensitivity.jl")

include("ExtInterface/ADRIA/Domain.jl")
include("ExtInterface/ReefMod/RMEDomain.jl")
include("ExtInterface/ReefMod/ReefModDomain.jl")

include("viz/viz.jl")

export
    growthODE,
    run_scenario, coral_spec,
    create_coral_struct, Intervention, SimConstants,
    SeedCriteriaWeights, FogCriteriaWeights,
    loc_area, site_k_area, loc_k_area,
    Domain, ADRIADomain,
    metrics, select, timesteps, env_stats, viz

# Interfaces for external models
export RMEDomain
export ReefModDomain

export RMEResultSet
# metric helper methods
# export dims, ndims

# List out compatible domain datapackages
const COMPAT_DPKG = ["0.7.0-rc", "0.7.0", "0.1.0-gbr"]
# This adds ~30 seconds to package load times
if ccall(:jl_generating_output, Cint, ()) == 1
    Base.precompile(Tuple{typeof(load_domain),String})   # time: 19.120537
    Base.precompile(Tuple{typeof(load_domain),String,String})
    Base.precompile(Tuple{typeof(setup_result_store!),Domain,DataFrame})   # time: 4.6720815
    Base.precompile(Tuple{typeof(combine_results),Vector{String}})   # time: 4.0178256
    Base.precompile(
        Tuple{
            typeof(growthODE),
            Matrix{Float64},
            Matrix{Float64},
            NamedTuple{
                (
                    :r,
                    :k,
                    :mb,
                    :comp,
                    :sm_comp,
                    :small_massives,
                    :small,
                    :mid,
                    :large,
                    :acr_5_11,
                    :acr_6_12,
                    :rec,
                    :sigma,
                    :M_sm,
                    :sXr,
                    :X_mb,
                    :cover
                ),
                Tuple{
                    Matrix{Float64},
                    Vector{Float64},
                    Matrix{Float64},
                    Float64,
                    Matrix{Float64},
                    SVector{3,Int64},
                    SVector{6,Int64},
                    SVector{19,Int64},
                    SVector{4,Int64},
                    SVector{2,Int64},
                    SVector{2,Int64},
                    Matrix{Float64},
                    Matrix{Float64},
                    Matrix{Float64},
                    Matrix{Float64},
                    Matrix{Float64},
                    Vector{Float64}
                }
            },
            Float64
        }
    )   # time: 1.4354926
    Base.precompile(
        Tuple{
            typeof(decision.rank_sites!),
            Matrix{Float64},
            Vector{Float64},
            Matrix{Int64},
            Int64,
            typeof(decision.adria_topsis),
            Int64
        }
    )   # time: 0.3518593
    Base.precompile(
        Tuple{
            typeof(decision.rank_sites!),
            Matrix{Float64},
            Vector{Float64},
            Matrix{Int64},
            Int64,
            typeof(decision.adria_vikor),
            Int64
        }
    )   # time: 0.3170264
    Base.precompile(
        Tuple{
            typeof(scenario_attributes),
            String,
            String,
            Vector{String},
            String,
            EnvLayer{String,Vector{Int64}},
            SimConstants,
            Vector{String},
            Vector{Float64},
            Vector{Float64},
            Vector{Tuple{Float64,Float64}}
        }
    )   # time: 0.2140636
    Base.precompile(Tuple{typeof(model_spec),Model})   # time: 0.1997914
    Base.precompile(
        Tuple{
            typeof(bleaching_mortality!),
            Matrix{Float64},
            Matrix{Float64},
            Vector{Float64},
            Int64,
            Vector{Float64},
            Vector{Float64},
            Vector{Float64},
            Vector{Float64},
            Float64
        }
    )   # time: 0.1940948
    Base.precompile(
        Tuple{
            typeof(decision.decision_matrix),
            Vector{Int64},
            Vector{Float64},
            Vector{Float64},
            Vector{Float64},
            Vector{Float64},
            Vector{Float64},
            Vector{Float64},
            Vector{Float64},
            Vector{Float64},
            Matrix{Float64},
            Matrix{Float64},
            Float64
        }
    )   # time: 0.1929096
    Base.precompile(
        Tuple{
            typeof(scenario_attributes),
            String,
            String,
            Vector{String},
            String,
            EnvLayer{String,Vector{Any}},
            Dict{String,Any},
            Vector{Any},
            Vector{Float64},
            Vector{Float64},
            Vector{Any}
        }
    )   # time: 0.1755622
    Base.precompile(Tuple{typeof(proportional_adjustment!),Matrix{Float64},Vector{Float64}})   # time: 0.1680073
    Base.precompile(Tuple{typeof(_remove_workers)})   # time: 0.1593244
    Base.precompile(Tuple{typeof(_setup_workers)})   # time: 0.1571776
    Base.precompile(Tuple{typeof(switch_RCPs!),Domain,String})   # time: 0.1284853
    # Base.precompile(Tuple{typeof(component_params),DataFrame,Type{CriteriaWeights}})   # time: 0.1223987
    Base.precompile(
        Tuple{
            Type{Domain},
            String,
            String,
            String,
            Vector{Int64},
            String,
            String,
            String,
            String,
            String,
            String,
            String
        }
    )   # time: 0.1113899
    Base.precompile(Tuple{typeof(setup_cache),Domain})   # time: 0.1060752
    Base.precompile(
        EnvLayer, (String, String, String, String, String, String, String, String, Any)
    )
    Base.precompile(load_results, (String,))
end

# @setup_workload begin
#     # Putting some things in `setup` can reduce the size of the
#     # precompile file and potentially make loading faster.
#     # ADRIA_DIR = pkgdir(ADRIA)
#     # EXAMPLE_DOMAIN_PATH = joinpath(ADRIA_DIR, "examples", "Example_domain")

#     @compile_workload begin

#         # Force precompiling of code that handles distributed infrastructure
#         addprocs(1)
#         @everywhere 1 + 1

#         # Compile progress bar
#         # f() = begin
#         #     @showprogress 1 for _ in 1:10
#         #     end
#         # end
#         # b = redirect_stdout(f, devnull)

#         # dom = ADRIA.load_domain(EXAMPLE_DOMAIN_PATH, "45")
#         # ADRIA.sample(dom, 16)
#         # ADRIA.model_spec(dom)

#         # ENV["ADRIA_DEBUG"] = "false"
#         # p_df = ADRIA.param_table(dom)
#         # rs1 = ADRIA.run_scenario(dom, p_df[1, :])
#         # delete!(ENV, "ADRIA_DEBUG")
#     end
# end

end
