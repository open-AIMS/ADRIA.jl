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
    FileIO,
    GeoInterface,
    Graphs,
    ImageIO,
    NetCDF,
    SimpleWeightedGraphs,
    Zarr

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

include("ecosystem/corals/coral_factors.jl")
include("ecosystem/corals/growth.jl")
include("ecosystem/corals/CoralGrowth.jl")
include("ecosystem/Ecosystem.jl")
include("ecosystem/DepthAttenuation.jl")
include("ecosystem/corals/Corals.jl")
include("ecosystem/corals/GrowthAcceleration.jl")
include("ecosystem/corals/flatten_overrides.jl")
include("ecosystem/cyclones.jl")
include("ecosystem/connectivity.jl")

include("Domain.jl")
include("io/datacubes.jl")
include("io/inputs.jl")  # Need to define input types before MCDA to make types available
include("io/initial_coral_cover.jl")
include("io/calib_params.jl")

include("decision/Decision.jl")
include("interventions/Interventions.jl")
include("interventions/seeding.jl")
include("interventions/fogging.jl")
include("interventions/moving_corals.jl")

include("io/ResultSet.jl")
include("spatial/spatial.jl")

# External model interfaces: domains and result sets.
# Must precede metrics/metrics.jl, which includes model-specific metric
# implementations (e.g. metrics/cscape.jl) that dispatch on these types.
include("ExtInterface/ADRIA/Domain.jl")
include("ExtInterface/ReefMod/RMEDomain.jl")
include("ExtInterface/ReefMod/ReefModDomain.jl")
include("ExtInterface/ReefMod/RMEResultSet.jl")
include("ExtInterface/CScape/CScapeResultSet.jl")

include("io/result_post_processing.jl")
include("io/result_io.jl")
include("io/sampling/sampling.jl")
include("metrics/metrics.jl")
include("metrics/performance.jl")

include("scenario.jl")
include("analysis/analysis.jl")

# Depends on metrics, so it is loaded after metrics/metrics.jl.
include("ExtInterface/ReefMod/RMEExport.jl")

include("viz/viz.jl")

export
    run_scenario, coral_spec, bin_edges,
    create_coral_struct, create_coral_instance,
    create_growth_acceleration_instance, Intervention, SimConstants, DepthAttenuation,
    SeedCriteriaWeights, FogCriteriaWeights, MCCriteriaWeights,
    loc_area, site_k_area, loc_k_area, loc_coral_cover, loc_recruits_cover,
    Domain, ADRIADomain,
    metrics, select, timesteps, env_stats, viz

using .analysis: AnnotatedOutcomes, attach_scenario_metadata
export AnnotatedOutcomes, attach_scenario_metadata

# Interfaces for external models
export RMEDomain
export ReefModDomain

export RMEResultSet
export CScapeResultSet
# metric helper methods
# export dims, ndims

# List out compatible domain datapackages
const COMPAT_DPKG = ["0.8.0"]

@compile_workload begin
    # DataCube / ZeroDataCube constructors — covers Dim{:timesteps/:scenarios/:locations} triples
    dc = DataCube(zeros(Float64, 5, 3, 10); timesteps=1:5, scenarios=1:3, locations=1:10)
    ZeroDataCube((:timesteps, :scenarios, :locations), (5, 3, 10))

    # Analysis helpers — pure computation, no filesystem I/O
    analysis.series_confint(randn(Float64, 10, 4))
    analysis.normalize(randn(Float64, 20))
    analysis.discretize_outcomes(randn(Float64, 20))

    # Metric dispatch paths — cover common array-input routing
    # scenario_total_cover: in_dims = (:timesteps, :locations, :scenarios)
    metrics.scenario_total_cover(rand(Float64, 5, 10, 3))
    # relative_cover: 5D plain-array path (timesteps×groups×sizes×locations×scenarios)
    metrics.relative_cover(rand(Float64, 5, 3, 10, 8, 2))
    # coral_evenness: in_dims = (:timesteps, :groups, :locations)
    metrics.coral_evenness(rand(Float64, 5, 3, 10))

    precompile(CSV.read, (String, Type{DataFrame}))
    # kwargs are dispatched through a generated "kwfunc" — precompile that too:
    precompile(
        Core.kwcall,
        (NamedTuple{(:comment,),Tuple{String}}, typeof(CSV.read), String, Type{DataFrame})
    )

    precompile(_build_dist_types, (DataFrame,))
    precompile(ADRIA.sample, (ADRIADomain, Int))
    precompile(ADRIA.sample_cf, (ADRIADomain, Int))
    precompile(ADRIA.sample_guided, (ADRIADomain, Int))
    precompile(ADRIA.sample_unguided, (ADRIADomain, Int))
    precompile(ADRIA.model_spec, (ADRIADomain, DataFrame))
    # `fix_factor!` writes the mutated spec back through `update!`; the workload block
    # below exercises the underlying reconstruct path against the assembled model.
    precompile(ADRIA.update!, (ADRIADomain, DataFrame))
    precompile(ADRIA.fix_factor!, (ADRIADomain, Symbol))
    precompile(ADRIA.fix_factor!, (ADRIADomain, Vector{Symbol}))

    redirect_stdio(stdout=devnull, stderr=devnull) do
        let
            n_locs, n_crit = 10, 4
            loc_names = Symbol.(:loc_, 1:n_locs)
            crit_names = [:in_conn, :out_conn, :heat_stress, :coral_cover]

            # Non-constant data (all-constant triggers filter_criteria, a separate path)
            dm_data = reshape(range(0.1, 1.0; length=n_locs * n_crit), n_locs, n_crit)

            # Mixed directions: MooraMethod requires at least one min and one max
            weights = fill(0.25, n_crit)
            directions = [maximum, maximum, minimum, maximum]

            dp = decision.DecisionPreferences(crit_names, weights, directions)
            dm = decision.decision_matrix(loc_names, crit_names, Matrix(dm_data))

            for method in decision.mcda_methods()
                try
                    decision.criteria_aggregated_scores(dp, dm, method)
                catch
                end
            end
        end
    end

    # Sampling dispatch — synthetic spec covers QMC + Distributions path without file IO
    # let
    #     spec = DataFrame(
    #         fieldname=[:p1, :p2],
    #         val=[0.5, 0.5],
    #         is_constant=[false, false],
    #         dist=[Uniform, Uniform],
    #         dist_params=[(0.0, 1.0), (0.0, 1.0)]
    #     )
    #     sample(spec, 4)
    # end

    # MCDA dispatch chain — 5 production methods, mixed directions for MooraMethod
    redirect_stdio(stdout=devnull, stderr=devnull) do
        let
            n_locs, n_crit = 10, 4
            loc_names = Symbol.(:loc_, 1:n_locs)
            crit_names = [:in_conn, :out_conn, :heat_stress, :coral_cover]

            # Non-constant data (all-constant triggers filter_criteria, a separate path)
            dm_data = reshape(range(0.1, 1.0; length=n_locs * n_crit), n_locs, n_crit)

            # Mixed directions: MooraMethod requires at least one min and one max
            weights = fill(0.25, n_crit)
            directions = [maximum, maximum, minimum, maximum]

            dp = decision.DecisionPreferences(crit_names, weights, directions)
            dm = decision.decision_matrix(loc_names, crit_names, Matrix(dm_data))

            for method in decision.mcda_methods()
                try
                    decision.criteria_aggregated_scores(dp, dm, method)
                catch
                end
            end
        end
    end

    # Full domain model assembly — covers the 9-component tuple form used by
    # load_domain() via _assemble_domain_model(). Without this workload entry,
    # Flatten.reconstruct for Tuple{EnvironmentalLayer,...,Coral,...,GrowthAcceleration}
    # is compiled on first interactive call (~20 s JIT). With this entry the cost
    # shifts to Pkg.precompile() time and interactive first-call is <1 s.
    redirect_stdio(stdout=devnull, stderr=devnull) do
        _arr = ones(Float32, 1, 1, 1)
        _model = _assemble_domain_model(
            EnvironmentalLayer(_arr, _arr, _arr),
            Intervention(),
            SeedCriteriaWeights(),
            FogCriteriaWeights(),
            MCCriteriaWeights(),
            DepthThresholds(),
            Coral(),
            GrowthAcceleration(),
            DepthAttenuation()
        )

        # Factor-fixing write-back path (see `update!` in Domain.jl for why the exact
        # `Float64` / `Vector{Float64}` column forms matter). Run twice: a freshly
        # assembled model has a `Float64`/`Int` `val` mix (`guided` etc. default to an
        # `Int` literal); the first write flips those to `Float64`, changing the `Param`
        # types handed to `Flatten.reconstruct`, and every `fix_factor!` after the first
        # sees that all-`Float64` shape. Values are written back unchanged -- compile
        # coverage only. Guarded so a failure degrades coverage rather than breaking
        # `using ADRIA`.
        try
            for _ in 1:2
                _ms = model_spec(_model)
                _model[:val] = Float64.(_ms.val)
                _model[:dist_params] = Vector{Float64}[
                    collect(Float64, _dp) for _dp in _ms.dist_params
                ]
            end
        catch err
            @debug "compile-workload factor write-back failed" err
        end
    end
end

# @compile_workload begin
#     if !isinteractive()
#         redirect_stdio(stdout=devnull, stderr=devnull) do
#             d = load_domain(joinpath(pkgdir(ADRIA), "test", "data", "Test_domain"), "45")
#             sample(d, 16)
#         end
#     end

#     # Cover the MCDA dispatch chain
#     # let
#     #     n_locs, n_crit = 10, 4
#     #     loc_names = Symbol.(:loc_, 1:n_locs)
#     #     crit_names = [:in_conn, :out_conn, :heat_stress, :coral_cover]

#     #     # Non-constant data (all-constant triggers filter_criteria, a separate path)
#     #     dm_data = reshape(range(0.1, 1.0; length=n_locs * n_crit), n_locs, n_crit)

#     #     # Mixed directions: Moora requires at least one min and one max
#     #     weights = fill(0.25, n_crit)
#     #     directions = [maximum, maximum, minimum, maximum]

#     #     dp = decision.DecisionPreferences(crit_names, weights, directions)
#     #     dm = decision.decision_matrix(loc_names, crit_names, Matrix(dm_data))

#     #     for method in decision.mcda_methods()
#     #         try
#     #             decision.criteria_aggregated_scores(dp, dm, method)
#     #         catch
#     #         end
#     #     end
#     # end
# end

end
