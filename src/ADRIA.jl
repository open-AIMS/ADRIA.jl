module ADRIA

# Utility packages
using
    PrecompileTools,
    RelocatableFolders,
    TOML,
    CpuId,
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
include("ecosystem/corals/Corals.jl")
include("ecosystem/corals/GrowthAcceleration.jl")
include("ecosystem/cyclones.jl")
include("ecosystem/connectivity.jl")

include("Domain.jl")
include("io/datacubes.jl")
include("io/inputs.jl")  # Need to define input types before MCDA to make types available
include("io/initial_coral_cover.jl")

include("decision/Decision.jl")
include("interventions/Interventions.jl")
include("interventions/seeding.jl")
include("interventions/fogging.jl")
include("interventions/moving_corals.jl")

include("io/ResultSet.jl")
include("spatial/spatial.jl")
include("io/result_io.jl")
include("io/rme_result_io.jl")
include("io/result_post_processing.jl")
include("io/sampling/sampling.jl")
include("metrics/metrics.jl")
include("metrics/performance.jl")

include("scenario.jl")
include("analysis/analysis.jl")

include("ExtInterface/ADRIA/Domain.jl")
include("ExtInterface/ReefMod/RMEDomain.jl")
include("ExtInterface/ReefMod/ReefModDomain.jl")

include("viz/viz.jl")

export
    run_scenario, coral_spec, bin_edges,
    create_coral_struct, create_coral_instance,
    create_growth_acceleration_instance, Intervention, SimConstants,
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
# metric helper methods
# export dims, ndims

# List out compatible domain datapackages
const COMPAT_DPKG = ["0.8.0"]
end
