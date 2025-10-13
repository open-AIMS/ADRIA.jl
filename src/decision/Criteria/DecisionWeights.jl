using ModelParameters
using Distributions, Statistics

abstract type DecisionWeights end
abstract type DecisionThresholds end

"""
    DepthThresholds <: DecisionThresholds

Properties used to filter out locations.
"""
Base.@kwdef struct DepthThresholds <: DecisionThresholds
    depth_min::Param = Factor(
        5.0;
        ptype="ordered categorical",
        dist=DiscreteOrderedUniformDist,
        dist_params=(2.0, 5.0, 0.5),  # Half metre intervals
        name="Minimum Depth",
        description="Minimum depth for a site to be included for consideration.\nNote: This value will be replaced with the shallowest depth value found if all sites are found to be deeper than `depth_min + depth_offset`."
    )
    depth_offset::Param = Factor(
        10.0;
        ptype="ordered categorical",
        dist=DiscreteOrderedUniformDist,
        dist_params=(10.0, 25.0, 0.5),  # Half metre intervals
        name="Depth Offset",
        description="Offset from minimum depth, used to indicate maximum depth."
    )
end

include("SeedCriteria.jl")
include("FogCriteria.jl")
