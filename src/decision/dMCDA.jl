module decision

using InteractiveUtils: subtypes
using StatsBase

using
    Combinatorics,
    DataFrames,
    JMcDM

using ADRIA: Factor, DiscreteOrderedUniformDist, component_params


# dummy dMCDA_vars() for dev
struct DMCDA_vars end

# dummy functions to allow precompilation
function guided_site_selection() end
function unguided_site_selection() end
function rank_sites!() end
function adria_topsis() end
function adria_vikor() end
function create_decision_matrix() end

function supported_jmcdm_methods()
    return [
        JMcDM.Topsis.TopsisMethod,
        JMcDM.VIKOR.VikorMethod,
        JMcDM.ARAS.ArasMethod,
        JMcDM.COCOSO.CocosoMethod,
        JMcDM.CODAS.CodasMethod,
        JMcDM.EDAS.EdasMethod,
        JMcDM.GREY.GreyMethod,
        JMcDM.MABAC.MabacMethod,
        JMcDM.MAIRCA.MaircaMethod,
        JMcDM.MARCOS.MarcosMethod,
        JMcDM.MOORA.MooraMethod,
        JMcDM.PIV.PIVMethod,
        JMcDM.PSI.PSIMethod,
        JMcDM.SAW.SawMethod,
        JMcDM.WASPAS.WaspasMethod,
        JMcDM.WPM.WPMMethod,
    ]
end

function mcda_methods()
    return supported_jmcdm_methods()
end

"""
    mcda_normalize(x::Vector)::Vector
    mcda_normalize(x::Matrix)::Matrix
    mcda_normalize(x::DataFrame)::DataFrame

Normalize a vector (wse/wsh), decision matrix, or weights for a scenario set for the
purpose of MCDA.
"""
function mcda_normalize(x::Vector)::Vector
    return x ./ sum(x)
end
function mcda_normalize(x::Matrix)::Matrix
    return x ./ sqrt.(sum(x .^ 2; dims = 1))
end
function mcda_normalize(x::DataFrame)::DataFrame
    return x ./ sum(Matrix(x); dims = 2)
end

"""
    align_rankings!(rankings::Array, s_order::Matrix, col::Int64)::Nothing

Align a vector of site rankings to match the indicated order in `s_order`.
"""
function align_rankings!(rankings::Array, s_order::Matrix, col::Int64)::Nothing
    # Fill target ranking column
    for (i, site_id) in enumerate(s_order[:, 1])
        rankings[rankings[:, 1] .== site_id, col] .= i
    end

    return nothing
end

"""
    within_depth_bounds(loc_depth::Vector{T}, depth_max::T, depth_min::T)::BitVector{T} where {T<:Float64}

Determines whether a location is within the min/max depth bounds.
Used to filter locations based on their depth for location selection.

# Arguments
- `loc_depth` : Depths of considered locations (typically the median depth)
- `depth_max` : Maximum depth for each considered location
- `depth_min` : Minimum depth for each considered location

# Returns
BitVector, of logical indices indicating locations which satisfy the depth criteria.
"""
function within_depth_bounds(
    loc_depth::Vector{T}, depth_max::T, depth_min::T
)::BitVector where {T <: Float64}
    return (loc_depth .<= depth_max) .& (loc_depth .>= depth_min)
end

"""
    summary_stat_env(env_layer::NamedDimsArray, dims::Union{Int64, Symbol, Tuple{Symbol, Symbol}}; w=0.5)::Vector{Float64}

Calculates weighted combinations of mean and standard deviation for a given environmental
factor.

# Arguments
- `env_layer` : Environmental data layer to calculate the mean of.
- `dims` : Dimensions to aggregate over.
- `w` : Weighting for std offset to mean.

# Returns
Weighted combination of mean and standard deviation of the projected environmental
conditions (e.g., DHWs, wave stress, etc):
    (μ * w) + (σ * (1 - w))
"""
function summary_stat_env(
    env_layer::AbstractArray,
    dims::Union{Int64, Symbol, Tuple{Symbol, Symbol}};
    w = 0.5,
)::Vector{Float64}
    return vec((mean(env_layer; dims = dims) .* w) .+ (std(env_layer; dims = dims) .* (1.0 - w)))
end


"""
    unguided_selection(
        n_iv_locs::Int64,
        k_area::Vector{Float64},
        depth::Vector{Int64}
    )

Randomly select intervention locations, constraining to locations that are able to support
corals and are within the desired depth range.

# Arguments
- `n_iv_locs` : Number of locations to seed
- `k_area` : Coral habitable available at each location (`k` value) in either relative or
             absolute units.
- `depth` : vector of location ids found to be within desired depth range

# Returns
Matrix, Name/IDs of selected locations, and their indices
"""
function unguided_selection(
    location_ids,
    n_iv_locs::Int64,
    k_area::Vector{Float64},
    depth::Vector{Int64}
)::Matrix
    # Filter down to site ids to be considered
    candidate_locs = depth[(k_area .> 0.0)[depth]]
    n_locs = length(candidate_locs)
    s_iv_locs = n_locs < n_iv_locs ? n_locs : n_iv_locs

    sel = StatsBase.sample(candidate_locs, s_iv_locs; replace = false)

    return [location_ids[sel] sel]
end


include("Criteria/DecisionPreferences.jl")
include("Criteria/DecisionWeights.jl")

end