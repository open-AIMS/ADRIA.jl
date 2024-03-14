module decision

using InteractiveUtils: subtypes
using StatsBase
using YAXArrays
using ADRIA: DataCube, Domain, EcoModel, n_locations, site_area, site_k_area

using ADRIA:
    Domain,
    EcoModel,
    n_locations,
    site_area,
    site_k_area,
    component_params

using ADRIA:
    DiscreteOrderedUniformDist,
    Factor

using
    Combinatorics,
    DataFrames,
    JMcDM

# dummy dMCDA_vars() for dev
struct DMCDA_vars end

# dummy functions to allow precompilation
function unguided_selection() end
function rank_sites!() end
function adria_topsis() end
function adria_vikor() end
function decision_matrix() end

function supported_jmcdm_methods()
    # Those commented out with failed a simple test where only heat stress is considered
    # and the method selected the hottest locations ("**")
    # or a nonsensical selection (e.g., 1, 2, 3, 4, 5; marked with "*^")
    return [
        #** JMcDM.Topsis.TopsisMethod,
        JMcDM.VIKOR.VikorMethod,
        #** JMcDM.ARAS.ArasMethod,
        #*^ JMcDM.COCOSO.CocosoMethod,
        #** JMcDM.CODAS.CodasMethod,
        #** JMcDM.EDAS.EdasMethod,
        #** JMcDM.GREY.GreyMethod,
        #** JMcDM.MABAC.MabacMethod,
        JMcDM.MAIRCA.MaircaMethod,
        #** JMcDM.MARCOS.MarcosMethod,
        #*^ JMcDM.MOORA.MooraMethod,
        JMcDM.PIV.PIVMethod,
        JMcDM.PSI.PSIMethod,
        #** JMcDM.SAW.SawMethod,
        #** JMcDM.WASPAS.WaspasMethod,
        #** JMcDM.WPM.WPMMethod
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
    return x ./ sqrt.(sum(x .^ 2; dims=1))
end
function mcda_normalize(x::DataFrame)::DataFrame
    return x ./ sum(Matrix(x); dims=2)
end

"""
    align_rankings!(rankings::Array, s_order::Matrix, col::Int64)::Nothing

Align a vector of site rankings to match the indicated order in `s_order`.
"""
function align_rankings!(rankings::Array, s_order::Matrix, col::Int64)::Nothing
    # Fill target ranking column
    for (i, site_id) in enumerate(s_order[:, 1])
        rankings[rankings[:, 1].==site_id, col] .= i
    end

    return nothing
end

"""
    identify_within_depth_bounds(
        loc_depths::Vector{Float64},
        min_depth::T,
        offset::T
    )::BitVector where {T<:Union{Int64,Float64}}

Identify locations that are in the desired depth bounds.

# Returns
BitVector, of logical indices indicating locations which satisfy the depth criteria.
"""
function identify_within_depth_bounds(
    loc_depths::Vector{Float64},
    min_depth::Float64,
    offset::Float64
)::BitVector
    max_depth = (min_depth + offset)
    @assert min_depth <= max_depth "Minimum depth must be lower than maximum depth"

    # Default to using all locations if constant, otherwise apply bounds.
    depth_criteria::BitArray{1} = fill(true, length(loc_depths))
    if .!all(loc_depths .== loc_depths[1])
        depth_criteria = (loc_depths .>= min_depth) .& (loc_depths .<= max_depth)
    end

    return depth_criteria
end

function weighted_projection(env_data, tstep, planning_horizon, decay, timeframe)
    # Determine subset of data to select data for planning horizon
    horizon::UnitRange{Int64} = tstep:min(tstep + planning_horizon, timeframe)
    d_s::UnitRange{Int64} = 1:length(horizon)

    # Put more weight on projected conditions closer to the decision point
    @views env_horizon = decay[d_s] .* env_data[horizon, :]
    return summary_stat_env(env_horizon, :timesteps)
end

"""
    summary_stat_env(env_layer::AbstractArray, dims::Union{Int64, Symbol, Tuple{Symbol, Symbol}}; w=0.5)::Vector{Float64}

Calculates weighted combinations of mean and standard deviation for a given environmental
factor.

# Arguments
- `env_layer` : Environmental data layer to calculate the mean of.
- `dims` : Dimensions to aggregate over.
- `w` : Weight for mean value (the complement will be used to weight stdev).

# Returns
If the time horizon > 1, returns the weighted combination of mean and standard deviation of the projected environmental
conditions (e.g., DHWs, wave stress, etc):
    (μ * w) + (σ * (1 - w))

Where the time horizon == 1, the original values are returned.
"""
function summary_stat_env(
    env_layer::AbstractArray,
    dims::Union{Int64,Symbol,Tuple{Symbol,Symbol}};
    w=0.5
)::Vector{Float64}
    if size(env_layer, 1) > 1
        return vec((mean(env_layer; dims=dims) .* w) .+ (std(env_layer; dims=dims) .* (1.0 - w)))
    end

    return vec(env_layer)
end

"""
    decision_frequency(
        start_year::Int64,
        timeframe::Int64,
        n_years::Int64,
        freq::Union{Float64,Int64}
    )::Vector{Bool}

# Arguments
- `start_year` : Year to begin deployments
- `timeframe` : Total length of simulation time
- `n_years` : Number of years to deploy for
- `freq` : Frequency (in years) of decision/deployments

# Returns
Vector of true/false indicating which years in simulation period to apply a decision.
"""
function decision_frequency(
    start_year::Int64, timeframe::Int64, n_years::Int64, freq::Union{Float64,Int64}
)::Vector{Bool}
    freq_timeframe = fill(false, timeframe)

    start_year = max(start_year, 2)
    if freq > 0
        max_consider = min(start_year + n_years - 1, timeframe)
        freq_timeframe[start_year:Int64(freq):max_consider] .= true
    else
        # Start at year 2 or the given specified start year
        freq_timeframe[start_year] = true
    end

    return freq_timeframe
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
    depth::BitVector
)::Vector{<:Union{Symbol,String,Int64}}
    # Filter down to site ids to be considered
    candidate_locs = findall((k_area .> 0.0) .& depth)
    n_locs = length(candidate_locs)
    s_iv_locs = n_locs < n_iv_locs ? n_locs : n_iv_locs

    sel = StatsBase.sample(candidate_locs, s_iv_locs; replace=false)

    return location_ids[sel]
end

include("Criteria/DecisionPreferences.jl")
include("Criteria/DecisionWeights.jl")
include("location_selection.jl")

export
    SeedPreferences,
    FogPreferences,
    SRMPreferences,
    decision_matrix,
    filter_criteria,
    update_criteria_values!,
    select_locations,
    unguided_selection,
    decision_frequency,
    weighted_projection,
    identify_within_depth_bounds,
    summary_stat_env,
    mcda_methods
end
