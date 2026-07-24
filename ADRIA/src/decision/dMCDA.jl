const COUNTERFACTUAL_SCEN_ENCODING::Int64 = -1
const UNGUIDED_SCEN_ENCODING::Int64 = 0

# dummy functions to allow precompilation
function unguided_selection() end
function rank_sites!() end
function adria_topsis() end
function adria_vikor() end
function decision_matrix() end

"""
    mcda_methods()

List of MCDA methods found to be suitable for use with ADRIA.

See "Validate included MCDA methods" testset in `test/mcda.jl` for details.
"""
function mcda_methods()
    return [
        JMcDM.COCOSO.CocosoMethod,
        JMcDM.MAIRCA.MaircaMethod,
        JMcDM.MOORA.MooraMethod,
        JMcDM.PIV.PIVMethod,
        JMcDM.VIKOR.VikorMethod
    ]
end

"""
    mcda_method_names()::Vector{String}

List of MCDA method names.
"""
function mcda_method_names()::Vector{String}
    return string.(nameof.(parentmodule.(mcda_methods())))
end

"""
    mcda_method_encoding(mcda_name::String)::Union{Int64, Nothing}

Find the scenario encoding of the given mcda method.
Throws `ArgumentError` if method not found.
"""
function mcda_method_encoding(mcda_name::String)::Int64
    method_idx = findfirst(mcda_method_names() .== uppercase(mcda_name))
    if isnothing(method_idx)
        msg = "Given MCDA method $(mcda_name) is not in list of mcda methods.\n"
        msg *= "Possible MCDA methods are $(mcda_method_names())"
        throw(ArgumentError(msg))
    end
    return method_idx
end

"""
    decision_method_encoding(method_name::String)::Int64

Get the integer encoding for the type of intervention to be used in the scenario dataframe.
"""
function decision_method_encoding(method_name::String)::Int64
    method_name = uppercase(method_name)
    if method_name == "COUNTERFACTUAL"
        return COUNTERFACTUAL_SCEN_ENCODING
    elseif method_name == "UNGUIDED"
        return UNGUIDED_SCEN_ENCODING
    end

    return mcda_method_encoding(method_name)
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
    for (i, loc_id) in enumerate(s_order[:, 1])
        rankings[rankings[:, 1] .== loc_id, col] .= i
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

"""
    build_decay(plan_horizon::Int64, projection_confidence::Float64; α::Float64=0.99)::Vector{Float64}

Per-timestep decay vector weighting environmental projections by lead time, for
use with `weighted_projection`. `α` is the per-step base decay rate (unchanged
default). `projection_confidence` (0.0-1.0) sets `exponent = 2.0 -
projection_confidence`: 0.0 gives exponent 2.0 (steep decay, historical
default), 1.0 gives exponent 1.0 (mild decay). `plan_horizon` controls the
window length (vector length), not decay shape.
"""
function build_decay(
    plan_horizon::Int64, projection_confidence::Float64; α::Float64=0.99
)::Vector{Float64}
    exponent = 2.0 - projection_confidence
    return α .^ (1:(plan_horizon + 1)) .^ exponent
end

"""
    weighted_projection(
        env_data::AbstractMatrix{Float64}, tstep::Int64, planning_horizon::Int64,
        decay::AbstractVector{Float64}, timeframe::Int64
    )::Vector{Float64}

Projection of environmental data for locations present in `env_data` for the next
`planning_horizon` timesteps weighted by `decay`.
"""
function weighted_projection(
    env_data::AbstractMatrix,
    tstep::Int64,
    planning_horizon::Int64,
    decay::AbstractVector{Float64},
    timeframe::Int64
)::Vector{Float64}
    horizon = tstep:min(tstep + planning_horizon, timeframe)
    n = length(horizon)
    # Materialize once as a plain Matrix (avoids YAXArrays broadcast dispatch),
    # then scale each row in-place so no second full-matrix allocation is needed.
    env_slice = Matrix{Float64}(env_data[horizon, :])
    @inbounds for i = 1:n
        @views env_slice[i, :] .*= decay[i]
    end
    return summary_stat_env(env_slice, 1)
end

"""
    summary_stat_env(env_layer::AbstractArray, dims::Union{Int64,Symbol,Tuple{Symbol,Symbol}}; w=0.5)::Vector{Float64}

Calculates weighted combinations of mean and standard deviation for a given environmental
factor.

# Arguments
- `env_layer` : Environmental data layer to calculate the mean of.
- `dims` : Dimensions to aggregate over.
- `w` : Weight for mean value (the complement will be used to weight stdev).

# Returns
If the time horizon > 1, returns the weighted combination of mean and standard deviation of
the projected environmental conditions (e.g., DHWs, wave stress, etc):
    (μ * w) + (σ * (1 - w))

Where the time horizon == 1, the original values are returned.
"""
function summary_stat_env(
    env_layer::AbstractArray,
    dims::Union{Int64,Symbol,Tuple{Symbol,Symbol}};
    w=0.5
)::Vector{Float64}
    if size(env_layer, 1) > 1
        w1 = 1.0 - w
        μ = mean(env_layer; dims=dims)
        σ = std(env_layer; mean=μ, dims=dims)
        @. μ = μ * w + σ * w1
        return vec(μ)
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
    unguided_selection(
        n_iv_locs::Int64,
        k_area::Vector{Float64}
    )

Randomly select intervention locations, constraining to locations that are able to support
corals and are within the desired depth range.

If no depth range is provided, then simply selects from reefs with available space.

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
    # Filter down to location ids to be considered
    candidate_locs = findall((k_area .> 0.0) .& depth)
    n_locs = length(candidate_locs)
    s_iv_locs = n_locs < n_iv_locs ? n_locs : n_iv_locs

    sel = StatsBase.sample(candidate_locs, s_iv_locs; replace=false)

    return location_ids[sel]
end
function unguided_selection(
    location_ids,
    n_iv_locs::Int64,
    k_area::Vector{Float64}
)::Vector{<:Union{Symbol,String,Int64}}
    # Filter down to location ids to be considered
    candidate_locs = findall(k_area .> 0.0)
    n_locs = length(candidate_locs)
    s_iv_locs = n_locs < n_iv_locs ? n_locs : n_iv_locs

    sel = StatsBase.sample(candidate_locs, s_iv_locs; replace=false)

    return location_ids[sel]
end
