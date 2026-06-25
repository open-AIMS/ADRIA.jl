using StatsBase
using YAXArrays
import YAXArrays.DD: At
using ADRIA: Factor, DiscreteOrderedUniformDist, component_params

abstract type DecisionPreference end

struct DecisionPreferences <: DecisionPreference
    names::Vector{Symbol}
    weights::Vector{Float64}
    directions::Vector{Function}
end

"""
    decision_matrix(loc_names::Vector{T}, criteria_names::Vector{T}, criteria_vals::Matrix)::YAXArray where {T<:Union{String,Symbol}}
    decision_matrix(loc_names::Vector{T}, criteria_names::Vector{T2}; kwargs...)::YAXArray where {T<:Union{String,Symbol}}

Construct a decision matrix.

# Arguments
- `loc_names` : location names
- `criteria_names` : name of criteria being considered
- `criteria_vals` : values for each criteria
- `kwargs` : Preset criteria values by their names (with prefix removed)
"""
function decision_matrix(
    loc_names::Vector{T},
    criteria_names::Vector{T2},
    criteria_vals::Matrix
)::YAXArray where {T<:Union{String,Symbol},T2<:Union{String,Symbol}}
    return DataCube(criteria_vals; location=loc_names, criteria=criteria_names)
end
function decision_matrix(
    loc_names::Vector{T},
    criteria_names::Vector{T2};
    kwargs...
)::YAXArray where {T<:Union{String,Symbol},T2<:Union{String,Symbol}}
    mat = decision_matrix(
        loc_names,
        criteria_names,
        zeros(length(loc_names), length(criteria_names))
    )

    if length(kwargs) > 0
        update_criteria_values!(mat; kwargs...)
    end

    return mat
end

"""
    update_criteria_values!(dm::YAXArray, values::Matrix)::Nothing
    update_criteria_values!(dm::YAXArray; kwargs...)::Nothing

# Arguments
- `dm` : Decision matrix
- `values` : New criteria values to update decision matrix with
- `kwargs` : New values to assign for specified criteria

# Returns
Nothing
"""
function update_criteria_values!(dm::YAXArray, values::Matrix)::Nothing
    dm.data .= values

    return nothing
end
function update_criteria_values!(dm::YAXArray; kwargs...)::Nothing
    for (criteria_name, value) in kwargs
        dm[criteria=At(criteria_name)] .= value
    end

    return nothing
end

"""
    filter_criteria(prefs::T, remove::Vector)::T where {T<:DecisionPreference}

Filter criteria marked to be removed according to vector of true/false.

# Arguments
- `prefs` : The DecisionPreference to update.
- `remove` : Boolean vector indicating which columns to remove
"""
function filter_criteria(prefs::T, remove::Vector)::T where {T<:DecisionPreference}
    return typeof(prefs)(
        prefs.names[.!remove],
        prefs.weights[.!remove],
        prefs.directions[.!remove]
    )
end

"""
    solve(dp::T, dm::YAXArray, method::DataType) where {T<:DecisionPreference}
    solve(dp::T, dm::YAXArray, method::Function) where {T<:DecisionPreference}

# Arguments
- `dp` : DecisionPreference
- `dm` : the decision matrix to assess
- `method` : An MCDA method provided by the JMcDM package

# Returns
JMcDM result type (to confirm)
"""
function solve(
    dp::T, dm::YAXArray, method::DataType
) where {T<:DecisionPreference}
    return mcdm(MCDMSetting(dm.data, dp.weights, dp.directions), method())
end
function solve(
    dp::T, dm::YAXArray, method::Function
) where {T<:DecisionPreference}
    return method(dm.data, dp.weights, dp.directions)
end

"""
    rank_by_index(dp::T, dm::YAXArray, method::Union{Function,DataType})::Vector{Int64} where {T<:DecisionPreference}

Default index rank method, returns location indices in order of their rank.

Note: Ignores constant criteria values.

# Arguments
- `dp` : DecisionPreferences
- `dm` : The decision matrix to assess
- `method` : An MCDA method provided by the JMcDM package

# Returns
Index of locations ordered by their rank
"""
function rank_by_index(
    dp::T, dm::YAXArray, method::Union{Function,DataType}
)::Vector{Int64} where {T<:DecisionPreference}
    decision_results = criteria_aggregated_scores(dp, dm, method)
    is_maximal = decision_results.bestIndex == argmax(decision_results.scores)
    return sortperm(decision_results.scores; rev=is_maximal)
end

"""
    criteria_aggregated_scores(dp::T, dm::YAXArray, method::Function)::NamedTuple where {T<:DecisionPreference}

Calculates raw aggreagted score for a set of locations, in order of the location indices.

Note: Ignores constant criteria values.

# Arguments
- `dp` : DecisionPreferences
- `dm` : The decision matrix to assess
- `method` : An MCDA method provided by the JMcDM package

# Returns
Returns raw aggreagted score for a set of locations
"""
function criteria_aggregated_scores(
    dp::T, dm::YAXArray, method::Union{Function,DataType}
)::NamedTuple where {T<:DecisionPreference}
    # Identify valid, non-constant, columns for use in MCDA
    is_const = Bool[length(x) == 1 for x in unique.(eachcol(dm.data))]

    # YAXArrays will throw error for all false boolean masks
    if all(is_const)
        throw(DomainError(is_const, "No Ranking Possible, all criteria are constant"))
    end

    # Recreate preferences, removing criteria that are constant for this scenario
    if any(is_const)
        dp = filter_criteria(dp, is_const)
    end

    # Assess decision matrix only using valid (non-constant) criteria
    res = solve(dp, dm[criteria=.!is_const], method)

    if all(isnan.(res.scores))
        # This may happen if there are constants in the decision matrix
        # or if the method fails for some reason...
        throw(DomainError(res.scores, "No ranking possible"))
    end
    return (scores=res.scores, bestIndex=res.bestIndex)
end

"""
    select_locations(dp::T, dm::YAXArray, method::Union{Function,DataType}, min_locs::Int64)::Vector{<:Union{String,Symbol,Int64}} where {T<:DecisionPreference}

Default location selection method.
Returns the selected location names ordered by their rank.

# Arguments
- `dp` : DecisionPreferences
- `dm` : The decision matrix to assess
- `method` : An MCDA method provided by the JMcDM package
- `min_locs` : Minimum number of locations to select

# Returns
Name of selected locations ordered by their rank
"""
function select_locations(
    dp::T, dm::YAXArray, method::Union{Function,DataType}, min_locs::Int64
)::Vector{<:Union{String,Symbol,Int64}} where {T<:DecisionPreference}
    local rank_idx
    try
        rank_idx = rank_by_index(dp, dm, method)
    catch err
        if err isa DomainError
            # Return empty matrix to signify no ranks
            return String[]
        end

        rethrow(err)
    end

    min_locs = min(min_locs, length(dm.location[rank_idx]))
    return collect(dm.location[rank_idx][1:min_locs]).data
end

"""
    apply_threshold(criteria_name::Symbol, threshold::Tuple, dm::YAXArray}

Apply a threshold filter to a given decision matrix, filtering out locations that are
outside the given bounds.

# Arguments
- `criteria_name` : Criteria to apply thresholds to
- `threshold` : Lower and upper bounds of values to retain (inclusive)
- `dm` : The decision matrix to apply thresholds to

# Returns
A new YAXArray
"""
function apply_threshold(
    criteria_name::Union{Symbol,String}, threshold::Tuple, dm::YAXArray
)
    if criteria_name ∉ dm.criteria
        msg = "Criteria $(criteria_name) was not found when applying thresholds."
        msg *= "\nNo changes occurred."
        @warn msg
        return dm
    end

    target_criteria = dm.criteria .== criteria_name
    target_vals = dm.data[:, target_criteria]

    valid_locs = vec(threshold[1] .<= target_vals .<= (threshold[1] + threshold[2]))

    return dm[location=valid_locs]
end

"""
    cluster_diversity(cluster_ids::Vector{<:Union{Int64,String,Symbol}})::Vector{Float64}

Calculate how selecting each location would improve cluster diversity.
Higher scores indicate locations from underrepresented clusters.
"""
function cluster_diversity(
    cluster_ids::Vector{<:Union{Int64,String,Symbol}}
)::Vector{Float64}
    # Count unique clusters and their frequencies
    clusters = unique(cluster_ids)
    n_clusters = length(clusters)
    n_locations = length(cluster_ids)

    # Calculate the ideal count per cluster (perfect distribution)
    ideal_count = n_locations / n_clusters

    # Create a map of cluster to its count
    cluster_counts = countmap(cluster_ids)

    # Calculate diversity score for each unique cluster
    cluster_scores = Dict{eltype(clusters),Float64}()

    for cluster in clusters
        # Calculate how underrepresented this cluster is
        count_ratio = ideal_count / cluster_counts[cluster]

        # Scale so the maximum possible value is 1.0
        # A fully represented cluster gets 0.0
        # An underrepresented cluster gets a score based on its deficit
        if count_ratio > 1.0  # Underrepresented
            cluster_scores[cluster] = min(1.0, (count_ratio - 1.0) / (n_clusters - 1.0))
        else  # Overrepresented
            cluster_scores[cluster] = 0.0
        end
    end

    # Map cluster scores back to each location
    diversity_scores = zeros(n_locations)
    for i in 1:n_locations
        diversity_scores[i] = cluster_scores[cluster_ids[i]]
    end

    return diversity_scores
end

"""
    geographic_separation(mean_distances::Matrix{Float64})::Vector{Float64}

Calculate how spatially separated each location is from others.
Higher scores indicate locations that are more distant from their neighbors.
"""
function geographic_separation(mean_distances::Vector{Float64})::Vector{Float64}
    # Normalize to [0,1]
    max_dist = maximum(mean_distances)
    if max_dist > 0
        return mean_distances ./ max_dist
    end

    return zeros(length(mean_distances))
end

# """
#     _update_state(cluster_ids::Vector, num_locs::Int64, max_members::Int64)

# Update list of:
# - selected clusters
# - number of locations selected for each cluster
# - the indices of clusters which exceed the membership rule
# - IDs of the above
# - and potential suitable alternative locations

# # Arguments
# - `cluster_ids` : ID of clusters
# - `num_locs` : Number of locations to select
# - `max_members` : Maximum number of members per cluster
# """
# function _update_state(cluster_ids::Vector, num_locs::Int64, max_members::Int64)
#     selected_clusters = cluster_ids[1:min(num_locs, length(cluster_ids))]
#     cluster_frequency = countmap(selected_clusters)
#     rule_violators_idx = findall(values(cluster_frequency) .> max_members)
#     exceeded_clusters = collect(keys(cluster_frequency))[rule_violators_idx]
#     potential_alternatives = cluster_ids[cluster_ids .∉ Ref(exceeded_clusters)]

#     return selected_clusters,
#     cluster_frequency,
#     rule_violators_idx,
#     exceeded_clusters,
#     potential_alternatives
# end
