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
    decision_matrix(loc_names::Vector{T}, criteria_names::Vector{T2}; kwargs...)::YAXArray where {T<:Union{String,Symbol}}
    decision_matrix(loc_names::Vector{T}, criteria_names::Vector{T}, criteria_vals::Matrix)::YAXArray where {T<:Union{String,Symbol}}

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

# Arguments
- `dm` : decision matrix
- `values` : new criteria values to update decision matrix with

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
    filter_criteria(prefs::T, is_const::Vector)::T where {T<:DecisionPreference}

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
    solve(dp::T, dm::YAXArray, method::Function) where {T<:DecisionPreference}

# Arguments
- `dp` : DecisionPreferences
- `dm` : the decision matrix to assess
- `method` : An MCDA method provided by the JMcDM package

# Returns
JMcDM result type (to confirm)
"""
function solve(
    dp::T, dm::YAXArray, method::Union{DataType}
) where {T<:DecisionPreference}
    return mcdm(MCDMSetting(dm.data, dp.weights, dp.directions), method())
end
function solve(
    dp::T, dm::YAXArray, method::Union{Function}
) where {T<:DecisionPreference}
    return method(dm.data, dp.weights, dp.directions)
end

"""
    rank_by_index(dp::T, dm::YAXArray, method::Function)::Vector{Int64} where {T<:DecisionPreference}

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
    # Identify valid, non-constant, columns for use in MCDA
    is_const = Bool[length(x) == 1 for x in unique.(eachcol(dm.data))]

    # Recreate preferences, removing criteria that are constant for this scenario
    _dp = filter_criteria(dp, is_const)

    # Assess decision matrix only using valid (non-constant) criteria
    res = solve(_dp, dm[criteria=.!is_const], method)

    scores = res.scores
    if all(isnan.(scores))
        # This may happen if there are constants in the decision matrix
        # or if the method fails for some reason...
        throw(DomainError(scores, "No ranking possible"))
    end

    is_maximal = res.bestIndex == argmax(res.scores)
    return sortperm(res.scores; rev=is_maximal)
end

"""
    select_locations(dp::T, dm::YAXArray, method::Function)::Vector{<:Union{String, Symbol}} where {T<:DecisionPreference}

Default location selection method.
Returns the selected location names ordered by their rank.

# Arguments
- `dp` : DecisionPreferences
- `dm` : The decision matrix to assess
- `method` : An MCDA method provided by the JMcDM package
- `min_locs` : Minimum number of locations to select

# Returns
Index of locations ordered by their rank
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

    return collect(dm.location[rank_idx][1:min_locs])
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

    valid_locs = vec(threshold[1] .<= target_vals .<= (threshold[1]+threshold[2]))

    return dm[location=valid_locs]
end
