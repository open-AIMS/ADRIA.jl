using ADRIA.analysis: col_normalize
using OrderedCollections: OrderedDict

const _ADRIA_ANALYSIS_PKG_ID = Base.PkgId(
    Base.UUID("bc4d0ea8-6565-4397-854d-28474bf8c6b3"), "ADRIAanalysis"
)

function relative_sensitivities(
    X, y::AbstractArray{<:Real}; S=10, stat=:median
)::Vector{Float64}
    if !haskey(Base.loaded_modules, _ADRIA_ANALYSIS_PKG_ID)
        error(
            "relative_sensitivities requires ADRIAanalysis: run `using ADRIAanalysis` first"
        )
    end
    sensitivity = getproperty(Base.loaded_modules[_ADRIA_ANALYSIS_PKG_ID], :sensitivity)
    return col_normalize(sensitivity.pawn(X, y; S=S)(; Si=stat))
end

"""
    outcome_probability(data::AbstractVector)::NamedTuple

Determine probability occurrence.
"""
function outcome_probability(data::AbstractVector)::NamedTuple
    p_outcomes = cdf.(fit(Distributions.Normal, data), data)

    n = length(data)
    return (
        values=[
            count(p_outcomes .> 0.80) / n,
            count((p_outcomes .> 0.70) .& (p_outcomes .<= 0.80)) / n,
            count((p_outcomes .>= 0.50) .& (p_outcomes .<= 0.70)) / n,
            count((p_outcomes .> 0.20) .& (p_outcomes .< 0.50)) / n,
            count(p_outcomes .< 0.20) / n],
        labels=[
            "Very High\n> 80%",
            "High\n70 - 80%",
            "Medium\n50 - 70%",
            "Low\n20 - 50%",
            "Very Low\n< 20%"
        ]
    )
end

"""
    _get_scenario_groups(ao::AnnotatedOutcomes; by_RCP=false) -> OrderedDict{Symbol,BitVector}

Extract scenario group masks from `ao.metadata`, keyed by RCP or scenario type.
"""
function _get_scenario_groups(
    ao::AnnotatedOutcomes; by_RCP::Bool=false
)::OrderedDict{Symbol,BitVector}
    if by_RCP
        groups = get(ao.metadata, :scenario_rcp_groups, nothing)
        isnothing(groups) && throw(ArgumentError(
            "RCP grouping is not available for RME-sourced outcomes."
        ))
        return groups
    end
    haskey(ao.metadata, :scenario_type_groups) || throw(
        ArgumentError(
            "AnnotatedOutcomes is missing :scenario_type_groups — was attach_scenario_metadata called?"
        )
    )
    return ao.metadata[:scenario_type_groups]
end
