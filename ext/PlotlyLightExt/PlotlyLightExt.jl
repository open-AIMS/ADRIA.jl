module PlotlyLightExt

using PlotlyLight
using YAXArrays
using ADRIA
using ADRIA: ResultSet

# TODO: All common helper methods should be moved into main `viz` module.

"""
    _calc_confint(ci)

Helper method to determine quantile-based confidence bounds.

# Arguments
- `ci` : Desired quantile-based confidence bound (e.g., 0.95)

# Returns
Tuple, lower and upper CI bounds.
"""
function _calc_confint(ci::T)::Tuple{T,T} where {T<:AbstractFloat}
    d = (1.0 - ci) / 2.0
    return (d, 1.0 - d)
end

"""
    _get_guided_labels()::Vector{String}

Returns labels for categories of the `guided` factor.
"""
function _get_guided_labels()::Vector{String}
    return [
        "Counterfactual",
        "Semi-guided",
        last.(split.(string.(ADRIA.decision.mcda_methods()), "."))...
    ]
end

function outcome_label(
    outcomes::YAXArray; metadata_key::Symbol=:metric_feature, label_case=titlecase
)::String
    outcome_metadata = outcomes.properties

    _outcome_label = if all(haskey.([outcome_metadata], [metadata_key, :metric_unit]))
        _metric_feature = label_case(outcome_metadata[metadata_key])
        _metric_unit = outcome_metadata[:metric_unit]
        _metric_label = !isempty(_metric_unit) ? "[$(_metric_unit)]" : ""
        "$(_metric_feature) $(_metric_label)"
    else
        ""
    end

    return _outcome_label
end

"""
    _guided_colors()

Unique colors for scenario types.

TODO: Move to common theme module and make colorblind friendly.
"""
function _guided_colors()
    # Extracted "Set1" colorway from Plotly
    return [
        "#e41a1c",  # red
        "#377eb8",  # blue
        "#4daf4a",  # green
        "#984ea3",  # purple
        "#ff7f00",  # orange
        "#ffff33",  # yellow
        "#a65628",  # brown
        "#f781bf",  # pink
        "#999999"   # gray
    ]
end

include("./viz/scenarios.jl")

end