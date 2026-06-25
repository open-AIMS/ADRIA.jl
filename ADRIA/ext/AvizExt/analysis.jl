using ADRIA.analysis: col_normalize
using ADRIA.sensitivity

function relative_sensitivities(
    X, y::AbstractArray{<:Real}; S=10, stat=:median
)::Vector{Float64}
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
