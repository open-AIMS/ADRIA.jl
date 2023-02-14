module plot

using GLMakie, DataFrames

using ADRIA: ResultSet, metrics.metric_label, timesteps, sensitivity, analysis.col_normalize
using ADRIA.NamedDims, ADRIA.AxisKeys
using Aviz: scenario_type, scenario_colors, COLORS
using GLMakie.Colors


"""
    _time_labels(labels)

Extract time step labels, ensuring last entry is always included.
"""
function _time_labels(labels)
    tick_pos = vcat(1:5:length(labels), [length(labels)])
    tick_label = vcat(string.(labels)[1:5:end], string.(labels)[end])

    return tick_pos, tick_label
end

function _dimkeys(outcomes::NamedDimsArray)
    return (; zip(dimnames(outcomes), axiskeys(outcomes))...)
end

"""
    timesteps(outcomes::NamedDimsArray)::Array{Int64}

Extract time step labels from outcome arrays.

# Arguments
- `outcomes` : Results to extract metadata from

# Returns
Array of time steps (years)
"""
function timesteps(outcomes::NamedDimsArray)::Array{Int64}
    axis_labels = _dimkeys(outcomes)

    if :timesteps in keys(axis_labels)
        return axis_labels.timesteps
    end

    return Int64[]
end


include("scenario.jl")
include("sensitivity.jl")

export scenario
export pawn

end
