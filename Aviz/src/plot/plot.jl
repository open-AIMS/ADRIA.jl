module plot

using GLMakie, DataFrames

using ADRIA: ResultSet, metrics.metric_label, timesteps, sensitivity, analysis.col_normalize
using ADRIA.NamedDims, ADRIA.AxisKeys
using Aviz: scenario_type, scenario_colors, COLORS
using GLMakie.Colors


function _time_labels(labels)
    tick_pos = vcat(1:5:length(labels), [length(labels)])
    tick_label = vcat(string.(labels)[1:5:end], string.(labels)[end])

    return tick_pos, tick_label
end

include("scenario.jl")
include("sensitivity.jl")

export scenario
export pawn

end
