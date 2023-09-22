using GLMakie, DataFrames

using ADRIA: ResultSet, metrics.metric_label, analysis.col_normalize, model_spec
using NamedDims, AxisKeys
using .AvizExt: scenario_type, scenario_colors, COLORS
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

"""
    _calc_gridsize(n_factors::Int64; max_cols::Int64=4)::Tuple{Int64,Int64}

Calculates a "nice" number of rows and columns from a given number of factors to display.
The number of rows for subplots are calculated based on the number of desired columns.

Note: `n_factors` == 1 is displayed as a single figure.
      `n_factors` <= 4 are always displayed with 2 columns.

# Arguments
- `n_factors` : Number of factors to organize in a grid.

# Returns
Number of rows and columns
"""
function _calc_gridsize(n_factors::Int64; max_cols::Int64=4)::Tuple{Int64,Int64}
    if n_factors <= 4
        if n_factors == 1
            return 1, 1
        end

        n_cols::Int64 = 2
    else
        n_cols = max_cols
    end

    n_rows::Int64 = ceil(Int64, n_factors / n_cols)

    return n_rows, n_cols
end

include("scenarios.jl")
include("sensitivity.jl")
include("clustering.jl")
include("rule_extraction.jl")
include("spatial.jl")
