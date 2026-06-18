using DataFrames, DimensionalData, YAXArrays
using ADRIA: axes_names, ResultSet, model_spec

const OPT_TYPE = Dict{Symbol,<:Any}
const DEFAULT_OPT_TYPE = Dict{Symbol,Any}

# Export shared spatial utilities for use in extensions
export _loc_id_col, _get_site_ids, _site_ids, _haversine_km, _nice_length
export GBR_COASTAL_PLACES, CoastalPlace, MapDecorationData, compute_map_decorations

function _no_backend_error()
    return error(
        "No visualization backend loaded. Load a backend before calling viz functions:\n" *
        "  using GLMakie      # interactive desktop\n" *
        "  using WGLMakie     # Pluto / browser\n" *
        "  using CairoMakie   # static PNG/SVG"
    )
end

"""
    _time_labels(labels)

Extract time step labels, ensuring last entry is always included.
"""
function _time_labels(labels; label_step=5)::Tuple{Vector{Int64},Vector{String}}
    labels_length = length(labels)
    labels_strings = string.(labels)

    tick_position = collect(1:label_step:labels_length)
    tick_label = collect(labels_strings[1:label_step:end])

    if (labels_length - 1) % label_step != 0
        return vcat(tick_position, labels_length), vcat(tick_label, labels_strings[end])
    end

    return tick_position, tick_label
end

"""
    timesteps(outcomes::YAXArray)::Vector{Int64}

Extract time step labels from outcome arrays. Delegates to `ADRIA.timesteps`.
"""
timesteps(outcomes::YAXArrays.YAXArray) = ADRIA.timesteps(outcomes)

function timesteps(outcomes::AbstractMatrix)::UnitRange{Int64}
    return 1:size(outcomes, 1)
end

"""
    _calc_gridsize(n_factors::Int64; max_cols::Int64=4)::Tuple{Int64,Int64}

Calculates a "nice" number of rows and columns from a given number of factors.
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

include("../outcome_metadata.jl")
include("spatial_utils.jl")
