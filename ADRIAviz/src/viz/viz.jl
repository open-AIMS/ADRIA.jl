module viz

using DataFrames, DimensionalData, Statistics, YAXArrays
using ADRIA
using ADRIA: axes_names, ResultSet, model_spec

const OPT_TYPE = Dict{Symbol,<:Any}
const DEFAULT_OPT_TYPE = Dict{Symbol,Any}

# Font sizing tier system for typography consistency across all Makie visualizations.
# 10pt is the minimum tick size for readability at 300 DPI (print). For large multi-panel
# layouts (5+), short numeric formats (e.g., ×10³ notation) are preferred over long decimals.
function set_typography_defaults!(axis_opts::OPT_TYPE; n_panels::Int=1)::OPT_TYPE
    if n_panels <= 1
        title_sz, label_sz, tick_sz = 16, 12, 10
    else
        title_sz, label_sz, tick_sz = 14, 12, 10
    end
    axis_opts[:titlesize] = get(axis_opts, :titlesize, title_sz)
    axis_opts[:xlabelsize] = get(axis_opts, :xlabelsize, label_sz)
    axis_opts[:ylabelsize] = get(axis_opts, :ylabelsize, label_sz)
    axis_opts[:xticklabelsize] = get(axis_opts, :xticklabelsize, tick_sz)
    axis_opts[:yticklabelsize] = get(axis_opts, :yticklabelsize, tick_sz)
    return axis_opts
end

function _no_backend_error()
    return ArgumentError(
        "No visualization backend loaded. Load a backend with `ADRIAviz.activate()` before calling viz functions:\n" *
        "  GLMakie      # interactive desktop\n" *
        "  WGLMakie     # Pluto / browser\n" *
        "  CairoMakie   # static PNG/SVG \n" *
        "  Plotly       # interactive browser-based"
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

Produces balanced, landscape-oriented grids:
- n=1 → (1,1), n=2 → (1,2), n=3–4 → (2,2), n=5–6 → (2,3), n≥7 → square-ish
"""
function _calc_gridsize(n_factors::Int64; max_cols::Int64=4)::Tuple{Int64,Int64}
    n_factors < 1 && error("n_factors must be >= 1")

    if n_factors == 1
        return 1, 1
    elseif n_factors == 2
        return 1, 2
    elseif n_factors <= 6
        n_rows = Int(ceil(sqrt(Float64(n_factors))))
        n_cols = Int(ceil(n_factors / n_rows))
        if n_rows > n_cols
            n_rows, n_cols = n_cols, n_rows
        end
        return n_rows, n_cols
    else
        n_side = Int(ceil(sqrt(Float64(n_factors))))
        return n_side, n_side
    end
end

include("../outcome_metadata.jl")
include("spatial_utils.jl")

"""
    _outcome_mask(outcome_threshold, y) -> BitVector

Convert the `outcome_threshold` keyword to a `BitVector` marking "behavioural" scenarios.
Dispatch on type keeps call sites free of branching.

- `nothing`                  : scenarios where `y > median(y)`
- `Real`                     : scenarios where `y > threshold`
- `AbstractVector{Bool}`     : passed through unchanged
- `AbstractVector{<:Integer}`: indices converted to a boolean mask
"""
_outcome_mask(::Nothing, y::AbstractVector{Float64}) = y .> median(y)
_outcome_mask(t::Real, y::AbstractVector{Float64}) = y .> Float64(t)
_outcome_mask(m::AbstractVector{Bool}, ::AbstractVector{Float64}) = m
function _outcome_mask(idx::AbstractVector{<:Integer}, y::AbstractVector{Float64})
    mask = falses(length(y))
    mask[idx] .= true
    return mask
end

"""
    _empirical_cdf(v) -> (sorted_values, cumulative_probs)

Return the sorted values and their empirical cumulative probabilities for a vector `v`.
"""
function _empirical_cdf(v::AbstractVector{<:Real})
    sv = sort(v)
    n = length(sv)
    return sv, collect((1:n) ./ n)
end

# Export shared spatial utilities for use in extensions (after they're defined)
export set_typography_defaults!
export _loc_id_col, _get_site_ids, _site_ids, _haversine_km, _nice_length, validate_extent
export GBR_COASTAL_PLACES, CoastalPlace, MapDecorationData, compute_map_decorations

# Export outcome metadata functions and types
export outcome_title, outcome_label, set_plot_opts!
export OPT_TYPE, DEFAULT_OPT_TYPE
export _time_labels, _calc_gridsize, timesteps
export _outcome_mask, _empirical_cdf

end  # module viz
