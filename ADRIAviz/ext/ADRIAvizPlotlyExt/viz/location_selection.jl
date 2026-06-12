using YAXArrays
using ADRIA: ResultSet

"""
    ADRIA.viz.ranks_to_frequencies(rs::ResultSet, frequencies::YAXArray,
                                   rank_ids::Union{Int,Vector{Int}};
                                   colorscale="Blues", title="", width=900, height=600)

Choropleth map of site-selection frequencies for one or more rank thresholds.
Each rank gets its own subplot panel. `frequencies` is a YAXArray with a
`ranks` dimension, as returned by `ADRIA.ranks_to_frequencies`.

# Arguments
- `rank_ids`  : Single rank integer or vector of rank integers to plot
- `colorscale`: Plotly colorscale name (default `"Blues"`)
- `title`     : Overall figure title
"""
function ADRIA.viz.ranks_to_frequencies(
    rs::ResultSet,
    frequencies::YAXArray,
    rank_ids::Union{Int64,Vector{Int64}};
    colorscale::String="Blues",
    title::String="Selection Frequency by Rank",
    width::Int=900,
    height::Int=600,
    kwargs...
)
    ids = rank_ids isa Int64 ? [rank_ids] : rank_ids
    n = length(ids)

    if n == 1
        freq_vec = collect(Float64, frequencies[ranks=ids[1]])
        return ADRIA.viz.map(
            rs,
            freq_vec;
            colorbar_label="Selection Frequency",
            colorscale=colorscale,
            title=isempty(title) ? "Rank $(ids[1])" : title,
            width=width,
            height=height
        )
    end

    # Multi-panel: one subplot per rank
    n_locs = length(rs.loc_data.site_id)
    outputs_matrix = hcat([collect(Float64, frequencies[ranks=r]) for r in ids]...)
    map_titles = ["Rank $(r)" for r in ids]

    return ADRIA.viz.map(
        rs,
        outputs_matrix,
        map_titles;
        colorscale=colorscale,
        colorbar_label="Selection Frequency",
        title=title,
        width=width,
        height=max(height, 400 * ceil(Int, n / 2))
    )
end

"""
    ADRIA.viz.selection_criteria_map(rs::Union{Domain,ResultSet},
                                     decision_matrix::YAXArray,
                                     scores::Vector{Float64};
                                     criteria=Array(decision_matrix.criteria),
                                     colorscale="Viridis", title="", width=1200, height=600)

Grid of choropleth maps showing the spatial distribution of each location
selection criterion, plus a final panel for the aggregate score.

# Arguments
- `decision_matrix`: YAXArray of size `(n_locs, n_criteria)`
- `scores`         : Aggregate per-location scores from `decision_matrix`
- `criteria`       : Criterion names (defaults to YAXArray dimension labels)
"""
function ADRIA.viz.selection_criteria_map(
    rs::Union{Domain,ResultSet},
    decision_matrix::YAXArray,
    scores::Vector{Float64};
    criteria::Vector{Symbol}=Array(decision_matrix.criteria),
    colorscale::String="Viridis",
    title::String="",
    width::Int=1200,
    height::Int=600,
    kwargs...
)
    if length(rs.loc_data.site_id) != size(decision_matrix, 1)
        error("Only unfiltered decision matrices can be plotted.")
    end

    # Append aggregate score as the last column
    outputs_matrix = hcat(Matrix(decision_matrix), scores)
    map_titles = vcat(
        [ADRIA.human_readable_name(string(c); title_case=true) for c in criteria],
        ["Aggregate Score"]
    )

    return ADRIA.viz.map(
        rs,
        outputs_matrix,
        map_titles;
        colorscale=colorscale,
        title=title,
        width=width,
        height=height
    )
end
