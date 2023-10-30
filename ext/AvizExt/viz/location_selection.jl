using AxisKeys, NamedDims
using ADRIA: ResultSet
"""
    ADRIA.viz.selection_frequency_map!(g::Union{GridLayout,GridPosition},
        rs::ResultSet, iv_type::String; scen_ids::Vector{Int64}=collect(1:size(rs.inputs, 1)),
        opts::Dict=Dict(:color_map => [:red, :blue], :colorbar_label => "Selection frequency"),
        axis_opts::Dict=Dict())
    ADRIA.viz.selection_frequency_map(rs::ResultSet, iv_type::String;
        scen_ids::Vector{Int64}=collect(1:size(rs.inputs, 1)),
        opts::Dict=Dict(:color_map => [:red, :blue], :colorbar_label => "Selection frequency"),
        fig_opts::Dict=Dict(), axis_opts::Dict=Dict())

Plot a spatial map of location selection frequencies.

# Arguments
- `rs` : Result set.
- `iv_type` : Intervention type (e.g. "seed" or "shade").
- `scen_ids` : Subset of scenarios to plot (could be robust scenarios, or all scenarios)
- `opts` : Aviz options
    - `colorbar_label`, label for colorbar. Defaults to "Relative Cover".
    - `color_map`, preferred colormap for plotting heatmaps.
- `axis_opts` : Additional options to pass to adjust Axis attributes
  See: https://docs.makie.org/v0.19/api/index.html#Axis

# Returns
Figure
"""
function ADRIA.viz.ranks_to_frequencies!(
    g::Union{GridLayout,GridPosition},
    rs::ResultSet,
    frequencies::NamedDimsArray,
    rank_ids::Vector{Int64};
    opts::Dict=Dict(),
    axis_opts::Dict=Dict(),
)
    rank_groups = Dict(Symbol(rank_grp) => rank_grp .== rank_ids for rank_grp in rank_ids)
    rank_colors = colors(rank_groups)
    for rr in rank_ids
        color_map = Dict(:color_map => [RGBf(1.0, 1.0, 1.0), rank_colors[Symbol(rr)]])
        merge!(opts, color_map)
        ADRIA.viz.map!(
            g,
            rs,
            AxisKeys.keyless(NamedDims.unname(frequencies[ranks=rr]));
            opts=opts,
            axis_opts=axis_opts,
        )
    end
    return g
end
function ADRIA.viz.ranks_to_frequencies!(
    g::Union{GridLayout,GridPosition},
    rs::ResultSet,
    frequencies::NamedDimsArray,
    rank_id::Int64;
    opts::Dict=Dict(:color_map => :CMRmap),
    axis_opts::Dict=Dict())

    return ADRIA.viz.map!(
        g,
        rs,
        AxisKeys.keyless(NamedDims.unname(frequencies[ranks=rank_id]));
        opts=opts,
        axis_opts=axis_opts,
    )

end
function ADRIA.viz.ranks_to_frequencies(
    rs::ResultSet,
    frequencies::NamedDimsArray,
    rank_ids::Union{Int64,Vector{Int64}};
    opts::Dict=Dict(),
    fig_opts::Dict=Dict(), axis_opts::Dict=Dict())

    f = Figure(; fig_opts...)
    g = f[1, 1] = GridLayout()

    return ADRIA.viz.ranks_to_frequencies!(
        g,
        rs,
        frequencies,
        rank_ids;
        opts=opts,
        axis_opts=axis_opts,
    )
end

function _default_colormap(rank_groups::Dict, alpha_vals::Dict)
    rank_colors = colors(rank_groups, alpha_vals)
    rank_ids = keys(rank_groups)
    return Dict(
        rank_grp =>
            [RGBA{Float32}(1.0, 1.0, 1.0, alpha_vals[rank_grp]), rank_colors[rank_grp]] for
        rank_grp in rank_ids
    )
end