using AxisKeys, NamedDims
using ADRIA: location_selection_frequencies
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
function ADRIA.viz.selection_frequency_map!(
    g::Union{GridLayout,GridPosition},
    rs::ResultSet, iv_type::String; scen_ids::Vector{Int64}=collect(1:size(rs.inputs, 1)),
    opts::Dict=Dict(:color_map => [:red, :blue], :colorbar_label => "Selection frequency"),
    axis_opts::Dict=Dict())

    loc_frequencies = location_selection_frequencies(rs, iv_type; n_loc_int=rs.sim_constants["n_site_int"], ind_metrics=scen_ids)
    ADRIA.viz.map!(g, rs, AxisKeys.keyless(NamedDims.unname(loc_frequencies)); opts=opts, axis_opts=axis_opts)
end
function ADRIA.viz.selection_frequency_map(
    rs::ResultSet,
    iv_type::String;
    scen_ids::Vector{Int64}=collect(1:size(rs.inputs, 1)),
    opts::Dict=Dict(:color_map => [:red, :blue], :colorbar_label => "Selection frequency"),
    fig_opts::Dict=Dict(), axis_opts::Dict=Dict())

    f = Figure(; fig_opts...)
    g = f[1, 1] = GridLayout()

    return ADRIA.viz.selection_frequency_map!(
        g, rs, iv_type; scen_ids=scen_ids, opts=opts, axis_opts=axis_opts
    )
end
