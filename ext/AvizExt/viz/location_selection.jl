using AxisKeys, NamedDims
using ADRIA: location_selection_frequencies
"""
    loc_selection_frequency_map!(g::Union{GridLayout,GridPosition},
        rs::ResultSet, iv_type::String; scen_ids::Vector{Int64}=collect(1:size(rs.inputs, 1)),
        fig_opts::Dict=Dict(), axis_opts::Dict=Dict())
    loc_selection_frequency_map(rs::ResultSet, iv_type::String;
        scen_ids::Vector{Int64}=collect(1:size(rs.inputs, 1)),
        fig_opts::Dict=Dict(), axis_opts::Dict=Dict())

Plot a spatial map of location selection frequencies.

# Arguments
- `rs` : Result set.
- `iv_type` : Intervention type (e.g. "seed" or "shade").
- `scen_ids` : Subset of scenarios to plot (could be robust scenarios, or all scenarios)

# Returns
Figure
"""
function ADRIA.viz.loc_selection_frequency_map!(g::Union{GridLayout,GridPosition},
    rs::ResultSet, iv_type::String; scen_ids::Vector{Int64}=collect(1:size(rs.inputs, 1)),
    fig_opts::Dict=Dict(), axis_opts::Dict=Dict())

    loc_frequencies = location_selection_frequencies(rs, iv_type; n_loc_int=rs.sim_constants["n_site_int"], ind_metrics=scen_ids)
    ADRIA.viz.map!(g, rs, AxisKeys.keyless(NamedDims.unname(loc_frequencies)); axis_opts=axis_opts,
        color_map=[:red, :blue])
end
function ADRIA.viz.loc_selection_frequency_map(rs::ResultSet, iv_type::String;
    scen_ids::Vector{Int64}=collect(1:size(rs.inputs, 1)),
    fig_opts::Dict=Dict(), axis_opts::Dict=Dict())

    f = Figure(; fig_opts...)
    g = f[1, 1] = GridLayout()

    ADRIA.viz.loc_selection_frequency_map!(g, rs, iv_type; scen_ids=scen_ids, fig_opts=fig_opts, axis_opts=axis_opts)
end
