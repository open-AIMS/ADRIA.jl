using Statistics

# ──────────────────────────────────────────────────────────────────────────────
# data_envelopment_analysis(rs, dea_output) — wrapper for API compatibility
# ──────────────────────────────────────────────────────────────────────────────

"""
    ADRIA.viz.data_envelopment_analysis(rs::Union{Domain,ResultSet}, dea_output; kwargs...)

Plot DEA results. The `rs` argument is accepted for API compatibility with the Makie backend
but is not used by the Plotly implementation.
"""
function ADRIA.viz.data_envelopment_analysis(
    rs::Union{ADRIA.Domain,ADRIA.ResultSet},
    dea_output;
    metrics_x_lab::String="",
    metrics_y_lab::String="",
    frontier_type::Symbol=:vrs_peers,
    title::String="Data Envelopment Analysis",
    kwargs...
)::PlotlyBase.Plot
    return ADRIA.viz.data_envelopment_analysis(
        dea_output;
        metrics_x_lab=metrics_x_lab,
        metrics_y_lab=metrics_y_lab,
        frontier_type=frontier_type,
        title=title,
        kwargs...
    )
end

# ──────────────────────────────────────────────────────────────────────────────
# data_envelopment_analysis(dea_output)
# ──────────────────────────────────────────────────────────────────────────────

"""
    ADRIA.viz.data_envelopment_analysis(dea_output; metrics_x_lab="", metrics_y_lab="",
                                         frontier_type=:vrs_peers, kwargs...)

Plot DEA results in 3 subplot rows:
1. Output metric scatter with efficiency frontier
2. VRS efficiency scores
3. CRS efficiency scores
"""
function ADRIA.viz.data_envelopment_analysis(
    dea_output;
    metrics_x_lab::String="Metric X",
    metrics_y_lab::String="Metric Y",
    frontier_type::Symbol=:vrs_peers,
    title::String="Data Envelopment Analysis",
    kwargs...
)::PlotlyBase.Plot
    Y = dea_output.Y
    vrs_vals = dea_output.vrs_vals
    crs_vals = dea_output.crs_vals
    peers = getproperty(dea_output, frontier_type)

    n_scens = size(Y, 1)
    scen_idx = 1:n_scens

    # ── Row 1: output metric scatter with frontier highlighted ────────────────
    data_cloud = PlotlyBase.scatter(;
        x=Y[:, 1], y=Y[:, 2], mode="markers",
        marker_color=_hex_to_rgb(ADRIAviz.COLORS[:scenarios]),
        marker_symbol="circle", marker_size=6,
        name="Scenario data cloud",
        xaxis="x", yaxis="y",
        type="scatter"
    )

    peer_idx = peers.J
    frontier_trace = PlotlyBase.scatter(;
        x=Y[peer_idx, 1], y=Y[peer_idx, 2], mode="markers",
        marker_color=_hex_to_rgb(ADRIAviz.COLORS[:guided]),
        marker_symbol="star", marker_size=10,
        name="Best practice frontier",
        xaxis="x", yaxis="y",
        type="scatter"
    )

    # ── Row 2: VRS efficiency bar ─────────────────────────────────────────────
    vrs_bar = PlotlyBase.scatter(;
        x=collect(scen_idx), y=vrs_vals, mode="markers",
        marker_color=_hex_to_rgb(ADRIAviz.COLORS[:guided]),
        name="VRS efficiency",
        xaxis="x2", yaxis="y2",
        type="scatter"
    )

    # ── Row 3: CRS efficiency bar ─────────────────────────────────────────────
    crs_bar = PlotlyBase.scatter(;
        x=collect(scen_idx), y=crs_vals, mode="markers",
        marker_color=_hex_to_rgb(ADRIAviz.COLORS[:unguided]),
        name="CRS efficiency",
        xaxis="x3", yaxis="y3",
        type="scatter"
    )

    row_h = 0.3
    gap = 0.05
    fsz = _plotly_font_sizes(1)
    layout = PlotlyBase.Layout(;
        ADRIA_LAYOUT_DEFAULTS...,
        title=PlotlyBase.attr(; text=title, font=PlotlyBase.attr(; size=fsz.title)),
        font=PlotlyBase.attr(; family="Open Sans, sans-serif", size=fsz.label),
        # Row 1
        xaxis=PlotlyBase.attr(;
            title_text=metrics_x_lab,
            domain=[0.0, 1.0],
            anchor="y",
            tickfont=PlotlyBase.attr(; size=fsz.tick)
        ),
        yaxis=PlotlyBase.attr(;
            title_text=metrics_y_lab,
            domain=[2 * (row_h + gap), 1.0],
            anchor="x",
            tickfont=PlotlyBase.attr(; size=fsz.tick)
        ),
        # Row 2
        xaxis2=PlotlyBase.attr(;
            title_text="Scenario",
            domain=[0.0, 1.0],
            anchor="y2",
            tickfont=PlotlyBase.attr(; size=fsz.tick)
        ),
        yaxis2=PlotlyBase.attr(;
            title_text="VRS Efficiency",
            domain=[row_h + gap, 2 * row_h + gap],
            anchor="x2",
            tickfont=PlotlyBase.attr(; size=fsz.tick)
        ),
        # Row 3
        xaxis3=PlotlyBase.attr(;
            title_text="Scenario",
            domain=[0.0, 1.0],
            anchor="y3",
            tickfont=PlotlyBase.attr(; size=fsz.tick)
        ),
        yaxis3=PlotlyBase.attr(;
            title_text="CRS Efficiency",
            domain=[0.0, row_h],
            anchor="x3",
            tickfont=PlotlyBase.attr(; size=fsz.tick)
        )
    )

    return PlotlyBase.Plot([data_cloud, frontier_trace, vrs_bar, crs_bar], layout)
end
