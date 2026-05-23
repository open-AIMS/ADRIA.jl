using YAXArrays
using ADRIA.analysis: series_confint

# ──────────────────────────────────────────────────────────────────────────────
# scenarios(matrix, clusters)
# ──────────────────────────────────────────────────────────────────────────────

"""
    ADRIA.viz.scenarios(outcomes::AbstractMatrix, clusters; summarize=true, kwargs...)

Plot scenario outcomes grouped by cluster assignment.
`outcomes` is `(timesteps × scenarios)`, `clusters` is an integer or BitVector.
"""
function ADRIA.viz.scenarios(
    outcomes::AbstractMatrix,
    clusters;
    summarize::Bool=true,
    alpha::Real=0.15,
    xlabel::String="Year",
    ylabel::String="",
    title::String="",
    kwargs...
)::PlotlyBase.Plot
    groups = _scenario_clusters(clusters)
    x_vals = 1:size(outcomes, 1)
    tickvals, ticktext = _year_ticks(x_vals)
    colors = _group_colors(groups)

    traces = PlotlyBase.AbstractTrace[]
    for (gname, mask) in groups
        any(mask) || continue
        subset = outcomes[:, mask]
        if ndims(subset) == 1
            subset = reshape(subset, :, 1)
        end
        color = get(colors, gname, "#4682b4")
        label = string(gname)

        if summarize
            ci = series_confint(subset)
            lo, mid, hi = ci[:, 1], ci[:, 2], ci[:, 3]
            ts = _confint_traces(
                collect(x_vals), lo, mid, hi; name=label, color=color, alpha=alpha
            )
            append!(traces, ts)
        else
            for j in axes(subset, 2)
                push!(
                    traces,
                    PlotlyBase.scatter(;
                        x=collect(x_vals), y=subset[:, j], mode="lines",
                        line_color=_hex_to_rgba(color, 0.5),
                        name=label, showlegend=j == 1, type="scatter"
                    )
                )
            end
        end
    end

    layout = PlotlyBase.Layout(;
        ADRIA_LAYOUT_DEFAULTS...,
        title_text=title,
        xaxis=PlotlyBase.attr(;
            title_text=xlabel, tickvals=collect(tickvals), ticktext=ticktext
        ),
        yaxis=PlotlyBase.attr(; title_text=ylabel)
    )
    return PlotlyBase.Plot(traces, layout)
end

# ──────────────────────────────────────────────────────────────────────────────
# clustered_scenarios(matrix, clusters)
# ──────────────────────────────────────────────────────────────────────────────

"""
    ADRIA.viz.clustered_scenarios(outcomes::AbstractMatrix, clusters; kwargs...)

Plot all individual scenario time-series coloured by cluster.
"""
function ADRIA.viz.clustered_scenarios(
    outcomes::AbstractMatrix,
    clusters;
    alpha::Real=0.4,
    xlabel::String="Year",
    ylabel::String="",
    title::String="Cluster Overview",
    kwargs...
)::PlotlyBase.Plot
    groups = _scenario_clusters(clusters)
    x_vals = 1:size(outcomes, 1)
    tickvals, ticktext = _year_ticks(x_vals)
    colors = _group_colors(groups)

    traces = PlotlyBase.AbstractTrace[]
    for (gname, mask) in groups
        any(mask) || continue
        subset = outcomes[:, mask]
        if ndims(subset) == 1
            subset = reshape(subset, :, 1)
        end
        color = get(colors, gname, "#4682b4")
        label = string(gname)

        for j in axes(subset, 2)
            push!(
                traces,
                PlotlyBase.scatter(;
                    x=collect(x_vals), y=subset[:, j], mode="lines",
                    line_color=_hex_to_rgba(color, alpha),
                    name=label, showlegend=j == 1, type="scatter"
                )
            )
        end
    end

    layout = PlotlyBase.Layout(;
        ADRIA_LAYOUT_DEFAULTS...,
        title_text=title,
        xaxis=PlotlyBase.attr(;
            title_text=xlabel, tickvals=collect(tickvals), ticktext=ticktext
        ),
        yaxis=PlotlyBase.attr(; title_text=ylabel)
    )
    return PlotlyBase.Plot(traces, layout)
end
