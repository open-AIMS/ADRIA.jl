using YAXArrays
using DataFrames
using Statistics
using ADRIA.analysis: series_confint
using ADRIA: functional_group_names, human_readable_name

# ──────────────────────────────────────────────────────────────────────────────
# taxonomy(rs::ResultSet)
# ──────────────────────────────────────────────────────────────────────────────

function ADRIA.viz.taxonomy(
    rs::ADRIA.ResultSet;
    kwargs...
)::PlotlyBase.Plot
    if !haskey(rs.outcomes, :relative_taxa_cover)
        throw(
            ArgumentError(
                "Unable to find relative_taxa_cover in outcomes. Pass it manually."
            )
        )
    end
    return ADRIA.viz.taxonomy(rs.inputs, rs.outcomes[:relative_taxa_cover]; kwargs...)
end

# ──────────────────────────────────────────────────────────────────────────────
# taxonomy(scenarios, relative_taxa_cover)
# ──────────────────────────────────────────────────────────────────────────────

"""
    ADRIA.viz.taxonomy(scenarios::DataFrame, relative_taxa_cover::YAXArray; kwargs...)

Plot relative coral cover by taxon group over time.
`relative_taxa_cover` has dims `(timesteps x groups x scenarios)`.
"""
function ADRIA.viz.taxonomy(
    scenarios::DataFrame,
    relative_taxa_cover::YAXArray;
    by_RCP::Bool=false,
    by_functional_groups::Bool=true,
    show_confints::Bool=true,
    alpha::Real=0.15,
    xlabel::String="Year",
    ylabel::String="Relative Cover",
    title::String="",
    kwargs...
)::PlotlyBase.Plot
    _scens = scenarios[collect(relative_taxa_cover.scenarios), :]
    scen_groups = if by_RCP
        _scenario_rcps(_scens)
    else
        _scenario_types(_scens)
    end

    x_vals = collect(relative_taxa_cover.timesteps)
    tickvals, ticktext = _year_ticks(x_vals)
    n_groups = length(relative_taxa_cover.groups)
    mat_3d = collect(relative_taxa_cover)  # (timesteps x groups x scenarios) -- avoids repeated cat
    n_scen_groups = length(scen_groups)

    # Colour palette for inner loop (groups or scen_groups)
    fallback_colors = ["#377eb8", "#ff7f00", "#4daf4a", "#984ea3", "#e41a1c",
        "#a65628", "#f781bf", "#999999"]

    fg_names = functional_group_names()
    taxa_names = if n_groups == length(fg_names)
        human_readable_name(fg_names; title_case=true)
    else
        ["Group $(g)" for g in 1:n_groups]
    end

    n_subplots = if by_functional_groups
        n_scen_groups
    else
        n_groups
    end

    fsz = _plotly_font_sizes(n_subplots)

    # Compute vertical-stack (1-col) domains for n_subplots panels.
    # Panel 1 is at the top; panels fill downward.
    y_gap = 0.08
    panel_h = (1.0 - y_gap * max(n_subplots - 1, 0)) / max(n_subplots, 1)

    legend_y = -0.05 - 0.03 * n_subplots

    layout = PlotlyBase.Layout(;
        ADRIA_LAYOUT_DEFAULTS...,
        title=PlotlyBase.attr(; text=title, font=PlotlyBase.attr(; size=fsz.title)),
        font=PlotlyBase.attr(; family="Open Sans, sans-serif", size=fsz.label),
        legend=PlotlyBase.attr(; orientation="h", x=0.5, xanchor="center", y=legend_y, yanchor="top"),
    )

    for i in 1:n_subplots
        sfx = i == 1 ? "" : string(i)
        y1 = 1.0 - (i - 1) * (panel_h + y_gap)
        y0 = y1 - panel_h
        layout[Symbol("xaxis$(sfx)")] = PlotlyBase.attr(;
            title_text=xlabel, domain=[0.0, 1.0], anchor="y$(sfx)",
            tickvals=tickvals, ticktext=ticktext
        )
        layout[Symbol("yaxis$(sfx)")] = PlotlyBase.attr(;
            title_text=ylabel, domain=[y0, y1], anchor="x$(sfx)"
        )
    end

    traces = PlotlyBase.AbstractTrace[]

    if by_functional_groups
        # One panel per scenario group; coloured by functional group index
        for (si, (_, smask)) in enumerate(scen_groups)
            any(smask) || continue
            sfx = si == 1 ? "" : string(si)

            for g = 1:n_groups
                subset = mat_3d[:, g, smask]
                if ndims(subset) == 1
                    subset = reshape(subset, :, 1)
                end
                color = fallback_colors[mod1(g, length(fallback_colors))]
                label = taxa_names[g]
                if show_confints
                    ci = series_confint(subset)
                    lo, mid, hi = ci[:, 1], ci[:, 2], ci[:, 3]
                    ts = _confint_traces(
                        x_vals,
                        lo,
                        mid,
                        hi;
                        name=label,
                        color=color,
                        alpha=alpha
                    )
                    for trace in ts
                        trace[:xaxis] = "x$(sfx)"
                        trace[:yaxis] = "y$(sfx)"
                    end
                    # Only show legend on first subplot's line trace
                    ts[2][:showlegend] = (si == 1)
                    ts[2][:legendgroup] = label
                    append!(traces, ts)
                else
                    push!(
                        traces,
                        PlotlyBase.scatter(;
                            x=x_vals, y=vec(mean(subset; dims=2)),
                            mode="lines", line_color=_hex_to_rgb(color),
                            name=label, xaxis="x$(sfx)", yaxis="y$(sfx)",
                            legendgroup=label, showlegend=(si == 1),
                            type="scatter"
                        )
                    )
                end
            end
        end
    else
        # One panel per functional group; coloured by scenario group
        colors = _group_colors(scen_groups)

        for (idx, g) in enumerate(1:n_groups)
            sfx = idx == 1 ? "" : string(idx)

            for (sname, smask) in scen_groups
                any(smask) || continue
                subset = mat_3d[:, g, smask]
                if ndims(subset) == 1
                    subset = reshape(subset, :, 1)
                end
                color = get(colors, sname, "#4682b4")
                label = string(sname)
                if show_confints
                    ci = series_confint(subset)
                    lo, mid, hi = ci[:, 1], ci[:, 2], ci[:, 3]
                    ts = _confint_traces(
                        x_vals, lo, mid, hi; name=label, color=color, alpha=alpha
                    )
                    for trace in ts
                        trace[:xaxis] = "x$(sfx)"
                        trace[:yaxis] = "y$(sfx)"
                    end
                    # Only show legend on first subplot's line trace
                    ts[2][:showlegend] = (idx == 1)
                    ts[2][:legendgroup] = label
                    append!(traces, ts)
                else
                    push!(
                        traces,
                        PlotlyBase.scatter(;
                            x=x_vals, y=vec(mean(subset; dims=2)),
                            mode="lines", line_color=_hex_to_rgb(color),
                            name=label, xaxis="x$(sfx)", yaxis="y$(sfx)",
                            legendgroup=label, showlegend=(idx == 1),
                            type="scatter"
                        )
                    )
                end
            end
        end
    end

    return PlotlyBase.Plot(traces, layout)
end
