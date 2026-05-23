using YAXArrays
using DataFrames
using Statistics
using ADRIA.analysis: series_confint

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
`relative_taxa_cover` has dims `(timesteps × groups × scenarios)`.
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
    mat_3d = collect(relative_taxa_cover)  # (timesteps × groups × scenarios) — avoids repeated cat
    n_scen_groups = length(scen_groups)

    # Colour palette for inner loop (groups or scen_groups)
    fallback_colors = ["#377eb8", "#ff7f00", "#4daf4a", "#984ea3", "#e41a1c",
        "#a65628", "#f781bf", "#999999"]

    traces = PlotlyBase.AbstractTrace[]

    if by_functional_groups
        # One panel per scenario group; coloured by functional group index
        for (si, (sname, smask)) in enumerate(scen_groups)
            any(smask) || continue
            for g in 1:n_groups
                subset = mat_3d[:, g, smask]
                if ndims(subset) == 1
                    subset = reshape(subset, :, 1)
                end
                color = _hex_to_rgb(fallback_colors[mod1(g, length(fallback_colors))])
                label = "$(sname) - group$(g)"
                if show_confints
                    ci = series_confint(subset)
                    lo, mid, hi = ci[:, 1], ci[:, 2], ci[:, 3]
                    ts = _confint_traces(
                        x_vals,
                        lo,
                        mid,
                        hi;
                        name=label,
                        color=fallback_colors[mod1(g, length(fallback_colors))],
                        alpha=alpha
                    )
                    append!(traces, ts)
                else
                    push!(
                        traces,
                        PlotlyBase.scatter(;
                            x=x_vals, y=vec(mean(subset; dims=2)),
                            mode="lines", line_color=color,
                            name=label, type="scatter"
                        )
                    )
                end
            end
        end
    else
        # Coloured by scenario group; separate line per group
        colors = _group_colors(scen_groups)
        for g in 1:n_groups
            for (sname, smask) in scen_groups
                any(smask) || continue
                subset = mat_3d[:, g, smask]
                if ndims(subset) == 1
                    subset = reshape(subset, :, 1)
                end
                color = get(colors, sname, "#4682b4")
                label = "$(sname) - group$(g)"
                if show_confints
                    ci = series_confint(subset)
                    lo, mid, hi = ci[:, 1], ci[:, 2], ci[:, 3]
                    ts = _confint_traces(
                        x_vals, lo, mid, hi; name=label, color=color, alpha=alpha
                    )
                    append!(traces, ts)
                else
                    push!(
                        traces,
                        PlotlyBase.scatter(;
                            x=x_vals, y=vec(mean(subset; dims=2)),
                            mode="lines", line_color=_hex_to_rgb(color),
                            name=label, type="scatter"
                        )
                    )
                end
            end
        end
    end

    layout = PlotlyBase.Layout(;
        ADRIA_LAYOUT_DEFAULTS...,
        title_text=title,
        xaxis=PlotlyBase.attr(;
            title_text=xlabel, tickvals=tickvals, ticktext=ticktext
        ),
        yaxis=PlotlyBase.attr(; title_text=ylabel)
    )
    return PlotlyBase.Plot(traces, layout)
end
