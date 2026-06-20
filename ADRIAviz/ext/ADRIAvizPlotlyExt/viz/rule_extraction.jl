using DataFrames

# ──────────────────────────────────────────────────────────────────────────────
# rules_scatter(scenarios, clusters, rules)
#
# `rules` is a Vector of `ADRIAanalysis.Rule`. Each rule's `.condition` is a
# Vector of clauses, where each clause is [feature_name, operator, threshold].
# Only rules with exactly 2 clauses are plotted.
# ──────────────────────────────────────────────────────────────────────────────

function ADRIA.viz.rules_scatter(
    scenarios::DataFrame,
    clusters,
    rules::AbstractVector;
    title::String="Rule Scatter",
    kwargs...
)::PlotlyBase.Plot
    # Filter to rules with exactly 2 clauses
    two_clause = filter(r -> length(r.condition) == 2, rules)

    if isempty(two_clause)
        return PlotlyBase.Plot(
            PlotlyBase.AbstractTrace[],
            PlotlyBase.Layout(; ADRIA_LAYOUT_DEFAULTS..., title_text=title)
        )
    end

    # Resolve clusters -> BitVector of "target" membership
    target_mask = if clusters isa BitVector
        clusters
    else
        Bool.(clusters .!= 0)
    end
    non_target_mask = .!target_mask

    target_color = get(ADRIAviz.COLORS, :target, "#1f78b4")
    non_target_color = get(ADRIAviz.COLORS, :non_target, "#ff7f00")

    # Pre-filter to rules whose features exist in scenarios
    valid_rules = filter(two_clause) do r
        feat1 = string(first(r.condition[1]))
        feat2 = string(first(r.condition[2]))
        feat1 in names(scenarios) && feat2 in names(scenarios)
    end

    if isempty(valid_rules)
        return PlotlyBase.Plot(
            PlotlyBase.AbstractTrace[],
            PlotlyBase.Layout(; ADRIA_LAYOUT_DEFAULTS..., title_text=title)
        )
    end

    n_valid = length(valid_rules)
    domains, n_rows, n_cols = _grid_domains(n_valid)
    fsz = _plotly_font_sizes(n_valid)

    layout = PlotlyBase.Layout(;
        ADRIA_LAYOUT_DEFAULTS...,
        title=PlotlyBase.attr(; text=title, font=PlotlyBase.attr(; size=fsz.title)),
        font=PlotlyBase.attr(; family="Open Sans, sans-serif", size=fsz.label),
        legend=PlotlyBase.attr(; orientation="v", x=1.02, xanchor="left", y=1.0, title_text="Clusters"),
        width=max(700, 360 * n_cols),
        height=max(400, 320 * n_rows),
    )

    traces = PlotlyBase.AbstractTrace[]

    for (idx, rule) in enumerate(valid_rules)
        feat1 = string(first(rule.condition[1]))
        feat2 = string(first(rule.condition[2]))

        x_all = scenarios[:, feat1]
        y_all = scenarios[:, feat2]

        sfx = idx == 1 ? "" : string(idx)
        xd, yd, _, _ = domains[idx]

        layout[Symbol("xaxis$(sfx)")] = PlotlyBase.attr(;
            title_text=feat1, domain=xd, anchor="y$(sfx)"
        )
        layout[Symbol("yaxis$(sfx)")] = PlotlyBase.attr(;
            title_text=feat2, domain=yd, anchor="x$(sfx)"
        )

        push!(
            traces,
            PlotlyBase.scatter(;
                x=x_all[non_target_mask], y=y_all[non_target_mask],
                mode="markers",
                marker_color=_hex_to_rgb(non_target_color),
                marker_symbol="circle", marker_size=5,
                name="Non-target",
                xaxis="x$(sfx)", yaxis="y$(sfx)",
                legendgroup="Non-target", showlegend=(idx == 1),
                type="scatter"
            )
        )
        push!(
            traces,
            PlotlyBase.scatter(;
                x=x_all[target_mask], y=y_all[target_mask],
                mode="markers",
                marker_color=_hex_to_rgb(target_color),
                marker_symbol="diamond", marker_size=7,
                name="Target",
                xaxis="x$(sfx)", yaxis="y$(sfx)",
                legendgroup="Target", showlegend=(idx == 1),
                type="scatter"
            )
        )
    end

    return PlotlyBase.Plot(traces, layout)
end
