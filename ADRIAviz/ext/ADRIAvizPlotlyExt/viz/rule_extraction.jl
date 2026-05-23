using DataFrames

# ──────────────────────────────────────────────────────────────────────────────
# rules_scatter(scenarios, clusters, rules)
#
# Rules format: Vector{Vector{Vector}} — each rule is Vector{Vector},
# each clause is [feature_name, operator, threshold].
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
    two_clause = filter(r -> length(r) == 2, rules)

    if isempty(two_clause)
        return PlotlyBase.Plot(
            PlotlyBase.AbstractTrace[],
            PlotlyBase.Layout(; ADRIA_LAYOUT_DEFAULTS..., title_text=title)
        )
    end

    # Resolve clusters → BitVector of "target" membership
    target_mask = if clusters isa BitVector
        clusters
    else
        Bool.(clusters .!= 0)
    end
    non_target_mask = .!target_mask

    target_color = get(ADRIAviz.COLORS, :target, "#1f78b4")
    non_target_color = get(ADRIAviz.COLORS, :non_target, "#ff7f00")

    traces = PlotlyBase.AbstractTrace[]

    for rule in two_clause
        clause1, clause2 = rule[1], rule[2]
        feat1 = string(first(clause1))
        feat2 = string(first(clause2))

        if !(feat1 in names(scenarios)) || !(feat2 in names(scenarios))
            continue
        end

        x_all = scenarios[:, feat1]
        y_all = scenarios[:, feat2]

        push!(
            traces,
            PlotlyBase.scatter(;
                x=x_all[non_target_mask], y=y_all[non_target_mask],
                mode="markers",
                marker_color=_hex_to_rgb(non_target_color),
                marker_symbol="circle", marker_size=5,
                name="Non-target",
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
                type="scatter"
            )
        )
    end

    layout = PlotlyBase.Layout(;
        ADRIA_LAYOUT_DEFAULTS...,
        title_text=title,
        xaxis=PlotlyBase.attr(; title_text="Feature 1"),
        yaxis=PlotlyBase.attr(; title_text="Feature 2")
    )
    return PlotlyBase.Plot(traces, layout)
end
