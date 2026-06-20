using Statistics
using ADRIA.analysis: col_normalize

# ──────────────────────────────────────────────────────────────────────────────
# pawn(Si) — heatmap of sensitivity indices
# ──────────────────────────────────────────────────────────────────────────────

function ADRIA.viz.pawn(
    Si::YAXArray;
    normalize::Bool=true,
    factors=:all,
    by::Symbol=:median,
    title::String="PAWN Sensitivity",
    kwargs...
)::PlotlyBase.Plot
    if normalize
        Si = YAXArray(Si.axes, col_normalize(Si.data))
    end

    if factors != :all
        Si = Si[factors = At(factors)]
    end

    # Sort rows by chosen statistic (descending)
    sort_vals = collect(Si[Si = At(by)])
    perm = sortperm(sort_vals; rev=true)
    Si = Si[perm, :]

    factor_names = string.(collect(Si.factors))
    si_names = string.(collect(Si.Si))
    z = collect(Si)  # (n_factors × n_Si)

    trace = PlotlyBase.heatmap(;
        z=z,
        x=si_names,
        y=factor_names,
        colorscale="Blues",
        type="heatmap"
    )

    # Scale height with factor count so labels stay legible for large factor sets.
    fig_height = max(400, 80 + 22 * length(factor_names))
    fsz = _plotly_font_sizes(1)

    layout = PlotlyBase.Layout(;
        ADRIA_LAYOUT_DEFAULTS...,
        title=PlotlyBase.attr(; text=title, font=PlotlyBase.attr(; size=fsz.title)),
        font=PlotlyBase.attr(; family="Open Sans, sans-serif", size=fsz.label),
        height=fig_height,
        xaxis=PlotlyBase.attr(;
            title_text="SI Statistic", tickfont=PlotlyBase.attr(; size=fsz.tick)
        ),
        yaxis=PlotlyBase.attr(;
            title_text="Factor", ticktext=factor_names, tickvals=factor_names,
            tickfont=PlotlyBase.attr(; size=fsz.tick),
            # Rows are sorted by `by` descending, so the most influential factor is
            # `factor_names[1]`. Plotly places y[1] at the bottom, so reverse the
            # axis to put the most influential factor at the top.
            autorange="reversed", automargin=true
        )
    )
    return PlotlyBase.Plot([trace], layout)
end

# ──────────────────────────────────────────────────────────────────────────────
# tsa(Si) — time-varying sensitivity indices
# ──────────────────────────────────────────────────────────────────────────────

function ADRIA.viz.tsa(
    Si::YAXArray;
    stat::Symbol=:median,
    xlabel::String="Year",
    ylabel::String="Sensitivity Index",
    title::String="Temporal Sensitivity",
    kwargs...
)::PlotlyBase.Plot
    factor_names = collect(Si.factors)
    timesteps_vals = collect(Si.timesteps)
    tickvals, ticktext = _year_ticks(timesteps_vals)

    fallback_colors = ["#377eb8", "#ff7f00", "#4daf4a", "#984ea3", "#e41a1c",
        "#a65628", "#f781bf", "#999999"]

    Si_stats = collect(Si.Si)
    stat_idx = findfirst(==(stat), Si_stats)
    mat = collect(Si)  # (n_factors × n_Si × n_timesteps) — avoids DimensionalData cat warnings

    traces = PlotlyBase.AbstractTrace[]
    for (i, fname) in enumerate(factor_names)
        y_vals = mat[i, stat_idx, :]
        color = fallback_colors[mod1(i, length(fallback_colors))]
        push!(
            traces,
            PlotlyBase.scatter(;
                x=timesteps_vals, y=y_vals, mode="lines+markers",
                line_color=_hex_to_rgb(color), name=string(fname), type="scatter"
            )
        )
    end

    fsz = _plotly_font_sizes(1)

    layout = PlotlyBase.Layout(;
        ADRIA_LAYOUT_DEFAULTS...,
        title=PlotlyBase.attr(; text=title, font=PlotlyBase.attr(; size=fsz.title)),
        font=PlotlyBase.attr(; family="Open Sans, sans-serif", size=fsz.label),
        xaxis=PlotlyBase.attr(;
            title_text=xlabel,
            tickvals=tickvals,
            ticktext=ticktext,
            tickfont=PlotlyBase.attr(; size=fsz.tick)
        ),
        yaxis=PlotlyBase.attr(;
            title_text=ylabel, tickfont=PlotlyBase.attr(; size=fsz.tick)
        )
    )
    return PlotlyBase.Plot(traces, layout)
end

# ──────────────────────────────────────────────────────────────────────────────
# rsa(Si, factor_values) — response surface for each factor
# ──────────────────────────────────────────────────────────────────────────────

function ADRIA.viz.rsa(
    Si::YAXArray,
    factor_values::AbstractMatrix;
    xlabel::String="Factor Value",
    ylabel::String="Sensitivity Index",
    title::String="Regional Sensitivity",
    kwargs...
)::PlotlyBase.Plot
    factor_names = collect(Si.factors)
    n_factors = length(factor_names)
    quantile_pos = collect(Si.si_quantile)
    mat = collect(Si)  # (n_factors × n_si_quantile) — avoids DimensionalData cat warnings

    fallback_colors = ["#377eb8", "#ff7f00", "#4daf4a", "#984ea3", "#e41a1c",
        "#a65628", "#f781bf", "#999999"]

    traces = PlotlyBase.AbstractTrace[]

    for (i, fname) in enumerate(factor_names)
        si_vals = mat[i, :]
        fv_col = factor_values[:, i]
        fv_range = quantile(fv_col, quantile_pos)
        color = fallback_colors[mod1(i, length(fallback_colors))]
        axis_sfx = i == 1 ? "" : string(i)

        push!(
            traces,
            PlotlyBase.scatter(;
                x=fv_range, y=si_vals, mode="lines+markers",
                line_color=_hex_to_rgb(color), name=string(fname),
                xaxis="x$(axis_sfx)", yaxis="y$(axis_sfx)",
                type="scatter"
            )
        )
    end

    fsz = _plotly_font_sizes(n_factors)
    layout = _grid_layout(
        factor_names,
        n_factors;
        title=title,
        xlabel=xlabel,
        ylabel=ylabel,
        fsz=fsz
    )
    return PlotlyBase.Plot(traces, layout)
end

# helper: compute the [x_domain, y_domain] rectangles for an n-panel subplot grid.
# Panels fill left-to-right, top-to-bottom (panel 1 is top-left). Reuses
# `_calc_gridsize` so the grid shape matches the rest of ADRIAviz.
function _grid_domains(n::Int; x_gap::Float64=0.08, y_gap::Float64=0.12)
    n_rows, n_cols = _calc_gridsize(n)
    w = (1.0 - x_gap * (n_cols - 1)) / n_cols
    h = (1.0 - y_gap * (n_rows - 1)) / n_rows
    domains = Tuple{Vector{Float64},Vector{Float64},Int,Int}[]
    for i = 1:n
        row = div(i - 1, n_cols) + 1   # row 1 is the top row
        col = mod(i - 1, n_cols) + 1
        x0 = (col - 1) * (w + x_gap)
        y1 = 1.0 - (row - 1) * (h + y_gap)
        y0 = y1 - h
        push!(domains, ([x0, x0 + w], [y0, y1], row, col))
    end
    return domains, n_rows, n_cols
end

# helper: build a subplot-grid Layout for the per-factor sensitivity plots.
# Each panel gets its own x/y axes (rotated x ticks + automargin so labels stay
# legible for large factor sets), a shared `xlabel`/`ylabel`, and a factor-name
# subplot title placed above each panel.
function _grid_layout(
    factor_names,
    n_factors::Int;
    title::String,
    xlabel::String,
    ylabel::String,
    fsz=(title=12, label=12, tick=12)
)
    domains, n_rows, n_cols = _grid_domains(n_factors)
    annotations = PlotlyBase.PlotlyAttribute[]
    layout = PlotlyBase.Layout(;
        ADRIA_LAYOUT_DEFAULTS...,
        title=PlotlyBase.attr(; text=title, font=PlotlyBase.attr(; size=fsz.title)),
        font=PlotlyBase.attr(; family="Open Sans, sans-serif", size=fsz.label),
        width=max(700, 360 * n_cols),
        height=max(400, 320 * n_rows),
        showlegend=false
    )
    for i = 1:n_factors
        xd, yd, _, col = domains[i]
        sfx = i == 1 ? "" : string(i)
        layout[Symbol("xaxis$(sfx)")] = PlotlyBase.attr(;
            title_text=xlabel, domain=xd, anchor="y$(sfx)",
            tickangle=-45, automargin=true, tickfont=PlotlyBase.attr(; size=fsz.tick)
        )
        layout[Symbol("yaxis$(sfx)")] = PlotlyBase.attr(;
            title_text=(col == 1 ? ylabel : ""), domain=yd, anchor="x$(sfx)",
            automargin=true, tickfont=PlotlyBase.attr(; size=fsz.tick)
        )
        push!(
            annotations,
            PlotlyBase.attr(;
                text=string(factor_names[i]), x=(xd[1] + xd[2]) / 2, y=yd[2],
                xref="paper", yref="paper", xanchor="center", yanchor="bottom",
                showarrow=false, font=PlotlyBase.attr(; size=fsz.label)
            )
        )
    end
    layout[:annotations] = annotations
    return layout
end

# ──────────────────────────────────────────────────────────────────────────────
# outcome_map(outcomes, factor_values)
# ──────────────────────────────────────────────────────────────────────────────

function ADRIA.viz.outcome_map(
    outcomes::YAXArray,
    factor_values::AbstractMatrix;
    xlabel::String="Factor Value",
    ylabel::String="Outcome",
    title::String="Outcome Map",
    alpha::Real=0.2,
    kwargs...
)::PlotlyBase.Plot
    factor_names = collect(outcomes.factors)
    n_factors = length(factor_names)
    quantile_pos = collect(outcomes.si_quantile)
    CI_vals = collect(outcomes.CI)
    lo_idx = findfirst(==(:lower), CI_vals)
    mid_idx = findfirst(==(:mean), CI_vals)
    hi_idx = findfirst(==(:upper), CI_vals)
    mat = collect(outcomes)  # (n_factors × n_CI × n_si_quantile) — avoids DimensionalData cat warnings

    fallback_colors = ["#377eb8", "#ff7f00", "#4daf4a", "#984ea3", "#e41a1c",
        "#a65628", "#f781bf", "#999999"]

    traces = PlotlyBase.AbstractTrace[]

    for (i, fname) in enumerate(factor_names)
        lo = mat[i, lo_idx, :]
        mid = mat[i, mid_idx, :]
        hi = mat[i, hi_idx, :]
        fv_col = factor_values[:, i]
        fv_range = quantile(fv_col, quantile_pos)
        color = fallback_colors[mod1(i, length(fallback_colors))]
        axis_sfx = i == 1 ? "" : string(i)

        fill_x = vcat(vec(fv_range), reverse(vec(fv_range)))
        fill_y = vcat(hi, reverse(lo))
        push!(
            traces,
            PlotlyBase.scatter(;
                x=fill_x, y=fill_y, fill="toself",
                fillcolor=_hex_to_rgba(color, alpha),
                line_color="rgba(0,0,0,0)",
                showlegend=false, hoverinfo="skip",
                name=string(fname) * "_band",
                xaxis="x$(axis_sfx)", yaxis="y$(axis_sfx)",
                type="scatter"
            )
        )
        push!(
            traces,
            PlotlyBase.scatter(;
                x=fv_range, y=mid, mode="lines",
                line_color=_hex_to_rgb(color), line_width=2,
                name=string(fname),
                xaxis="x$(axis_sfx)", yaxis="y$(axis_sfx)",
                type="scatter"
            )
        )
    end

    fsz = _plotly_font_sizes(n_factors)
    layout = _grid_layout(
        factor_names,
        n_factors;
        title=title,
        xlabel=xlabel,
        ylabel=ylabel,
        fsz=fsz
    )
    return PlotlyBase.Plot(traces, layout)
end

# ──────────────────────────────────────────────────────────────────────────────
# convergence(Si_conv, factors) — sensitivity convergence lines
# ──────────────────────────────────────────────────────────────────────────────

function ADRIA.viz.convergence(
    Si_conv::YAXArray,
    factors::AbstractVector;
    plot_overlay::Bool=true,
    stat::Symbol=:median,
    xlabel::String="N Scenarios",
    ylabel::String="Sensitivity Index",
    title::String="Sensitivity Convergence",
    kwargs...
)::PlotlyBase.Plot
    sample_sizes = collect(Si_conv.n_scenarios)
    Si_stats = collect(Si_conv.Si)
    stat_idx = findfirst(==(stat), Si_stats)
    mat = collect(Si_conv)  # (n_factors × n_Si × n_samples) — avoids DimensionalData cat warnings
    fallback_colors = ["#377eb8", "#ff7f00", "#4daf4a", "#984ea3", "#e41a1c",
        "#a65628", "#f781bf", "#999999"]

    traces = PlotlyBase.AbstractTrace[]
    for (i, fname) in enumerate(factors)
        y_vals = mat[i, stat_idx, :]
        color = fallback_colors[mod1(i, length(fallback_colors))]
        push!(
            traces,
            PlotlyBase.scatter(;
                x=sample_sizes, y=y_vals, mode="lines+markers",
                line_color=_hex_to_rgb(color), name=string(fname), type="scatter"
            )
        )
    end

    fsz = _plotly_font_sizes(1)

    layout = PlotlyBase.Layout(;
        ADRIA_LAYOUT_DEFAULTS...,
        title=PlotlyBase.attr(; text=title, font=PlotlyBase.attr(; size=fsz.title)),
        font=PlotlyBase.attr(; family="Open Sans, sans-serif", size=fsz.label),
        xaxis=PlotlyBase.attr(;
            title_text=xlabel, tickfont=PlotlyBase.attr(; size=fsz.tick)
        ),
        yaxis=PlotlyBase.attr(;
            title_text=ylabel, tickfont=PlotlyBase.attr(; size=fsz.tick)
        )
    )
    return PlotlyBase.Plot(traces, layout)
end
