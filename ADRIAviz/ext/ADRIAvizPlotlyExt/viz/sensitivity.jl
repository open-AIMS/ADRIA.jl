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
# rsa(X, y, foi) — 2D scatter of two factors colored by outcome
# ──────────────────────────────────────────────────────────────────────────────

"""
    ADRIA.viz.rsa(X, y, foi; with_contour=true, title, kwargs...) -> PlotlyBase.Plot

2D scatter of two input factors colored by outcome metric, with optional binned-average
contour overlay (histogram2dcontour with histfunc="avg").

# Arguments
- `X`            : Feature matrix (DataFrame, columns = factors)
- `y`            : Outcome values per scenario
- `foi`          : NTuple{2,Symbol} -- (x-axis factor, y-axis factor)
- `with_contour` : Add a histogram2dcontour overlay showing mean outcome per bin.
"""
function ADRIA.viz.rsa(
    X::DataFrame,
    y::AbstractVector{<:Real},
    foi::NTuple{2,Symbol};
    with_contour::Bool=true,
    title::String="Regional Sensitivity",
    kwargs...
)::PlotlyBase.Plot
    x_vals = Float64.(X[!, foi[1]])
    y_vals = Float64.(X[!, foi[2]])
    z_vals = Float64.(y)

    fsz = _plotly_font_sizes(1)

    traces = PlotlyBase.AbstractTrace[]

    if with_contour
        push!(
            traces,
            PlotlyBase.histogram2dcontour(;
                x=x_vals,
                y=y_vals,
                z=z_vals,
                histfunc="avg",
                colorscale="Viridis",
                showscale=false,
                opacity=0.4,
                ncontours=15,
                contours=PlotlyBase.attr(; coloring="fill", showlines=false),
                type="histogram2dcontour"
            )
        )
    end

    push!(
        traces,
        PlotlyBase.scatter(;
            x=x_vals,
            y=y_vals,
            mode="markers",
            marker=PlotlyBase.attr(;
                color=z_vals,
                colorscale="Viridis",
                showscale=true,
                opacity=0.6,
                size=5,
                colorbar=PlotlyBase.attr(; title="outcome")
            ),
            type="scatter"
        )
    )

    layout = PlotlyBase.Layout(;
        ADRIA_LAYOUT_DEFAULTS...,
        title=PlotlyBase.attr(; text=title, font=PlotlyBase.attr(; size=fsz.title)),
        font=PlotlyBase.attr(; family="Open Sans, sans-serif", size=fsz.label),
        xaxis=PlotlyBase.attr(;
            title_text=string(foi[1]), tickfont=PlotlyBase.attr(; size=fsz.tick)
        ),
        yaxis=PlotlyBase.attr(;
            title_text=string(foi[2]), tickfont=PlotlyBase.attr(; size=fsz.tick)
        )
    )
    return PlotlyBase.Plot(traces, layout)
end
function ADRIA.viz.rsa(
    X::DataFrame,
    y::AbstractVector{<:Real},
    factor_pairs::AbstractVector{<:NTuple{2,Symbol}};
    with_contour::Bool=true,
    title::String="Regional Sensitivity",
    kwargs...
)::PlotlyBase.Plot
    n_pairs = length(factor_pairs)
    domains, n_rows, n_cols = _grid_domains(n_pairs)
    fsz = _plotly_font_sizes(n_pairs)
    z_vals = Float64.(y)

    traces = PlotlyBase.AbstractTrace[]
    layout = PlotlyBase.Layout(;
        ADRIA_LAYOUT_DEFAULTS...,
        title=PlotlyBase.attr(; text=title, font=PlotlyBase.attr(; size=fsz.title)),
        font=PlotlyBase.attr(; family="Open Sans, sans-serif", size=fsz.label),
        width=max(600, 500 * n_cols),
        height=max(400, 450 * n_rows),
        showlegend=false
    )

    for (i, pair) in enumerate(factor_pairs)
        x_vals = Float64.(X[!, pair[1]])
        y_vals = Float64.(X[!, pair[2]])
        sfx = i == 1 ? "" : string(i)
        xd, yd, _, _ = domains[i]

        if with_contour
            push!(
                traces,
                PlotlyBase.histogram2dcontour(;
                    x=x_vals,
                    y=y_vals,
                    z=z_vals,
                    histfunc="avg",
                    colorscale="Viridis",
                    showscale=false,
                    opacity=0.4,
                    ncontours=15,
                    contours=PlotlyBase.attr(; coloring="fill", showlines=false),
                    xaxis="x$(sfx)",
                    yaxis="y$(sfx)",
                    type="histogram2dcontour"
                )
            )
        end

        push!(
            traces,
            PlotlyBase.scatter(;
                x=x_vals,
                y=y_vals,
                mode="markers",
                marker=PlotlyBase.attr(;
                    color=z_vals,
                    colorscale="Viridis",
                    showscale=(i == n_pairs),
                    opacity=0.6,
                    size=5,
                    colorbar=PlotlyBase.attr(; title="outcome")
                ),
                showlegend=false,
                xaxis="x$(sfx)",
                yaxis="y$(sfx)",
                type="scatter"
            )
        )

        layout[Symbol("xaxis$(sfx)")] = PlotlyBase.attr(;
            title_text=string(pair[1]), anchor="y$(sfx)", domain=xd,
            tickfont=PlotlyBase.attr(; size=fsz.tick)
        )
        layout[Symbol("yaxis$(sfx)")] = PlotlyBase.attr(;
            title_text=string(pair[2]), anchor="x$(sfx)", domain=yd,
            tickfont=PlotlyBase.attr(; size=fsz.tick)
        )
    end

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

# ──────────────────────────────────────────────────────────────────────────────
# outcome_map(X, y, factor/factors) — scatter of factor value vs outcome
# ──────────────────────────────────────────────────────────────────────────────

"""
    ADRIA.viz.outcome_map(X, y, factor; with_density=true, title, kwargs...) -> PlotlyBase.Plot
    ADRIA.viz.outcome_map(X, y, factors; with_density=true, title, kwargs...) -> PlotlyBase.Plot

Per-scenario scatter of factor value vs outcome, with optional 2D density contour overlay.

The multi-factor overload produces one subplot per factor using make_subplots.

# Arguments
- `X`            : Feature matrix (DataFrame, columns = factors)
- `y`            : Outcome values per scenario
- `factor`       : Symbol naming the factor (single-factor form)
- `factors`      : Vector of Symbols (multi-factor form)
- `with_density` : Add a histogram2dcontour density fill.
"""
function ADRIA.viz.outcome_map(
    X::DataFrame,
    y::AbstractVector{<:Real},
    factor::Symbol;
    with_density::Bool=true,
    title::String="Outcome Map",
    kwargs...
)::PlotlyBase.Plot
    x_f = Float64.(X[!, factor])
    y_f = Float64.(y)

    fsz = _plotly_font_sizes(1)

    traces = PlotlyBase.AbstractTrace[]

    if with_density && length(unique(x_f)) > 3
        push!(
            traces,
            PlotlyBase.histogram2dcontour(;
                x=x_f,
                y=y_f,
                colorscale="Blues",
                showscale=false,
                opacity=0.5,
                ncontours=12,
                contours=PlotlyBase.attr(; coloring="fill", showlines=false),
                type="histogram2dcontour"
            )
        )
    end

    push!(
        traces,
        PlotlyBase.scatter(;
            x=x_f,
            y=y_f,
            mode="markers",
            marker=PlotlyBase.attr(; color="#4477aa", opacity=0.4, size=5),
            type="scatter"
        )
    )

    layout = PlotlyBase.Layout(;
        ADRIA_LAYOUT_DEFAULTS...,
        title=PlotlyBase.attr(; text=title, font=PlotlyBase.attr(; size=fsz.title)),
        font=PlotlyBase.attr(; family="Open Sans, sans-serif", size=fsz.label),
        xaxis=PlotlyBase.attr(;
            title_text=string(factor), tickfont=PlotlyBase.attr(; size=fsz.tick)
        ),
        yaxis=PlotlyBase.attr(;
            title_text="outcome", tickfont=PlotlyBase.attr(; size=fsz.tick)
        )
    )
    return PlotlyBase.Plot(traces, layout)
end
function ADRIA.viz.outcome_map(
    X::DataFrame,
    y::AbstractVector{<:Real},
    factors::AbstractVector{Symbol};
    with_density::Bool=true,
    title::String="Outcome Map",
    kwargs...
)::PlotlyBase.Plot
    n_factors = length(factors)
    domains, n_rows, n_cols = _grid_domains(n_factors)
    fsz = _plotly_font_sizes(n_factors)

    traces = PlotlyBase.AbstractTrace[]
    layout = PlotlyBase.Layout(;
        ADRIA_LAYOUT_DEFAULTS...,
        title=PlotlyBase.attr(; text=title, font=PlotlyBase.attr(; size=fsz.title)),
        font=PlotlyBase.attr(; family="Open Sans, sans-serif", size=fsz.label),
        width=max(600, 380 * n_cols),
        height=max(400, 350 * n_rows),
        showlegend=false
    )

    for (i, factor) in enumerate(factors)
        x_f = Float64.(X[!, factor])
        sfx = i == 1 ? "" : string(i)
        xd, yd, _, _ = domains[i]

        if with_density && length(unique(x_f)) > 3
            push!(
                traces,
                PlotlyBase.histogram2dcontour(;
                    x=x_f,
                    y=Float64.(y),
                    colorscale="Blues",
                    showscale=false,
                    opacity=0.5,
                    ncontours=12,
                    contours=PlotlyBase.attr(; coloring="fill", showlines=false),
                    xaxis="x$(sfx)",
                    yaxis="y$(sfx)",
                    type="histogram2dcontour"
                )
            )
        end
        push!(
            traces,
            PlotlyBase.scatter(;
                x=x_f,
                y=Float64.(y),
                mode="markers",
                marker=PlotlyBase.attr(; color="#4477aa", opacity=0.4, size=5),
                showlegend=false,
                xaxis="x$(sfx)",
                yaxis="y$(sfx)",
                type="scatter"
            )
        )
        layout[Symbol("xaxis$(sfx)")] = PlotlyBase.attr(;
            title_text=string(factor), anchor="y$(sfx)", domain=xd,
            tickfont=PlotlyBase.attr(; size=fsz.tick)
        )
        layout[Symbol("yaxis$(sfx)")] = PlotlyBase.attr(;
            title_text="outcome", anchor="x$(sfx)", domain=yd,
            tickfont=PlotlyBase.attr(; size=fsz.tick)
        )
    end

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
