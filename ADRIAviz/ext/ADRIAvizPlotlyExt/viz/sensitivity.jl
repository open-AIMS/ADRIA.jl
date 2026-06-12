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
        Si = Si[factors=At(factors)]
    end

    # Sort rows by chosen statistic (descending)
    sort_vals = collect(Si[Si=At(by)])
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

    layout = PlotlyBase.Layout(;
        ADRIA_LAYOUT_DEFAULTS...,
        title_text=title,
        xaxis=PlotlyBase.attr(; title_text="SI Statistic"),
        yaxis=PlotlyBase.attr(;
            title_text="Factor", ticktext=factor_names, tickvals=factor_names
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

    layout = PlotlyBase.Layout(;
        ADRIA_LAYOUT_DEFAULTS...,
        title_text=title,
        xaxis=PlotlyBase.attr(; title_text=xlabel, domain=_subplot_domain(1, n_factors)),
        yaxis=PlotlyBase.attr(; title_text=ylabel)
    )
    for i in 2:n_factors
        fname = factor_names[i]
        axis_sfx = string(i)
        layout[Symbol("xaxis$(axis_sfx)")] = PlotlyBase.attr(;
            title_text="$(fname) value", domain=_subplot_domain(i, n_factors)
        )
        layout[Symbol("yaxis$(axis_sfx)")] = PlotlyBase.attr(; title_text=ylabel)
    end

    return PlotlyBase.Plot(traces, layout)
end

# helper: compute subplot x-domain for panel i of n total
function _subplot_domain(i::Int, n::Int)
    gap = 0.05
    w = (1.0 - gap * (n - 1)) / n
    lo = (i - 1) * (w + gap)
    return [lo, lo + w]
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

    layout = PlotlyBase.Layout(;
        ADRIA_LAYOUT_DEFAULTS...,
        title_text=title,
        xaxis=PlotlyBase.attr(; title_text=xlabel, domain=_subplot_domain(1, n_factors)),
        yaxis=PlotlyBase.attr(; title_text=ylabel)
    )
    for i in 2:n_factors
        fname = factor_names[i]
        axis_sfx = string(i)
        layout[Symbol("xaxis$(axis_sfx)")] = PlotlyBase.attr(;
            title_text=string(fname), domain=_subplot_domain(i, n_factors)
        )
        layout[Symbol("yaxis$(axis_sfx)")] = PlotlyBase.attr(; title_text=ylabel)
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

    layout = PlotlyBase.Layout(;
        ADRIA_LAYOUT_DEFAULTS...,
        title_text=title,
        xaxis=PlotlyBase.attr(; title_text=xlabel),
        yaxis=PlotlyBase.attr(; title_text=ylabel)
    )
    return PlotlyBase.Plot(traces, layout)
end
