using Statistics
using Printf
using ADRIA: _is_discrete_factor

"""
    _get_guided_labels()::Vector{String}

Returns labels for categories of the `guided` factor.
"""
function _get_guided_labels()::Vector{String}
    return [
        "cf",
        "unguided",
        last.(split.(string.(ADRIA.decision.mcda_methods()), "."))...
    ]
end

"""
    ADRIA.viz.pawn(Si::YAXArray; opts::Dict{Symbol,<:Any}=Dict{Symbol,Any}(), fig_opts::Dict{Symbol,<:Any}=Dict{Symbol,Any}(), axis_opts::Dict{Symbol,<:Any}=Dict{Symbol,Any}())
    ADRIA.viz.pawn!(f::Union{GridLayout,GridPosition}, Si::YAXArray; opts::Dict{Symbol,<:Any}=Dict{Symbol,Any}(), axis_opts::Dict{Symbol,<:Any}=Dict{Symbol,Any}())

Display heatmap of sensitivity analysis.

# Arguments
- `Si` : Sensitivity analysis results from `pawn()`
- `opts` : Additional figure customization options
    - `normalize` : Normalize each column ∈ [0, 1] to obtain relative sensitivity
    - `factors` : List of factors to display (factors are filtered after normalization)
- `axis_opts` : Additional options to pass to adjust Axis attributes
  See: https://docs.makie.org/v0.19/api/index.html#Axis

# Arguments
- `Si` : Results from sensitivity analysis
- `opts` : Additional figure customization options
    - `normalize` : Normalize each column ∈ [0, 1] to obtain relative sensitivity
    - `factors` : List of factors to display (factors are filtered after normalization)
    - `by` : Symbol of index to sort by in descending order (defaults to `:median`)
- `fig_opts` : Additional options to pass to adjust Figure creation
See: https://docs.makie.org/v0.19/api/index.html#Figure
- `axis_opts` : Additional options to pass to adjust Axis attributes
See: https://docs.makie.org/v0.19/api/index.html#Axis

# Returns
Makie figure
"""
function ADRIA.viz.pawn!(
    g::Union{GridLayout,GridPosition},
    Si::YAXArray;
    opts::OPT_TYPE=DEFAULT_OPT_TYPE(),
    axis_opts::OPT_TYPE=DEFAULT_OPT_TYPE()
)
    set_typography_defaults!(axis_opts)
    xtick_rot = get(axis_opts, :xticklabelrotation, 2.0 / π)

    norm = get(opts, :normalize, true)
    if norm
        Si = YAXArray(Si.axes, col_normalize(Si.data))
    end

    foi = get(opts, :factors, :all)
    if foi != :all
        Si = Si[factors = At(foi)]
    end

    # Sort by
    sort_by = get(opts, :by, :median)
    Si = Si[sortperm(collect(Si[Si = At(sort_by)]); rev=true), :]

    y, x = Si.axes
    ax = Axis(
        g[1, 1];
        xticks=(1:length(x), string.(x)),
        yticks=(1:length(y), string.(y)),
        xticklabelrotation=xtick_rot,
        axis_opts...
    )
    ax.yreversed = true

    heatmap!(ax, Matrix(Si'))

    return g
end
function ADRIA.viz.pawn(
    Si::YAXArray;
    opts::OPT_TYPE=DEFAULT_OPT_TYPE(),
    fig_opts::OPT_TYPE=DEFAULT_OPT_TYPE(),
    axis_opts::OPT_TYPE=DEFAULT_OPT_TYPE()
)
    set_figure_defaults(fig_opts)
    f = Figure(; fig_opts...)
    g = f[1, 1] = GridLayout()
    ADRIA.viz.pawn!(g, Si; opts, axis_opts)

    return f
end

"""
    ADRIA.viz.tsa(rs::ResultSet, si::YAXArray; opts::Dict{Symbol,<:Any}=Dict{Symbol,Any}(), fig_opts::Dict=Dict(), axis_opts::Dict=Dict())
    ADRIA.viz.tsa!(f::Union{GridLayout,GridPosition}, rs::ResultSet, si::YAXArray; opts::Dict{Symbol,<:Any}=Dict{Symbol,Any}(), axis_opts::Dict{Symbol,<:Any}=Dict{Symbol,Any}())

Display temporal sensitivity analysis

# Arguments
- `rs` : ResultSet
- `si` : Sensitivity indices from ADRIA temporal sensitivity analysis
- `opts` : Additional figure customization options
    - `stat` : Sensitivity statistic to display (defaults to :median)
- `fig_opts` : Additional options to pass to adjust Figure creation
  See: https://docs.makie.org/v0.19/api/index.html#Figure
- `axis_opts` : Additional options to pass to adjust Axis attributes
  See: https://docs.makie.org/v0.19/api/index.html#Axis

# Returns
Makie figure
"""
function ADRIA.viz.tsa!(
    g::Union{GridLayout,GridPosition},
    rs::ResultSet,
    si::YAXArray;
    opts::OPT_TYPE=DEFAULT_OPT_TYPE(),
    axis_opts::OPT_TYPE=DEFAULT_OPT_TYPE()
)
    set_typography_defaults!(axis_opts)
    stat = get(opts, :stat, :median)

    xtick_rot = get(axis_opts, :xticklabelrotation, 2.0 / π)
    xlabel = get(axis_opts, :xlabel, "Years")
    ylabel = get(axis_opts, :ylabel, L"\text{PAWN}_\text{%$(stat)}")

    factors, Si, timesteps = si.axes
    x_tickpos, x_ticklabel = _time_labels(timesteps)
    ax = Axis(
        g[1, 1];
        xticks=(x_tickpos, x_ticklabel),
        xticklabelrotation=xtick_rot,
        xlabel=xlabel,
        ylabel=ylabel,
        axis_opts...
    )

    # Some factors in the model_spec related to CB_CALIB_PARAMS are not inputs to the model
    # so they need to be filtered before selecting `:component` col
    all_comps = model_spec(rs, collect(si.factors.val))[:, :component]

    # Hacky special case handling for SSP/RCP
    if :RCP in factors || :SSP in factors
        pushfirst!(all_comps, "EnvironmentalLayer")
    end

    # min_step = (1 / 0.05)
    # color_weight = min((1.0 / (length(factors) / min_step)), 0.6)
    comps = unique(all_comps)
    dc = distinguishable_colors(length(comps), [RGB(1, 1, 1), RGB(0, 0, 0)]; dropseed=true)
    lns = Plot[
        series!(
            ax,
            si[Si = At(stat)][findall(all_comps .== _cmp), :].data;
            labels=repeat([_cmp], count(all_comps .== _cmp)),
            solid_color=(dc[i], 0.2)
        )
        for (i, _cmp) in enumerate(comps)
    ]

    Legend(g[1, 2], lns, comps)

    return g
end
function ADRIA.viz.tsa(
    rs::ResultSet,
    si::YAXArray;
    opts::OPT_TYPE=DEFAULT_OPT_TYPE(),
    fig_opts::OPT_TYPE=DEFAULT_OPT_TYPE(),
    axis_opts::OPT_TYPE=DEFAULT_OPT_TYPE()
)
    fig_opts[:size] = get(fig_opts, :size, (1200, 600))
    set_figure_defaults(fig_opts)
    f = Figure(; fig_opts...)
    g = f[1, 1] = GridLayout()
    ADRIA.viz.tsa!(g, rs, si; opts, axis_opts)

    return f
end

"""
    ADRIA.viz.rsa!(g, X, y, foi; with_contour=true, opts=..., axis_opts=...)
    ADRIA.viz.rsa(X, y, foi; with_contour=true, opts=..., fig_opts=..., axis_opts=...) -> Figure

2D scatter of two input factors colored by outcome metric.

# Arguments
- `g`            : GridLayout or GridPosition to draw into
- `X`            : Feature matrix (DataFrame, columns = factors)
- `y`            : Outcome values per scenario
- `foi`          : NTuple{2,Symbol} — (x-axis factor, y-axis factor)
- `with_contour` : Add a tricontourf overlay. Exploratory only -- Delaunay-based
                   interpolation can produce artefacts in sparse or irregularly sampled
                   scenario spaces and is meaningless for discrete factors.
"""
function ADRIA.viz.rsa!(
    g::Union{GridLayout,GridPosition},
    X::DataFrame,
    y::AbstractVector{<:Real},
    foi::NTuple{2,Symbol};
    with_contour::Bool=true,
    _add_colorbar::Bool=true,
    opts::OPT_TYPE=DEFAULT_OPT_TYPE(),
    axis_opts::OPT_TYPE=DEFAULT_OPT_TYPE()
)
    set_typography_defaults!(axis_opts)
    xlabel = get(axis_opts, :xlabel, string(foi[1]))
    ylabel = get(axis_opts, :ylabel, string(foi[2]))

    ax = Axis(g[1, 1]; xlabel=xlabel, ylabel=ylabel, axis_opts...)

    x_vals = Float64.(X[!, foi[1]])
    y_vals = Float64.(X[!, foi[2]])

    if with_contour
        tricontourf!(ax, x_vals, y_vals, y; colormap=:viridis, alpha=0.3)
    end

    scatter!(
        ax, x_vals, y_vals;
        color=y, colormap=:viridis,
        strokewidth=0, alpha=0.5, markersize=10
    )

    if _add_colorbar
        Colorbar(g[1, 2]; colormap=:viridis, label="outcome")
    end

    return g
end
function ADRIA.viz.rsa(
    X::DataFrame,
    y::AbstractVector{<:Real},
    foi::NTuple{2,Symbol};
    with_contour::Bool=true,
    opts::OPT_TYPE=DEFAULT_OPT_TYPE(),
    fig_opts::OPT_TYPE=DEFAULT_OPT_TYPE(),
    axis_opts::OPT_TYPE=DEFAULT_OPT_TYPE()
)::Figure
    get!(fig_opts, :size, (700, 500))
    set_figure_defaults(fig_opts)
    f = Figure(; fig_opts...)
    g = f[1, 1] = GridLayout()
    ADRIA.viz.rsa!(g, X, y, foi; with_contour=with_contour, opts=opts, axis_opts=axis_opts)
    return f
end
function ADRIA.viz.rsa!(
    g::Union{GridLayout,GridPosition},
    X::DataFrame,
    y::AbstractVector{<:Real},
    factor_pairs::AbstractVector{<:NTuple{2,Symbol}};
    with_contour::Bool=true,
    opts::OPT_TYPE=DEFAULT_OPT_TYPE(),
    axis_opts::OPT_TYPE=DEFAULT_OPT_TYPE()
)
    n_pairs = length(factor_pairs)
    n_rows, n_cols = _calc_gridsize(n_pairs)
    for (idx, pair) in enumerate(factor_pairs)
        row = div(idx - 1, n_cols) + 1
        col = mod(idx - 1, n_cols) + 1
        ax_opts = copy(axis_opts)
        set_typography_defaults!(ax_opts; n_panels=n_pairs)
        ADRIA.viz.rsa!(
            g[row, col], X, y, pair;
            with_contour=with_contour, _add_colorbar=false, opts=opts, axis_opts=ax_opts
        )
    end
    Colorbar(g[1:n_rows, n_cols + 1]; colormap=:viridis, label="outcome")
    rg, cg = _adaptive_gap(n_pairs)
    if g isa GridLayout
        rowgap!(g, rg)
        colgap!(g, cg)
    end
    return g
end
function ADRIA.viz.rsa(
    X::DataFrame,
    y::AbstractVector{<:Real},
    factor_pairs::AbstractVector{<:NTuple{2,Symbol}};
    with_contour::Bool=true,
    opts::OPT_TYPE=DEFAULT_OPT_TYPE(),
    fig_opts::OPT_TYPE=DEFAULT_OPT_TYPE(),
    axis_opts::OPT_TYPE=DEFAULT_OPT_TYPE()
)::Figure
    n_pairs = length(factor_pairs)
    n_rows, n_cols = _calc_gridsize(n_pairs)
    set_figure_defaults(fig_opts; n_rows=n_rows, n_cols=n_cols)
    f = Figure(; fig_opts...)
    g = f[1, 1] = GridLayout()
    ADRIA.viz.rsa!(
        g, X, y, factor_pairs; with_contour=with_contour, opts=opts, axis_opts=axis_opts
    )
    return f
end
function ADRIA.viz.rsa!(
    g::Union{GridLayout,GridPosition},
    X::DataFrame,
    y::AbstractVector{<:Real},
    factors::AbstractVector{Symbol};
    with_density::Bool=true,
    opts::OPT_TYPE=DEFAULT_OPT_TYPE(),
    axis_opts::OPT_TYPE=DEFAULT_OPT_TYPE()
)
    n_factors = length(factors)
    n_rows, n_cols = _calc_gridsize(n_factors)
    y_f = Float64.(y)
    for (idx, factor) in enumerate(factors)
        row = div(idx - 1, n_cols) + 1
        col = mod(idx - 1, n_cols) + 1
        ax_opts = copy(axis_opts)
        set_typography_defaults!(ax_opts; n_panels=n_factors)
        ax_opts[:xlabel] = string(factor)
        ax_opts[:ylabel] = get(axis_opts, :ylabel, "outcome")
        ax = Axis(g[row, col]; ax_opts...)
        x_raw = X[!, factor]
        x_f = Float64.(x_raw)
        is_discrete = eltype(x_raw) <: Integer || length(unique(x_f)) <= 10
        if with_density && !is_discrete
            hexbin!(ax, x_f, y_f; bins=20, colormap=:Blues, threshold=1, strokewidth=0)
        end
        scatter!(
            ax, x_f, y_f;
            color=y_f, colormap=:viridis, strokewidth=0, alpha=0.5, markersize=8
        )
    end
    Colorbar(g[1:n_rows, n_cols + 1]; colormap=:viridis, label="outcome")
    rg, cg = _adaptive_gap(n_factors)
    if g isa GridLayout
        rowgap!(g, rg)
        colgap!(g, cg)
    end
    return g
end
function ADRIA.viz.rsa(
    X::DataFrame,
    y::AbstractVector{<:Real},
    factors::AbstractVector{Symbol};
    with_density::Bool=true,
    opts::OPT_TYPE=DEFAULT_OPT_TYPE(),
    fig_opts::OPT_TYPE=DEFAULT_OPT_TYPE(),
    axis_opts::OPT_TYPE=DEFAULT_OPT_TYPE()
)::Figure
    n_rows, n_cols = _calc_gridsize(length(factors))
    set_figure_defaults(fig_opts; n_rows=n_rows, n_cols=n_cols)
    f = Figure(; fig_opts...)
    g = f[1, 1] = GridLayout()
    ADRIA.viz.rsa!(
        g, X, y, factors; with_density=with_density, opts=opts, axis_opts=axis_opts
    )
    return f
end

"""
    ADRIA.viz.outcome_map!(g, X, y, factor; with_density=true, opts=..., axis_opts=...)
    ADRIA.viz.outcome_map(X, y, factor; ...) -> Figure
    ADRIA.viz.outcome_map!(g, X, y, factors; with_density=true, opts=..., axis_opts=...)
    ADRIA.viz.outcome_map(X, y, factors; ...) -> Figure

Per-scenario scatter of factor value vs outcome, colored by outcome value.

Points are colored by their outcome so that outcome quality remains visible even in
densely populated regions of the factor space (e.g. when no-intervention scenarios
dominate). Alpha transparency of overlapping points provides density context.

The multi-factor overloads produce one subplot per factor (one row per factor, single
column). All kwargs are forwarded unchanged to each single-factor call.

- `with_density` : Unused; kept for API compatibility.
"""
function ADRIA.viz.outcome_map!(
    g::Union{GridLayout,GridPosition},
    X::DataFrame,
    y::AbstractVector{<:Real},
    factor::Symbol;
    with_density::Bool=true,
    _add_colorbar::Bool=true,
    opts::OPT_TYPE=DEFAULT_OPT_TYPE(),
    axis_opts::OPT_TYPE=DEFAULT_OPT_TYPE()
)
    set_typography_defaults!(axis_opts)
    xlabel = get(axis_opts, :xlabel, string(factor))
    ylabel = get(axis_opts, :ylabel, "outcome")

    x_raw = X[!, factor]
    x_f = Float64.(x_raw)
    y_f = Float64.(y)

    is_discrete = eltype(x_raw) <: Integer || length(unique(x_f)) <= 10
    if is_discrete && length(unique(x_f)) > 1
        span = maximum(x_f) - minimum(x_f)
        jitter_scale = span / max(1, length(unique(x_f))) * 0.3
        x_plot = x_f .+ jitter_scale .* (rand(length(x_f)) .- 0.5)
    else
        x_plot = x_f
    end

    ax = Axis(g[1, 1]; xlabel=xlabel, ylabel=ylabel, axis_opts...)

    scatter!(
        ax,
        x_plot,
        y_f;
        color=y_f,
        colormap=:viridis,
        alpha=0.5,
        markersize=6,
        strokewidth=0
    )
    if _add_colorbar
        Colorbar(g[1, 2]; colormap=:viridis, colorrange=extrema(y_f), label="outcome")
    end

    if factor == :guided
        fv_labels = _get_guided_labels()
        unique_vals = sort(unique(x_f))
        n_labels = min(length(fv_labels), length(unique_vals))
        ax.xticks = (unique_vals[1:n_labels], fv_labels[1:n_labels])
        ax.xticklabelrotation = pi / 4
    end

    return g
end
function ADRIA.viz.outcome_map(
    X::DataFrame,
    y::AbstractVector{<:Real},
    factor::Symbol;
    with_density::Bool=true,
    opts::OPT_TYPE=DEFAULT_OPT_TYPE(),
    fig_opts::OPT_TYPE=DEFAULT_OPT_TYPE(),
    axis_opts::OPT_TYPE=DEFAULT_OPT_TYPE()
)::Figure
    set_figure_defaults(fig_opts)
    f = Figure(; fig_opts...)
    g = f[1, 1] = GridLayout()
    ADRIA.viz.outcome_map!(
        g, X, y, factor; with_density=with_density, opts=opts, axis_opts=axis_opts
    )
    return f
end
function ADRIA.viz.outcome_map!(
    g::Union{GridLayout,GridPosition},
    X::DataFrame,
    y::AbstractVector{<:Real},
    factors::AbstractVector{Symbol};
    with_density::Bool=true,
    opts::OPT_TYPE=DEFAULT_OPT_TYPE(),
    axis_opts::OPT_TYPE=DEFAULT_OPT_TYPE()
)
    n_factors = length(factors)
    n_rows, n_cols = _calc_gridsize(n_factors)
    for (i, factor) in enumerate(factors)
        row = div(i - 1, n_cols) + 1
        col = mod(i - 1, n_cols) + 1
        ADRIA.viz.outcome_map!(
            g[row, col], X, y, factor;
            with_density=with_density, _add_colorbar=false, opts=opts,
            axis_opts=axis_opts
        )
    end
    y_f = Float64.(y)
    Colorbar(
        g[1:n_rows, n_cols + 1]; colormap=:viridis, colorrange=extrema(y_f), label="outcome"
    )
    return g
end
function ADRIA.viz.outcome_map(
    X::DataFrame,
    y::AbstractVector{<:Real},
    factors::AbstractVector{Symbol};
    with_density::Bool=true,
    opts::OPT_TYPE=DEFAULT_OPT_TYPE(),
    fig_opts::OPT_TYPE=DEFAULT_OPT_TYPE(),
    axis_opts::OPT_TYPE=DEFAULT_OPT_TYPE()
)::Figure
    n_rows, n_cols = _calc_gridsize(length(factors))
    panel_w, panel_h, gap = 400, 350, 20
    auto_w = n_cols * panel_w + (n_cols - 1) * gap
    auto_h = n_rows * panel_h + (n_rows - 1) * gap
    get!(fig_opts, :size, (auto_w, auto_h))
    set_figure_defaults(fig_opts)
    f = Figure(; fig_opts...)
    g = f[1, 1] = GridLayout()
    ADRIA.viz.outcome_map!(
        g, X, y, factors; with_density=with_density, opts=opts, axis_opts=axis_opts
    )
    return f
end

"""
    _series_convergence(g::GridPosition, Si_conv::YAXArray, factors::Vector{Symbol};
        opts::Dict{Symbol,<:Any}=Dict{Symbol,Any}(:plot_overlay => true), axis_opts::Dict{Symbol,<:Any}=Dict{Symbol,Any}())

Plot sensitivity values for an increasing number of scenarios as a series, with each member
of the series representing a factor or model component.

# Arguments
- `Si_conv` : Produced using ADRIA.analysis.convergence()
- `factors` : Factors/model components to plot
- `opts` : Additional figure customization options
        - `plot_overlay` : true, to plot overlaid series (the default), false to plot series as grid of subplots
- `axis_opts` : Additional options to pass to adjust Axis attributes
      See: https://docs.makie.org/v0.19/api/index.html#Axis

# Returns
Makie figure
"""
function _series_convergence(
    g::GridPosition,
    Si_conv::YAXArray,
    factors::Vector{Symbol};
    opts::OPT_TYPE=DEFAULT_OPT_TYPE(:plot_overlay => true),
    axis_opts::OPT_TYPE=DEFAULT_OPT_TYPE()
)
    plot_overlay = get(opts, :plot_overlay, true)
    n_scenarios = collect(lookup(Si_conv, :n_scenarios))
    grps = Dict(Symbol(foi_grp) => foi_grp .== factors for foi_grp in factors)

    xlabel = pop!(axis_opts, :xlabel, "N scenarios")
    ylabel = pop!(axis_opts, :ylabel, "Sensitivity")

    if :title in keys(axis_opts)
        title_val = pop!(axis_opts, :title)
    end

    _colors::Dict{Symbol,COLOR_TYPE} = colors(grps)
    if plot_overlay
        ax = Axis(g; axis_opts...)
        scenarios_confint!(
            ax,
            permutedims(Si_conv[Si = At([:lb, :median, :ub])], (3, 1, 2)).data,
            collect(keys(grps)),
            _colors;
            x_vals=n_scenarios
        )
        ax.xlabel = xlabel
        ax.ylabel = ylabel

        if @isdefined(title_val)
            ax.title = title_val
        end
    else
        _alphas::Dict{Symbol,Float64} = alphas(grps)

        n_factors::Int64 = length(factors)
        if n_factors > 30
            ArgumentError("Too many factors to plot. Maximum number supported is 30.")
        end

        n_rows, n_cols = _calc_gridsize(n_factors)
        factor_names = String.(factors)
        axs = Axis[]
        step::Int64 = 1

        for row = 1:n_rows, col = 1:n_cols
            ax::Axis = Axis(g[row, col]; title=factor_names[step], axis_opts...)

            lines!(
                ax,
                n_scenarios,
                Si_conv[Si = At(:median)][factors = At(factors[step])].data;
                color=(_colors[factors[step]], _alphas[factors[step]])
            )

            band!(
                ax,
                n_scenarios,
                Si_conv[Si = At(:lb), factors = At(factors[step])].data,
                Si_conv[Si = At(:ub), factors = At(factors[step])].data;
                color=(_colors[factors[step]], _alphas[factors[step]])
            )
            step += 1
            push!(axs, ax)
            if step > n_factors
                break
            end
        end

        if n_factors > 1
            linkyaxes!(axs...)
            x_lbl_sz = axis_opts[:xlabelsize]
            y_lbl_sz = get(axis_opts, :ylabelsize, x_lbl_sz)
            Label(g[n_rows + 1, :]; text=xlabel, fontsize=x_lbl_sz + 2)
            Label(g[:, 0]; text=ylabel, fontsize=y_lbl_sz + 2, rotation=pi / 2)

            if @isdefined(title_val)
                Label(g[0, :]; text=title_val, fontsize=axis_opts[:titlesize])
            end
        else
            axs[1].xlabel = xlabel
            axs[1].ylabel = ylabel

            if @isdefined(title_val)
                axs[1].title = title_val
            end
        end

        try
            # Clear empty figures
            trim!(g)
        catch err
            if !(err isa MethodError)
                # GridPosition plots a single figure so does
                # not need empty figures to be cleared
                # If any other error is encountered, something else happened.
                rethrow(err)
            end
        end
    end

    return g
end

"""
    _heatmap_convergence(g::GridPosition, Si_conv::YAXArray, factors::Vector{Symbol};
        opts::Dict{Symbol,<:Any}=Dict{Symbol,Any}(), axis_opts::Dict{Symbol,<:Any}=Dict{Symbol,Any}())

Plot sensitivity values for an increasing number of scenarios as a heatmap, with each row
of the heatmap representing a factor or model component.

# Arguments
- `Si_conv` : Produced using ADRIA.analysis.convergence()
- `factors` : Factors/model components to plot
- `opts` : Additional figure customization options
    - `colorbar_label` : string indicating how to label colorbar in heatmap
    - `color_map` : colormap to use for heatmap
- `axis_opts` : Additional options to pass to adjust Axis attributes
      See: https://docs.makie.org/v0.19/api/index.html#Axis

# Returns
Makie figure
"""
function _heatmap_convergence(
    g::GridPosition,
    Si_conv::YAXArray,
    factors::Vector{Symbol};
    opts::OPT_TYPE=DEFAULT_OPT_TYPE(),
    axis_opts::OPT_TYPE=DEFAULT_OPT_TYPE()
)
    set_typography_defaults!(axis_opts)
    y_label = get(axis_opts, :ylabel, "Factors")
    x_label = get(axis_opts, :xlabel, "N scenarios")

    z = Array(Si_conv[Si = At(:median)])
    xtick_vals = (1:length(Si_conv.n_scenarios), string.(Si_conv.n_scenarios))
    ytick_vals = (1:length(factors), string.(factors))

    ax = Axis(
        g[1, 1];
        xticks=xtick_vals,
        yticks=ytick_vals,
        xlabel=x_label,
        ylabel=y_label,
        axis_opts...
    )
    heatmap!(ax, z')
    colorbar_label = get(opts, :colorbar_label, "Relative Sensitivity")
    color_map = get(opts, :color_map, :viridis)

    Colorbar(g[1, 2]; colormap=color_map, label=colorbar_label, height=Relative(0.65))
    return g
end

"""
    ADRIA.viz.convergence(Si_conv::YAXArray, factors::Vector{Symbol}; series_opts::Dict{Symbol,<:Any}=Dict{Symbol,Any}(),
        axis_opts::Dict{Symbol,<:Any}=Dict{Symbol,Any}())
    ADRIA.viz.convergence!(f::Figure, Si_conv::YAXArray, factors::Vector{Symbol}; series_opts::Dict{Symbol,<:Any}=Dict{Symbol,Any}(),
        axis_opts::Dict{Symbol,<:Any}=Dict{Symbol,Any}())

Plot sensitivity metric for increasing number of scenarios to illustrate convergence.

# Arguments
- `Si_conv` : Produced using ADRIA.analysis.convergence()
- `factors` : Factors/model components to plot
- `opts` : Additional figure customization options
    - `viz_type` : :heat_map to plot heatmap, :series to plot as series
    - `plot_overlay` : true, to plot overlaid series (the default), false to plot series as grid of subplots
    - `colorbar_label` : string indicating how to label colorbar in heatmap
    - `color_map` : colormap to use for heatmap
- `fig_opts` : Additional options to pass to adjust Figure creation
  See: https://docs.makie.org/v0.19/api/index.html#Figure
- `axis_opts` : Additional options to pass to adjust Axis attributes
  See: https://docs.makie.org/v0.19/api/index.html#Axis

# Returns
Makie figure
"""
function ADRIA.viz.convergence(
    Si_conv::YAXArray,
    factors::Vector{Symbol};
    opts::OPT_TYPE=DEFAULT_OPT_TYPE(:viz_type => :series),
    fig_opts::OPT_TYPE=DEFAULT_OPT_TYPE(),
    axis_opts::OPT_TYPE=DEFAULT_OPT_TYPE()
)
    fig_opts[:size] = get(fig_opts, :size, (1200, 600))
    set_figure_defaults(fig_opts)
    f = Figure(; fig_opts...)
    ADRIA.viz.convergence!(
        f[1, 1],
        Si_conv,
        factors;
        opts=opts,
        axis_opts=axis_opts
    )
    return f
end
function ADRIA.viz.convergence!(
    g::Union{GridLayout,GridPosition},
    Si_conv::YAXArray,
    factors::Vector{Symbol};
    opts::OPT_TYPE=DEFAULT_OPT_TYPE(:viz_type => :series),
    axis_opts::OPT_TYPE=DEFAULT_OPT_TYPE()
)
    _n =
        (get(opts, :viz_type, :series) == :series && !get(opts, :plot_overlay, true)) ?
        length(factors) : 1
    set_typography_defaults!(axis_opts; n_panels=_n)
    viz_type = get(opts, :viz_type, :series)
    if viz_type == :series
        return _series_convergence(g, Si_conv, factors; opts=opts, axis_opts=axis_opts)
    elseif viz_type == :heatmap
        return _heatmap_convergence(g, Si_conv, factors; opts=opts, axis_opts=axis_opts)
    else
        error("Convergence plot $(viz_type) is not expected.")
    end
end
