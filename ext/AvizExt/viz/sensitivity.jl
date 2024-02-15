using Statistics
using Printf
using ADRIA.sensitivity: _get_cat_quantile

"""
    ADRIA.viz.pawn(Si::NamedDimsArray; opts::Dict=Dict(), fig_opts::Dict=Dict(), axis_opts::Dict=Dict())
    ADRIA.viz.pawn!(f::Union{GridLayout,GridPosition}, Si::NamedDimsArray; opts::Dict=Dict(), axis_opts::Dict=Dict())

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
    Si::NamedDimsArray;
    opts::Dict=Dict(),
    axis_opts::Dict=Dict(),
)
    xtick_rot = get(axis_opts, :xticklabelrotation, 2.0 / π)

    norm = get(opts, :normalize, true)
    if norm
        Si = col_normalize(Si)
    end

    foi = get(opts, :factors, :all)
    if foi != :all
        Si = Si(foi, :)
    end

    # Sort by
    sort_by = get(opts, :by, :median)
    Si = Si[sortperm(Si(Si=sort_by), rev=true), :]

    y, x = axiskeys(Si)
    ax = Axis(
        g[1, 1],
        xticks=(1:length(x), string.(x)),
        yticks=(1:length(y), string.(y)),
        xticklabelrotation=xtick_rot;
        axis_opts...
    )
    ax.yreversed = true

    heatmap!(ax, Matrix(Si'))

    return g
end
function ADRIA.viz.pawn(
    Si::NamedDimsArray; opts::Dict=Dict(), fig_opts::Dict=Dict(), axis_opts::Dict=Dict()
)
    f = Figure(; fig_opts...)
    g = f[1, 1] = GridLayout()
    ADRIA.viz.pawn!(g, Si; opts, axis_opts)

    return f
end

"""
    ADRIA.viz.tsa(rs::ResultSet, si::NamedDimsArray; opts::Dict=Dict(), fig_opts::Dict=Dict(), axis_opts::Dict=Dict())
    ADRIA.viz.tsa!(f::Union{GridLayout,GridPosition}, rs::ResultSet, si::NamedDimsArray; opts, axis_opts)

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
    g::Union{GridLayout,GridPosition}, rs::ResultSet, si::NamedDimsArray; opts, axis_opts
)
    stat = get(opts, :stat, :median)

    xtick_rot = get(axis_opts, :xticklabelrotation, 2.0 / π)
    xlabel = get(axis_opts, :xlabel, "Years")
    ylabel = get(axis_opts, :ylabel, L"\text{PAWN}_\text{%$(stat)}")

    factors, Si, timesteps = axiskeys(si)
    x_tickpos, x_ticklabel = _time_labels(timesteps)
    ax = Axis(
        g[1, 1],
        xticks=(x_tickpos, x_ticklabel),
        xticklabelrotation=xtick_rot,
        xlabel=xlabel,
        ylabel=ylabel;
        axis_opts...
    )

    all_comps = model_spec(rs)[:, :component]

    # Hacky special case handling for SSP/RCP
    if :RCP in factors || :SSP in factors
        pushfirst!(all_comps, "EnvironmentalLayer")
    end

    # min_step = (1 / 0.05)
    # color_weight = min((1.0 / (length(factors) / min_step)), 0.6)
    comps = unique(all_comps)
    dc = distinguishable_colors(length(comps), [RGB(1, 1, 1), RGB(0, 0, 0)], dropseed=true)
    lns = Plot[
        series!(
            ax,
            si(Si=stat)[findall(all_comps .== _cmp), :],
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
    si::NamedDimsArray;
    opts::Dict=Dict(),
    fig_opts::Dict=Dict(),
    axis_opts::Dict=Dict(),
)
    f = Figure(; fig_opts...)
    g = f[1, 1] = GridLayout()
    ADRIA.viz.tsa!(g, rs, si; opts, axis_opts)

    return f
end

"""
    ADRIA.viz.rsa(rs::ResultSet, si::Dataset, factors::Vector{String}; opts::Dict=Dict(), fig_opts::Dict=Dict(), axis_opts::Dict=Dict())
    ADRIA.viz.rsa!(f::Union{GridLayout,GridPosition}, rs::ResultSet, si::Dataset, factors::Vector{String}; opts::Dict=Dict(), axis_opts::Dict=Dict())

Plot regional sensitivities of up to 30 factors.

# Arguments
- `rs` : ResultSet
- `si` : Results from ADRIA regional sensitivity analysis
- `opts` : Additional figure customization options
- `fig_opts` : Additional options to pass to adjust Figure creation
  See: https://docs.makie.org/v0.19/api/index.html#Figure
- `axis_opts` : Additional options to pass to adjust Axis attributes
  See: https://docs.makie.org/v0.19/api/index.html#Axis

# Returns
Makie figure
"""
function ADRIA.viz.rsa!(
    g::Union{GridLayout,GridPosition},
    rs::ResultSet,
    si::Dataset,
    factors::Vector{Symbol};
    opts::Dict=Dict(),
    axis_opts::Dict=Dict(),
)
    n_factors::Int64 = length(factors)
    if n_factors > 30
        ArgumentError("Too many factors to plot. Maximum number supported is 30.")
    end

    n_rows, n_cols = _calc_gridsize(n_factors)

    xlabel = get(axis_opts, :xlabel, "Factor Value")
    ylabel = get(axis_opts, :ylabel, L"\text{Relative } S_{i}")

    if :title in keys(axis_opts)
        title_val = pop!(axis_opts, :title)
    end

    # min_step = (1 / 0.05)
    # color_weight = min((1.0 / (length(factors) / min_step)), 0.6)

    # Color by component
    ms = model_spec(rs)
    foi = ms.fieldname .∈ [factors]
    all_comps = ms[foi, :component]
    f_names = ms[foi, :fieldname]
    h_names = ms[foi, :name]
    dist_params = ms[foi, :dist_params]
    f_types = ms[foi, :ptype]

    if any(f_names .== :guided)
        fv_labels = [
            "unguided", "cf", last.(split.(string.(ADRIA.decision.mcda_methods()), "."))...
        ]
    end

    # Hacky special case handling for SSP/RCP
    if :RCP in factors || :SSP in factors
        loc = first(findall((factors .== :RCP) .|| (factors .== :SSP)))
        insert!(all_comps, loc, "EnvironmentalLayer")
        insert!(f_names, loc, "RCP")
        insert!(h_names, loc, "SSP/RCP")
        insert!(dist_params, loc, (1, length(unique(rs.inputs.RCP))))
    end

    # comps = unique(all_comps)
    # dc = distinguishable_colors(length(comps), [RGB(1, 1, 1), RGB(0, 0, 0)], dropseed=true)
    curr::Int64 = 1
    axs = Axis[]
    for r in 1:n_rows
        for c in 1:n_cols
            f_name = factors[curr]
            f_type = f_types[curr]
            f_vals = rs.inputs[:, f_name]

            if f_type == "unordered categorical"
                fv_s = _get_cat_quantile(
                    ms[ms.fieldname .== f_name, :], f_name, collect(si[f_name].axes[1])
                )
            else
                fv_s = round.(quantile(f_vals, collect(si[f_name].axes[1])), digits=2)
            end
            ax::Axis = Axis(
                g[r, c];
                title=h_names[f_names .== f_name][1],
                axis_opts...,
            )

            scatterlines!(ax, fv_s, collect(si[f_name]); markersize=15)
            if f_name == :guided
                ax.xticks = (fv_s, fv_labels)
                ax.xticklabelrotation = pi / 4
            end
            push!(axs, ax)
            curr += 1

            if curr > n_factors
                break
            end
        end
    end

    linkyaxes!(axs...)
    Label(g[end+1, :], text=xlabel, fontsize=32)
    Label(g[1:(end - 1), 0]; text=ylabel, fontsize=32, rotation=π / 2.0)

    if :title in keys(axis_opts)
        Label(g[0, :], text=title_val, fontsize=40)
    end

    # Clear empty figures
    trim!(g)

    return g
end
function ADRIA.viz.rsa(
    rs::ResultSet,
    si::Dataset,
    factors::Vector{Symbol};
    opts::Dict=Dict(),
    fig_opts::Dict=Dict(),
    axis_opts::Dict=Dict(),
)
    f = Figure(; fig_opts...)
    g = f[1, 1] = GridLayout()
    ADRIA.viz.rsa!(g, rs, si, factors; axis_opts)

    return f
end

"""
    ADRIA.viz.outcome_map(rs::ResultSet, outcomes::NamedDimsArray, factors::Vector{String}; opts::Dict=Dict(), fig_opts::Dict=Dict(), axis_opts::Dict=Dict())
    ADRIA.viz.outcome_map!(f::Union{GridLayout,GridPosition}, rs::ResultSet, outcomes::NamedDimsArray, factors::Vector{String}; opts, axis_opts)

Plot outcomes mapped to factor regions for up to 30 factors.

# Arguments
- `rs` : ResultSet
- `outcomes` : ADRIA Outcome Mapping results
- `factors` : The factors of interest to display
- `opts` : Additional figure customization options
- `fig_opts` : Additional options to pass to adjust Figure creation
  See: https://docs.makie.org/v0.19/api/index.html#Figure
- `axis_opts` : Additional options to pass to adjust Axis attributes
  See: https://docs.makie.org/v0.19/api/index.html#Axis

# Returns
Makie figure
"""
function ADRIA.viz.outcome_map!(
    g::Union{GridLayout,GridPosition},
    rs::ResultSet,
    outcomes::NamedDimsArray,
    factors::Vector{Symbol};
    opts::Dict=Dict(),
    axis_opts::Dict=Dict(),
)
    # TODO: Clean up and compartmentalize as a lot of code here are duplicates of those
    #       found in `rsa()`
    n_factors::Int64 = length(factors)
    if n_factors > 30
        ArgumentError("Too many factors to plot. Maximum number supported is 30.")
    end

    n_rows, n_cols = _calc_gridsize(n_factors)

    xlabel = pop!(axis_opts, :xlabel, "Factor Value")
    ylabel = pop!(axis_opts, :ylabel, "Outcome")

    if :title in keys(axis_opts)
        title_val = pop!(axis_opts, :title)
    end

    # min_step = (1 / 0.05)
    # color_weight = min((1.0 / (length(factors) / min_step)), 0.6)

    # Color by component
    ms = model_spec(rs)
    foi = ms.fieldname .∈ [factors]
    all_comps = ms[foi, :component]
    f_names = ms[foi, :fieldname]
    h_names = ms[foi, :name]
    dist_params = ms[foi, :dist_params]

    # Hacky special case handling for SSP/RCP
    if :RCP in factors || :SSP in factors
        loc = first(findall((factors .== :RCP) .|| (factors .== :SSP)))
        insert!(all_comps, loc, "EnvironmentalLayer")
        insert!(f_names, loc, "RCP")
        insert!(h_names, loc, "SSP/RCP")
        insert!(dist_params, loc, (1, length(unique(rs.inputs.RCP))))
    end

    bin_slices, factor_list, CIs = axiskeys(outcomes)
    b_slices = parse.(Float64, bin_slices)

    if any(f_names .== :guided)
        fv_labels = [
            "unguided", "cf", last.(split.(string.(ADRIA.decision.mcda_methods()), "."))...
        ]
    end
    curr::Int64 = 1
    axs = Axis[]
    for r in 1:n_rows
        for c in 1:n_cols
            f_name = Symbol(factors[curr])
            f_vals = rs.inputs[:, f_name]

            if f_name == :guided
                fv_s = collect(1:length(fv_labels))
            else
                fv_s = round.(quantile(f_vals, b_slices), digits=2)
            end

            ax::Axis = Axis(
                g[r, c]; title=h_names[f_names .== factors[curr]][1], axis_opts...
            )

            band!(
                ax,
                fv_s[.!ismissing.(outcomes(; factors=f_name, CI=:lower))],
                collect(skipmissing(outcomes(; factors=f_name, CI=:lower))),
                collect(skipmissing(outcomes(; factors=f_name, CI=:upper))),
            )
            scatterlines!(ax, fv_s, outcomes(; factors=f_name, CI=:mean); markersize=15)

            if f_name == :guided
                ax.xticks = (fv_s, fv_labels)
                ax.xticklabelrotation = pi / 4
            end

            push!(axs, ax)
            curr += 1

            if curr > n_factors
                break
            end
        end
    end

    if n_factors > 1
        linkyaxes!(axs...)
        Label(g[n_rows + 1, :]; text=xlabel, fontsize=24)
        Label(g[:, 0]; text=ylabel, fontsize=24, rotation=pi / 2)

        if @isdefined(title_val)
            Label(g[0, :]; text=title_val, fontsize=32)
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

    return g
end
function ADRIA.viz.outcome_map(
    rs::ResultSet,
    si::NamedDimsArray,
    factors::Vector{Symbol};
    opts::Dict=Dict(),
    fig_opts::Dict=Dict(),
    axis_opts::Dict=Dict(),
)
    f = Figure(; fig_opts...)
    g = f[1, 1] = GridLayout()
    ADRIA.viz.outcome_map!(g, rs, si, factors; opts=opts, axis_opts=axis_opts)

    return f
end

"""
    _series_convergence(g::GridPosition, Si_conv::NamedDimsArray, factors::Vector{Symbol};
        opts::Dict=Dict(:plot_overlay => true), axis_opts::Dict=Dict())

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
    Si_conv::NamedDimsArray,
    factors::Vector{Symbol};
    opts::Dict=Dict(:plot_overlay => true),
    axis_opts::Dict=Dict(),
)
    plot_overlay = get(opts, :plot_overlay, true)
    n_scenarios = Si_conv.n_scenarios
    grps = Dict(Symbol(foi_grp) => foi_grp .== factors for foi_grp in factors)

    xlabel = pop!(axis_opts, :xlabel, "N scenarios")
    ylabel = pop!(axis_opts, :ylabel, "Sensitivity")

    if :title in keys(axis_opts)
        title_val = pop!(axis_opts, :title)
    end

    if plot_overlay
        ax = Axis(g; axis_opts...)
        _colors = colors(grps)
        scenarios_confint!(
            ax,
            permutedims(Si_conv(; Si=[:lb, :median, :ub]), (3, 1, 2)),
            collect(keys(grps)),
            _colors;
            x_vals=n_scenarios,
        )
        ax.xlabel = xlabel
        ax.ylabel = ylabel

        if @isdefined(title_val)
            ax.title = title_val
        end
    else
        _colors::Dict{Symbol,Union{Symbol,RGBA{Float32}}} = colors(grps)
        _alphas::Dict{Symbol,Float64} = alphas(grps)

        n_factors::Int64 = length(factors)
        if n_factors > 30
            ArgumentError("Too many factors to plot. Maximum number supported is 30.")
        end

        n_rows, n_cols = _calc_gridsize(n_factors)
        factor_names = String.(factors)
        axs = Axis[]
        step::Int64 = 1

        for row in 1:n_rows, col in 1:n_cols
            ax::Axis = Axis(g[row, col]; title=factor_names[step], axis_opts...)

            lines!(
                ax,
                n_scenarios,
                Si_conv(; Si=:median)(; factors=factors[step]);
                color=(_colors[factors[step]], _alphas[factors[step]]),
            )
            band!(
                ax,
                n_scenarios,
                Si_conv(; Si=:lb)(; factors=factors[step]),
                Si_conv(; Si=:ub)(; factors=factors[step]);
                color=(_colors[factors[step]], _alphas[factors[step]]),
            )
            step += 1
            push!(axs, ax)
            if step > n_factors
                break
            end
        end

        if n_factors > 1
            linkyaxes!(axs...)
            Label(g[n_rows + 1, :]; text=xlabel, fontsize=24)
            Label(g[:, 0]; text=ylabel, fontsize=24, rotation=pi / 2)

            if @isdefined(title_val)
                Label(g[0, :]; text=title_val, fontsize=32)
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
    _heatmap_convergence(g::GridPosition, Si_conv::NamedDimsArray, factors::Vector{Symbol};
        opts::Dict=Dict(), axis_opts::Dict=Dict())

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
    Si_conv::NamedDimsArray,
    factors::Vector{Symbol};
    opts::Dict=Dict(),
    axis_opts::Dict=Dict(),
)
    y_label = get(axis_opts, :ylabel, "Factors")
    x_label = get(axis_opts, :xlabel, "N scenarios")
    y_labelsize = get(axis_opts, :ylabelsize, 22)
    x_labelsize = get(axis_opts, :xlabelsize, 22)

    z = Array(Si_conv(; Si=:median))
    xtick_vals = (1:length(Si_conv.n_scenarios), string.(Si_conv.n_scenarios))
    ytick_vals = (1:length(factors), string.(factors))

    ax = Axis(
        g[1, 1];
        xticks=xtick_vals,
        yticks=ytick_vals,
        xlabel=x_label,
        ylabel=y_label,
        xlabelsize=x_labelsize,
        ylabelsize=y_labelsize,
        axis_opts...,
    )
    heatmap!(ax, z')
    colorbar_label = get(opts, :colorbar_label, "Relative Sensitivity")
    color_map = get(opts, :color_map, :viridis)

    Colorbar(g[1, 2]; colormap=color_map, label=colorbar_label, height=Relative(0.65))
    return g
end

"""
    ADRIA.viz.convergence(Si_conv::NamedDimsArray, factors::Vector{Symbol}; series_opts::Dict=Dict(),
        axis_opts::Dict=Dict())
    ADRIA.viz.convergence!(f::Figure, Si_conv::NamedDimsArray, factors::Vector{Symbol}; series_opts::Dict=Dict(),
        axis_opts::Dict=Dict())

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
    Si_conv::NamedDimsArray,
    factors::Vector{Symbol};
    opts::Dict=Dict(:viz_type => :series),
    fig_opts::Dict=Dict(),
    axis_opts::Dict=Dict(),
)
    f = Figure(; fig_opts...)
    ADRIA.viz.convergence!(
        f[1, 1],
        Si_conv,
        factors;
        opts=opts,
        axis_opts=axis_opts,
    )
    return f
end
function ADRIA.viz.convergence!(
    g::Union{GridLayout,GridPosition},
    Si_conv::NamedDimsArray,
    factors::Vector{Symbol};
    opts::Dict=Dict(:viz_type => :series),
    axis_opts::Dict=Dict(),
)
    viz_type = get(opts, :viz_type, :series)
    if viz_type == :series
        return _series_convergence(g, Si_conv, factors; opts=opts, axis_opts=axis_opts)
    elseif viz_type == :heatmap
        return _heatmap_convergence(g, Si_conv, factors; opts=opts, axis_opts=axis_opts)
    else
        error("Convergence plot $(viz_type) is not expected.")
    end

end
