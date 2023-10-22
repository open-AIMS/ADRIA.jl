using Statistics
using Printf


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
GLMakie figure
"""
function ADRIA.viz.pawn!(g::Union{GridLayout,GridPosition}, Si::NamedDimsArray; opts::Dict=Dict(), axis_opts::Dict=Dict())
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

    heatmap!(ax, Matrix(Si'); colorrange=(0.0, maximum(Si)))
    # Colorbar(g[:, end+1], ax)

    return g
end
function ADRIA.viz.pawn(Si::NamedDimsArray; opts::Dict=Dict(), fig_opts::Dict=Dict(), axis_opts::Dict=Dict())
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
GLMakie figure
"""
function ADRIA.viz.tsa!(g::Union{GridLayout,GridPosition}, rs::ResultSet, si::NamedDimsArray; opts, axis_opts)
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
    lns = Combined[
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
function ADRIA.viz.tsa(rs::ResultSet, si::NamedDimsArray; opts::Dict=Dict(), fig_opts::Dict=Dict(), axis_opts::Dict=Dict())
    f = Figure(; fig_opts...)
    g = f[1, 1] = GridLayout()
    ADRIA.viz.tsa!(g, rs, si; opts, axis_opts)

    return f
end

"""
    ADRIA.viz.rsa(rs::ResultSet, si::NamedDimsArray, factors::Vector{String}; opts::Dict=Dict(), fig_opts::Dict=Dict(), axis_opts::Dict=Dict())
    ADRIA.viz.rsa!(f::Union{GridLayout,GridPosition}, rs::ResultSet, si::NamedDimsArray, factors::Vector{String}; opts, axis_opts)

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
GLMakie figure
"""
function ADRIA.viz.rsa!(g::Union{GridLayout,GridPosition}, rs::ResultSet, si::NamedDimsArray, factors::Vector{String}; opts, axis_opts)
    n_factors::Int64 = length(factors)
    if n_factors > 30
        ArgumentError("Too many factors to plot. Maximum number supported is 30.")
    end

    n_rows, n_cols = _calc_gridsize(n_factors)

    xtick_rot = get(axis_opts, :xticklabelrotation, 2.0 / π)
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
    bounds = ms[foi, :bounds]

    # Hacky special case handling for SSP/RCP
    if :RCP in factors || :SSP in factors
        loc = first(findall((factors .== :RCP) .|| (factors .== :SSP)))
        insert!(all_comps, loc, "EnvironmentalLayer")
        insert!(f_names, loc, "RCP")
        insert!(h_names, loc, "SSP/RCP")
        insert!(bounds, loc, (1, length(unique(rs.inputs.RCP))))
    end

    # comps = unique(all_comps)
    # dc = distinguishable_colors(length(comps), [RGB(1, 1, 1), RGB(0, 0, 0)], dropseed=true)
    bin_slices, factor_list = axiskeys(si)
    b_slices = parse.(Float64, bin_slices)
    curr::Int64 = 1
    axs = Axis[]
    for r in 1:n_rows
        for c in 1:n_cols
            f_vals = rs.inputs[:, Symbol(factors[curr])]
            fv_s = quantile(f_vals, b_slices)
            # fv_s = String[(i == 1) || iseven(i) ? @sprintf("%.1f", fv) : "" for (i, fv) in enumerate(quantile(f_vals, b_slices))]
            # xtick_labels = (1:length(bin_slices), fv_s)

            ax::Axis = Axis(
                g[r, c],
                title=h_names[f_names.==factors[curr]][1];
                axis_opts...
            )

            scatterlines!(ax, fv_s, si(factors=Symbol(factors[curr])), markersize=15)
            push!(axs, ax)
            curr += 1

            if curr > n_factors
                break
            end
        end
    end

    linkyaxes!(axs...)
    Label(g[end+1, :], text=xlabel, fontsize=32)
    Label(g[1:end-1, 0], text=ylabel, fontsize=32, rotation=pi / 2)

    if :title in keys(axis_opts)
        Label(g[0, :], text=title_val, fontsize=40)
    end

    # Clear empty figures
    trim!(g)

    return g
end
function ADRIA.viz.rsa(rs::ResultSet, si::NamedDimsArray, factors::Vector{String}; opts::Dict=Dict(), fig_opts::Dict=Dict(), axis_opts::Dict=Dict())
    f = Figure(; fig_opts...)
    g = f[1, 1] = GridLayout()
    ADRIA.viz.rsa!(g, rs, si, factors; opts, axis_opts)

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
GLMakie figure
"""
function ADRIA.viz.outcome_map!(g::Union{GridLayout,GridPosition}, rs::ResultSet, outcomes::NamedDimsArray, factors::Vector{String}; opts::Dict=Dict(), axis_opts::Dict=Dict())
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
    bounds = ms[foi, :bounds]

    # Hacky special case handling for SSP/RCP
    if :RCP in factors || :SSP in factors
        loc = first(findall((factors .== :RCP) .|| (factors .== :SSP)))
        insert!(all_comps, loc, "EnvironmentalLayer")
        insert!(f_names, loc, "RCP")
        insert!(h_names, loc, "SSP/RCP")
        insert!(bounds, loc, (1, length(unique(rs.inputs.RCP))))
    end

    bin_slices, factor_list, CIs = axiskeys(outcomes)
    b_slices = parse.(Float64, bin_slices)

    if any(f_names .== :guided)
        fv_labels = ["unguided", "cf", last.(split.(string.(ADRIA.methods_mcda), "."))...]
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
                g[r, c],
                title=h_names[f_names.==factors[curr]][1];
                axis_opts...
            )

            band!(ax,
                fv_s[.!ismissing.(outcomes(factors=f_name, CI=:lower))],
                collect(skipmissing(outcomes(factors=f_name, CI=:lower))),
                collect(skipmissing(outcomes(factors=f_name, CI=:upper)))
            )
            scatterlines!(ax,
                fv_s,
                outcomes(factors=f_name, CI=:mean), markersize=15
            )

            if f_name == :guided
                ax.xticks = (fv_s,fv_labels)
                ax.xticklabelrotation = pi/4
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
        Label(g[n_rows+1, :], text=xlabel, fontsize=24)
        Label(g[:, 0], text=ylabel, fontsize=24, rotation=pi / 2)

        if @isdefined(title_val)
            Label(g[0, :], text=title_val, fontsize=32)
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
function ADRIA.viz.outcome_map(rs::ResultSet, si::NamedDimsArray, factors::Vector{String}; opts::Dict=Dict(), fig_opts::Dict=Dict(), axis_opts::Dict=Dict())
    f = Figure(; fig_opts...)
    g = f[1, 1] = GridLayout()
    ADRIA.viz.outcome_map!(g, rs, si, factors; opts, axis_opts)

    return f
end
