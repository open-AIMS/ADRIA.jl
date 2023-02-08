"""
    pawn!(f::GridPosition, Si::NamedDimsArray; opts::Dict=Dict(), axis_opts::Dict=Dict())

# Arguments
- `f` : Figure GridPosition
- `Si` : Sensitivity analysis results from `pawn()`
- `opts` : Additional figure customization options
    - `normalize` : Normalize each column ∈ [0, 1] to obtain relative sensitivity
    - `factors` : List of factors to display (factors are filtered after normalization)
- `axis_opts` : Additional options to pass to adjust Axis attributes  
  See: https://docs.makie.org/v0.19/api/index.html#Axis

# Returns
GLMakie figure
"""
function pawn!(f::GridPosition, Si::NamedDimsArray; opts::Dict=Dict(), axis_opts::Dict=Dict())
    xtick_rot = get(axis_opts, :xticklabelrotation, 2.0 / π)

    norm = get(opts, :normalize, true)
    if norm
        Si = col_normalize(Si)
    end

    foi = get(opts, :factors, :all)
    if foi != :all
        Si = Si[string.(foi), :]
    end

    y, x = axiskeys(Si)
    ax = Axis(
        f,
        xticks=(1:length(x), string.(x)),
        yticks=(1:length(y), string.(y)),
        xticklabelrotation=xtick_rot;
        axis_opts...
    )
    ax.yreversed = true

    heatmap!(ax, Matrix(Si'))

    return f
end

"""
    pawn(Si::NamedDimsArray; opts::Dict=Dict(), fig_opts::Dict=Dict(), axis_opts::Dict=Dict())

Display heatmap of sensitivity analysis.

# Arguments
- `Si` : Results from sensitivity analysis
- `opts` : Additional figure customization options  
    - `normalize` : Normalize each column ∈ [0, 1] to obtain relative sensitivity
    - `factors` : List of factors to display (factors are filtered after normalization)
- `fig_opts` : Additional options to pass to adjust Figure creation  
  See: https://docs.makie.org/v0.19/api/index.html#Figure
- `axis_opts` : Additional options to pass to adjust Axis attributes  
  See: https://docs.makie.org/v0.19/api/index.html#Axis

# Returns
GLMakie figure
"""
function pawn(Si::NamedDimsArray; opts::Dict=Dict(), fig_opts::Dict=Dict(), axis_opts::Dict=Dict())
    f = Figure(; fig_opts...)
    pawn!(f[1, 1], Si; opts, axis_opts)

    return f
end


function tsa!(f::GridPosition, tsa::NamedDimsArray, stat::Symbol=:median; opts, axis_opts)
    xtick_rot = get(axis_opts, :xticklabelrotation, 2.0 / π)
    xlabel = get(axis_opts, :xlabel, "Years")
    ylabel = get(axis_opts, :ylabel, L"\text{PAWN}_\text{%$(stat)}")

    factors, Si, timesteps = axiskeys(tsa)
    x_tickpos, x_ticklabel = _time_labels(timesteps)
    ax = Axis(
        f,
        xticks=(x_tickpos, x_ticklabel),
        xticklabelrotation=xtick_rot,
        xlabel=xlabel,
        ylabel=ylabel;
        axis_opts...
    )

    series!(ax, tsa(Si=stat), labels=factors, color=distinguishable_colors(length(factors)))

    return f
end
function tsa(tsa::NamedDimsArray, stat=:median; opts::Dict=Dict(), fig_opts::Dict=Dict(), axis_opts::Dict=Dict())
    f = Figure(; fig_opts...)
    tsa!(f[1, 1], tsa, stat; opts, axis_opts)

    return f
end
