"""
    pawn!(f::GridPosition, Si::NamedMatrix; opts::Dict=Dict(), axis_opts::Dict=Dict())

# Arguments
- `f` : Figure GridPosition
- `Si` : Sensitivity analysis results from `pawn()`
- `opts` : Additional figure customization options
- `axis_opts` : Additional options to pass to adjust Axis attributes  
  See: https://docs.makie.org/v0.19/api/index.html#Axis

# Returns
GLMakie figure
"""
function pawn!(f::GridPosition, Si::NamedMatrix; opts::Dict=Dict(), axis_opts::Dict=Dict())
    xtick_rot = get(axis_opts, :xticklabelrotation, 2.0 / π)

    norm = get(opts, :normalize, true)
    if norm
        Si = col_normalize(Si)
    end

    foi = get(opts, :factors, :all)
    if foi != :all
        Si = Si[string.(foi), :]
    end

    y, x = names(Si)
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
    pawn(Si::NamedMatrix; opts::Dict=Dict(), fig_opts::Dict=Dict(), axis_opts::Dict=Dict())

Display heatmap of sensitivity analysis.

# Arguments
- `Si` : Results from sensitivity analysis
- `opts` : Additional figure customization options  
           normalize : column-wise normalization of values prior to figure display
           factors   : list (of strings) naming factors of interest to display from `Si`
- `fig_opts` : Additional options to pass to adjust Figure attributes  
  See: https://docs.makie.org/v0.19/api/index.html#Figure
- `axis_opts` : Additional options to pass to adjust Axis attributes  
  See: https://docs.makie.org/v0.19/api/index.html#Axis

# Returns
GLMakie figure
"""
function pawn(Si::NamedMatrix; opts::Dict=Dict(), fig_opts::Dict=Dict(), axis_opts::Dict=Dict())
    res = get(fig_opts, :resolution, (800, 500))
    f = Figure(resolution=res)
    pawn!(f[1, 1], Si; opts, axis_opts)

    return f
end
