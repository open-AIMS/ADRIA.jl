"""
    

# Arguments
- `rs` : ResultSet
- `y` : scenario outcomes
- `opts` : Aviz options  
    - `by_RCP`, color by RCP otherwise color by scenario type. Defaults to false.
- `fig_opts` : Additional options to pass to adjust Figure creation  
  See: https://docs.makie.org/v0.19/api/index.html#Figure
- `axis_opts` : Additional options to pass to adjust Axis attributes  
  See: https://docs.makie.org/v0.19/api/index.html#Axis
- `series_opts` : Additional options to pass to adjust Series attributes  
  See: https://docs.makie.org/v0.19/api/index.html#series!

# Returns
Figure
"""


"""
    scenario(rs::ADRIA.ResultSet, y::NamedDimsArray; opts=Dict(by_RCP => false), fig_opts=Dict(), axis_opts=Dict(), series_opts=Dict())    
    scenario!(f::Union{GridLayout,GridPosition}, rs::ADRIA.ResultSet, y::NamedDimsArray; opts=Dict(by_RCP => false), axis_opts=Dict(), series_opts=Dict())

Plot scenario outcomes over time.

# Arguments
- `rs` : ResultSet
- `y` : results of scenario metric
- `opts` : Aviz options 
    - `by_RCP`, color by RCP otherwise color by scenario type. Defaults to false.
- `axis_opts` : Additional options to pass to adjust Axis attributes  
  See: https://docs.makie.org/v0.19/api/index.html#Axis
- `series_opts` : Additional options to pass to adjust Series attributes  
  See: https://docs.makie.org/v0.19/api/index.html#series!

# Returns
GridPosition
"""
function scenario!(g::Union{GridLayout,GridPosition}, rs::ResultSet, y::NamedDimsArray;
    opts::Dict=Dict(:by_RCP => false), axis_opts::Dict=Dict(), series_opts::Dict=Dict())

    # Ensure last year is always shown in x-axis
    xtick_vals = get(axis_opts, :xticks, _time_labels(timesteps(y)))
    xtick_rot = get(axis_opts, :xticklabelrotation, 2 / π)

    ax = Axis(
        g[1, 1],
        xticks=xtick_vals,
        xticklabelrotation=xtick_rot;
        axis_opts...
    )

    min_step = (1 / 0.05)
    color_weight = min((1.0 / (size(y, 1) / min_step)), 0.6)

    if :color ∉ keys(series_opts)
        hide_idx = get(series_opts, :hide_series, BitVector())

        if get(opts, :by_RCP, false) == false
            series_opts = merge(series_opts, Dict(:color => scenario_colors(rs, color_weight, hide_idx)))

            cf = LineElement(color=COLORS[:counterfactual], linestyle=nothing)
            ug = LineElement(color=COLORS[:unguided], linestyle=nothing)
            gu = LineElement(color=COLORS[:guided], linestyle=nothing)

            eles = [cf, ug, gu]
            labels = ["No Intervention", "Unguided", "Guided"]
        else
            rcp_ids = sort(Int.(unique(rs.inputs[:, :RCP])))
            c = [COLORS[_r] for _r in [Symbol("RCP$(r_id)") for r_id in rcp_ids]]
            r_s = Symbol[Symbol("RCP$(r_id)") for r_id in rcp_ids]

            eles = LineElement[
                LineElement(color=_c, linestyle=nothing)
                for _c in [COLORS[_r] for _r in r_s]
            ]

            series_opts = merge(series_opts, Dict(:color => map(x -> (COLORS[Symbol("RCP$(Int(x))")], color_weight), rs.inputs[:, :RCP])))
            labels = String.(r_s)
        end

        # Add legend
        Legend(g[1, 2], eles, labels, halign=:left, valign=:top, margin=(10, 10, 10, 10))
    end

    ls = series!(ax, y'; series_opts...)
    # ax.ylabel = metric_label(metric)
    ax.xlabel = "Year"

    return g
end
function scenario(rs::ResultSet, y::NamedDimsArray; opts::Dict=Dict(:by_RCP => false), fig_opts::Dict=Dict(), axis_opts::Dict=Dict(), series_opts::Dict=Dict())
    f = Figure(; fig_opts...)
    g = f[1, 1] = GridLayout()
    scenario!(g, rs, y; opts, axis_opts, series_opts)

    return f
end


# """
#     scenario!(f::GridPosition, rs::ADRIA.ResultSet, metric, metric_args=Dict(); opts=Dict(by_RCP => false), axis_opts=Dict(), series_opts=Dict())

# Add figure to a given figure GridPosition.

# # Arguments
# - `f` : Makie figure grid position
# - `rs` : ResultSet
# - `metric` : Callable, metric to display
# - `metric_args` : Additional arguments to pass to `metric`
# - `opts` : Aviz options 
#     - `by_RCP`, color by RCP otherwise color by scenario type. Defaults to false.
# - `axis_opts` : Additional options to pass to adjust Axis attributes  
#   See: https://docs.makie.org/v0.19/api/index.html#Axis
# - `series_opts` : Additional options to pass to adjust Series attributes  
#   See: https://docs.makie.org/v0.19/api/index.html#series!

# # Returns
# GridPosition
# """
# function scenario!(f::GridPosition, rs::ResultSet, metric, metric_args::Dict=Dict();
#     opts::Dict=Dict(:by_RCP => false), axis_opts::Dict=Dict(), series_opts::Dict=Dict())

#     # Ensure last year is always shown in x-axis
#     ts = timesteps(rs)
#     xtick_vals = (
#         push!(ts[1:5:end], ts[end]),  # xtick positions
#         push!(string.(ts[1:5:end]), string(ts[end]))  # displayed vals
#     )
#     xtick_vals = get(axis_opts, :xticks, xtick_vals)
#     xtick_rot = get(axis_opts, :xticklabelrotation, 2 / π)

#     ax = Axis(
#         f,
#         xticks=xtick_vals,
#         xticklabelrotation=xtick_rot;
#         axis_opts...
#     )
#     met_vals = metric(rs; metric_args...)

#     min_step = (1 / 0.05)
#     color_weight = min((1.0 / (size(met_vals, 1) / min_step)), 0.6)

#     if :color ∉ keys(series_opts)
#         hide_idx = get(series_opts, :hide_series, BitVector())

#         if get(opts, :by_RCP, false) == false
#             series_opts = merge(series_opts, Dict(:color => scenario_colors(rs, color_weight, hide_idx)))

#             cf = LineElement(color=COLORS[:counterfactual], linestyle=nothing)
#             ug = LineElement(color=COLORS[:unguided], linestyle=nothing)
#             g = LineElement(color=COLORS[:guided], linestyle=nothing)

#             eles = [cf, ug, g]
#             labels = ["No Intervention", "Unguided", "Guided"]
#         else
#             rcp_ids = sort(Int.(unique(rs.inputs[:, :RCP])))
#             c = [COLORS[_r] for _r in [Symbol("RCP$(r_id)") for r_id in rcp_ids]]
#             r_s = [Symbol("RCP$(r_id)") for r_id in rcp_ids]

#             eles = [
#                 LineElement(color=_c, linestyle=nothing)
#                 for _c in [COLORS[_r] for _r in r_s]
#             ]

#             series_opts = merge(series_opts, Dict(:color => map(x -> (COLORS[Symbol("RCP$(Int(x))")], color_weight), rs.inputs[:, :RCP])))
#             labels = String.(r_s)
#         end

#         # Add legend
#         Legend(f[1, 2], eles, labels)
#     end

#     ls = series!(ax, ts, met_vals'; series_opts...)
#     ax.ylabel = metric_label(metric)
#     ax.xlabel = "Year"

#     return f
# end

# """
#     scenario(rs::ADRIA.ResultSet, metric, metric_args=Dict(); opts=Dict(by_RCP => false), fig_opts=Dict(), axis_opts=Dict(), series_opts=Dict())

# # Arguments
# - `rs` : ResultSet
# - `metric` : Callable, metric to display
# - `metric_args` : Additional arguments to pass to `metric`
# - `opts` : Aviz options  
#     - `by_RCP`, color by RCP otherwise color by scenario type. Defaults to false.
# - `fig_opts` : Additional options to pass to adjust Figure creation  
#   See: https://docs.makie.org/v0.19/api/index.html#Figure
# - `axis_opts` : Additional options to pass to adjust Axis attributes  
#   See: https://docs.makie.org/v0.19/api/index.html#Axis
# - `series_opts` : Additional options to pass to adjust Series attributes  
#   See: https://docs.makie.org/v0.19/api/index.html#series!

# # Returns
# Figure
# """
# function scenario(rs::ResultSet, metric, metric_args::Dict=Dict();
#     opts::Dict=Dict(:by_RCP => false), fig_opts::Dict=Dict(), axis_opts::Dict=Dict(), series_opts::Dict=Dict())
#     f = Figure(; fig_opts...)
#     scenario!(f[1, 1], rs, metric, metric_args; opts, axis_opts, series_opts)

#     return f
# end
