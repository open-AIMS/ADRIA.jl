using ADRIA.analysis: col_normalize

function pairplot!(display, outcomes::NamedTuple)
    n_outcomes = length(outcomes)
    for (row, (r_m, r_out)) in zip(1:n_outcomes, pairs(outcomes))  # rows
        for (col, (c_m, c_out)) in zip(1:n_outcomes, pairs(outcomes))  # cols
            t = Axis(display[row, col])
            if col == 1
                t.ylabel = "$r_m"
            end

            if 1 <= col <= n_outcomes
                if row == 1 & col != 2
                    # show y-axis on second plot on first row
                    hideydecorations!(t; grid=false, ticks=false)
                elseif row > 1 && col > 1
                    # show y-axis only in first column
                    hideydecorations!(t; grid=false, ticks=false)
                end
            end

            if row < n_outcomes
                # Hide x-axis for all except for the last row
                # and second-last row of the final column
                if col == n_outcomes
                    if row != n_outcomes - 1
                        hidexdecorations!(t; grid=false, ticks=false)
                    end
                else
                    hidexdecorations!(t; grid=false, ticks=false)
                end
            end

            if row == n_outcomes
                t.xlabel = "$c_m"
            end

            if r_m == c_m
                hidedecorations!(t; label=false)
                hidespines!(t)
                continue
            end

            scatter!(t, vec(c_out), vec(r_out))
        end
    end
end

function pairplot!(display, data, names)
    n_outcomes = size(data, 2)

    # Plot each column against each other
    for (row, r_m) in enumerate(names)
        r_out = data[:, row]
        for (col, c_m) in enumerate(names)
            c_out = data[:, col]
            t = Axis(display[row, col])

            if col == 1
                t.ylabel = "$r_m"
            end

            if row == 1 && col != 2
                # show y-axis on second plot on first row
                if col > 1
                    hideydecorations!(t; grid=false, ticks=false)
                end
            elseif row > 1 && col > 1
                # show y-axis only in first column
                hideydecorations!(t; grid=false, ticks=false)
            end

            if row < n_outcomes
                # Hide x-axis for all except for the last row
                # and second-last row of the final column
                if col == n_outcomes
                    if row != n_outcomes - 1
                        hidexdecorations!(t; grid=false, ticks=false)
                    end
                else
                    hidexdecorations!(t; grid=false, ticks=false)
                end
            end

            if row == n_outcomes
                t.xlabel = "$c_m"
            end

            if r_m == c_m
                hidedecorations!(t; label=false)
                hidespines!(t)
                continue
            end

            scatter!(t, vec(c_out), vec(r_out); markersize=1)
        end
    end

    rowgap!(display, 20)
    colgap!(display, 20)
    return nothing
end

"""
    plot_poly!(disp, xy)

Display map elements (hacky fix until GeoMakie is updated)
"""
function plot_poly!(disp, xy)
    if length(xy[1]) == 2 && typeof(xy[1][1]) <: Real
        poly!(disp, Point2f0[xy...])  # need to work out colors based on metric/score
    else
        for _xy in xy
            plot_poly!(disp, _xy)
        end
    end
end

# """
# Parallel Coordinate Plot
# """
# function pcp!(disp, outcomes::NamedTuple; c_scheme=:inferno)
#     data = hcat(values(outcomes)...)
#     out_names = keys(outcomes)
#     k, n = size(data)

#     # limits = extrema.(eachcol(data))
#     # scaled = hcat([(d .- mi) ./ (ma - mi) for (d, (mi, ma)) in zip(eachcol(data), limits)]...)

#     for i in 1:k
#         # get(c_scheme, (i - 1) / (k - 1))
#         # show_axis = false,
#         # color = :inferno

#         # (d[i] - l[1]) ./ (l[2] - l[1])
#         lines!(disp, scaled[i, :], color=(:blue, 0.1))
#     end

#     vlines!(disp, 1:n; color=:black)
#     disp.xticks = (1:n, [string.(out_names)...])

#     return disp
# end

"""
Parallel Coordinate Plot

# Arguments
- ax : Axis
- data : Matrix
- names : Vector or Tuple of names (for x-axis)
- color : Color tuple or Vector of colors
"""
function pcp!(ax, data, names::Union{Vector,Tuple}; color=(:blue, 0.1))
    n = size(data, 2)

    vlines!(ax, 1:n; color=:black)
    series!(ax, 1:n, data; color=color)

    ax.xticks = (1:n, [string.(names)...])
    ax.xticklabelrotation = 0.45
    hideydecorations!(ax)

    return ax
end
