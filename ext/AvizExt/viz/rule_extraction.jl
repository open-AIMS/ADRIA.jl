using SIRUS
import ADRIA.analysis: Rule

"""
    ADRIA.viz.rules_scatter(rs::ResultSet, scenarios::DataFrame, clusters::Vector{Int64}, outcomes::NamedDimsArray, rules::Vector{Rule{Vector{Vector},Vector{Float64}}}; fig_opts::Dict=Dict(), axis_opts::Dict=Dict())

# Description
For each condition Rule, plots a scatter showing one condition clause at each axis.
The target area of each condition is algo highlited, and the positive and negative
class data have different colors. This maybe be usefull for bumphunting.

# Arguments
- `rs` : ResultSet
- `scenarios` : Scenarios used to generate ResultSet
- `clusters` : Clusters used to extract Rule using SIRUS
- `outcomes` : Results of scenario metric
- `rules` : Rules extracted from scenarios and clusters

# Returns
Figure with condition Rule scatter plots
"""
function ADRIA.viz.rules_scatter(rs::ResultSet, scenarios::DataFrame, clusters::Vector{Int64}, outcomes::NamedDimsArray, rules::Vector{Rule{Vector{Vector},Vector{Float64}}}; fig_opts::Dict=Dict(), axis_opts::Dict=Dict())
    f = Figure(; fig_opts...)
    g = f[1, 1] = GridLayout()

    # For now we are only plotting conditions with two clauses
    rules = filter(r -> length(r.condition) == 2, rules)

    ADRIA.viz.rules_scatter!(g, rs, scenarios, clusters, outcomes, rules; axis_opts=axis_opts)

    return f
end
function ADRIA.viz.rules_scatter!(g::Union{GridLayout,GridPosition}, rs::ResultSet, scenarios::DataFrame, clusters::Vector{Int64}, outcomes::NamedDimsArray, rules::Vector{Rule{Vector{Vector},Vector{Float64}}}; axis_opts::Dict=Dict())
    sub_g = g[1, 1] = GridLayout()

    # Target cluster index
    target_index = clusters .== ADRIA.analysis.target_cluster(clusters, outcomes)

    # To get human readable names later
    ms = model_spec(rs)

    n_factors = length(rules)
    n_rows, n_cols = _calc_gridsize(n_factors)

    for r in 1:n_rows
        for c in 1:n_cols
            # Get condition clauses to be shown in x and y axis
            index = c + (r - 1) * n_cols
            if length(rules) < index
                break
            end
            condition = rules[index].condition

            # Feature names
            condition_features = first.(condition)
            foi = ms.fieldname .∈ [condition_features]
            h_names = ms[foi, :name]

            ax::Axis = Axis(
                sub_g[r, c],
                xlabel=h_names[1],
                ylabel=h_names[2],
                title=_sub_title(h_names, condition);
                axis_opts...
            )

            # Features to be plotted
            x_features, y_features = [scenarios[:, first(c)] for c in condition]

            # Set inferior and superior x and y limits
            x_min, x_max = _find_limits(x_features)
            y_min, y_max = _find_limits(y_features)
            xlims!(x_min, x_max), ylims!(y_min, y_max)

            for t in unique(target_index)
                x = x_features[t.==target_index]
                y = y_features[t.==target_index]
                cat_color = t == 1 ? :blue : :red
                scatter!(ax, x, y, color=cat_color, marker=:circle, markersize=4)
            end

            # Draw lines at clause breakpoints
            vlines!(ax, [last(condition[1])], color=:black)
            hlines!(ax, [last(condition[2])], color=:black)

            # Highlight target area
            poly!(ax, _target_rect(scenarios, condition), color=(:black, 0.2))
        end
    end

    g
end

function _find_limits(features)
    delta = (maximum(features) - minimum(features)) * 0.1
    [minimum(features) - delta, maximum(features) + delta]
end

function _sub_title(h_names, condition)
    inequalities = [c[2] == :L ? " < " : " ≥ " for c in condition]
    values = string.([round(c[3]; digits=2) for c in condition])
    join(h_names .* inequalities .* values, "\n")
end

function _target_rect(features, condition)
    # Inferior and superior limits for x and y features (in that order)
    lims = [_find_limits(features[:, first(c)]) for c in condition]

    # Compute target areas parameters
    x, y = [c[2] == :L ? l[1] : c[3] for (c, l) in zip(condition, lims)]
    w, h = [c[2] == :L ? c[3] - l[1] : l[2] - c[3] for (c, l) in zip(condition, lims)]

    Rect(x, y, w, h)
end
