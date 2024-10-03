using SIRUS
import ADRIA.analysis: Rule

"""
    ADRIA.viz.rules_scatter(rs::ResultSet, scenarios::DataFrame, clusters::Vector{Int64}, rules::Vector{Rule{Vector{Vector}, Vector{Float64}}}; opts::Dict{Symbol,<:Any}=Dict{Symbol,Any}(), fig_opts::Dict{Symbol,<:Any}=Dict{Symbol,Any}(), axis_opts::Dict{Symbol,<:Any}=Dict{Symbol,Any}())::Figure
    ADRIA.viz.rules_scatter(rs::ResultSet, scenarios::DataFrame, clusters::BitVector, rules::Vector{Rule{Vector{Vector},Vector{Float64}}}; opts::Dict{Symbol,<:Any}=Dict{Symbol,Any}(), fig_opts::Dict{Symbol,<:Any}=Dict{Symbol,Any}(), axis_opts::Dict{Symbol,<:Any}=Dict{Symbol,Any}())::Figure
    ADRIA.viz.rules_scatter!(g::Union{GridLayout,GridPosition}, rs::ResultSet, scenarios::DataFrame, clusters::Vector{Int64}, rules::Vector{Rule{Vector{Vector},Vector{Float64}}}; opts::Dict{Symbol,<:Any}=Dict{Symbol,Any}(), axis_opts::Dict{Symbol,<:Any}=Dict{Symbol,Any}())

# Description
For each condition Rule, plots a scatter showing one condition clause at each axis.
The target area of each condition is algo highlited, and the positive and negative
class data have different colors. This maybe be usefull for bumphunting.

# Arguments
- `rs` : ResultSet
- `scenarios` : Scenarios used to generate ResultSet
- `clusters` : Vector indicating targeted (1) and non targeted (0) scenarios
- `rules` : Rules extracted from scenarios and clusters
- `opts` : Additional figure customization options
- `fig_opts` : Additional options to pass to adjust Figure creation
  See: https://docs.makie.org/v0.19/api/index.html#Figure
- `axis_opts` : Additional options to pass to adjust Axis attributes
  See: https://docs.makie.org/v0.19/api/index.html#Axis

# Returns
Figure with condition Rule scatter plots
"""
function ADRIA.viz.rules_scatter(
    rs::ResultSet,
    scenarios::DataFrame,
    clusters::Vector{Int64},
    rules::Vector{Rule{Vector{Vector},Vector{Float64}}};
    opts::OPT_TYPE=DEFAULT_OPT_TYPE(),
    fig_opts::OPT_TYPE=DEFAULT_OPT_TYPE(),
    axis_opts::OPT_TYPE=DEFAULT_OPT_TYPE()
)::Figure
    # For now we only plot conditions with two clauses
    rules = filter(r -> length(r.condition) == 2, rules)

    # Setup figure size
    n_rows = ceil(10 / 4)
    base_width = 1200
    base_height = n_rows * base_width / 5
    fig_opts[:size] = get(fig_opts, :size, (base_width, base_height))

    f = Figure(; fig_opts...)
    g = f[1, 1] = GridLayout()

    ADRIA.viz.rules_scatter!(
        g, rs, scenarios, clusters, rules; opts=opts, axis_opts=axis_opts
    )

    return f
end
function ADRIA.viz.rules_scatter(
    rs::ResultSet,
    scenarios::DataFrame,
    clusters::BitVector,
    rules::Vector{Rule{Vector{Vector},Vector{Float64}}};
    opts::OPT_TYPE=DEFAULT_OPT_TYPE(),
    fig_opts::OPT_TYPE=DEFAULT_OPT_TYPE(),
    axis_opts::OPT_TYPE=DEFAULT_OPT_TYPE()
)::Figure
    return ADRIA.viz.rules_scatter(
        rs,
        scenarios,
        Int64.(clusters),
        rules;
        opts=opts,
        fig_opts=fig_opts,
        axis_opts=axis_opts
    )
end
function ADRIA.viz.rules_scatter!(
    g::Union{GridLayout,GridPosition},
    rs::ResultSet,
    scenarios::DataFrame,
    clusters::Vector{Int64},
    rules::Vector{Rule{Vector{Vector},Vector{Float64}}};
    opts::OPT_TYPE=DEFAULT_OPT_TYPE(),
    axis_opts::OPT_TYPE=DEFAULT_OPT_TYPE()
)
    sub_g = g[1, 1] = GridLayout()

    # Colors and Labels Setup
    colors = ["#1f78b4", "#ff7f00"]
    labels = ["Target", "Non-target"]
    marker = :circle
    title_size = 11
    labels_size = 11

    n_factors = length(rules)
    n_rows, n_cols = _calc_gridsize(n_factors)
    spec = model_spec(rs)
    for r in 1:n_rows
        for c in 1:n_cols
            # Get condition clauses to be shown in x and y axis
            index = c + (r - 1) * n_cols
            length(rules) < index && break
            condition = rules[index].condition

            # Human readable feature names
            fieldnames::Vector{String} = first.(condition)
            feature_names = _feature_names(fieldnames, spec)

            ax::Axis = Axis(
                sub_g[r, c];
                xlabel=feature_names[1],
                ylabel=feature_names[2],
                title=_readable_condition(condition, feature_names),
                titlesize=title_size,
                xlabelsize=labels_size,
                ylabelsize=labels_size,
                axis_opts...
            )

            # Features to be plotted
            x_features, y_features = [scenarios[:, first(c)] for c in condition]

            # Set inferior and superior x and y limits
            x_min, x_max = _find_limits(x_features)
            y_min, y_max = _find_limits(y_features)
            xlims!(x_min, x_max), ylims!(y_min, y_max)

            for c in unique(clusters)
                x = x_features[c .== clusters]
                y = y_features[c .== clusters]
                cat_color = c == 1 ? colors[1] : colors[2]
                scatter!(ax, x, y; color=cat_color, marker=marker, markersize=4)
            end

            if get(opts, :target_area, true)
                _highlight_target_area(ax, condition, scenarios)
            end
        end
    end

    # Create Legend
    markers = MarkerElement[MarkerElement(; color=_c, marker=marker) for _c in colors]
    Legend(
        sub_g[1, n_cols + 1],
        markers,
        labels,
        "Clusters";
        halign=:left,
        valign=:top,
        margin=(5, 5, 5, 5)
    )

    Label(
        g[1, :, Top()], "SIRUS extracted rules"; padding=(0, 0, 50, 0), font=:bold,
        valign=:bottom
    )

    return g
end

function _highlight_target_area(ax::Axis, condition::Vector{Vector}, scenarios::DataFrame)
    # Draw lines at clause breakpoints
    vlines!(ax, [last(condition[1])]; color=(:black, 0.4))
    hlines!(ax, [last(condition[2])]; color=(:black, 0.4))

    # Highlight target area
    poly!(ax, _target_area(scenarios, condition); color=(:black, 0.08))
    return nothing
end

function _find_limits(features::Vector{<:Real})
    delta = (maximum(features) - minimum(features)) * 0.1
    return [minimum(features) - delta, maximum(features) + delta]
end

function _readable_condition(condition::Vector{Vector}, feature_names::Vector{String})
    inequalities = [c[2] == :L ? " < " : " ≥ " for c in condition]
    values = string.([round(c[3]; digits=2) for c in condition])
    return join(feature_names .* inequalities .* values, "\n")
end

function _target_area(scenarios::DataFrame, condition::Vector{Vector})
    # Inferior and superior limits for x and y features (in that order)
    lims = [_find_limits(scenarios[:, first(c)]) for c in condition]

    # Compute target areas parameters
    x, y = [c[2] == :L ? l[1] : c[3] for (c, l) in zip(condition, lims)]
    w, h = [c[2] == :L ? c[3] - l[1] : l[2] - c[3] for (c, l) in zip(condition, lims)]

    return Rect(x, y, w, h)
end

function _feature_names(fieldnames::Vector{String}, spec::DataFrame)::Vector{String}
    return [spec[spec.fieldname .== f, :name][1] for f in Symbol.(fieldnames)]
end
