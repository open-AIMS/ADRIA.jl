"""
    ts_cluster(data::AbstractMatrix, clusters::Vector{Int64}; fig_opts::Dict=Dict(), axis_opts::Dict=Dict())
    ts_cluster!(g::Union{GridLayout,GridPosition}, data::AbstractMatrix, clusters::Vector{Int64}; axis_opts::Dict=Dict())

Visualize clustered scenarios.

- `data` : Matrix of scenario data
- `clusters` :

# Returns
Figure
"""
function ADRIA.viz.ts_cluster(data::AbstractMatrix, clusters::Vector{Int64}; fig_opts::Dict=Dict(), axis_opts::Dict=Dict())
    f = Figure(; fig_opts...)
    g = f[1, 1] = GridLayout()

    ADRIA.viz.ts_cluster!(g, data, clusters; axis_opts=axis_opts)

    return f
end
function ADRIA.viz.ts_cluster!(g::Union{GridLayout,GridPosition}, data::AbstractMatrix, clusters::Vector{Int64}; axis_opts::Dict=Dict())
    # Ensure last year is always shown in x-axis
    xtick_vals = get(axis_opts, :xticks, _time_labels(timesteps(data)))
    xtick_rot = get(axis_opts, :xticklabelrotation, 2 / π)
    ax = Axis(
        g[1, 1],
        xticks=xtick_vals,
        xticklabelrotation=xtick_rot;
        axis_opts...
    )

    # Filter clusters and data for non-zero clusters
    filtered_clusters = filter(c -> c != 0, clusters)
    filtered_data = data[:, clusters .> 0]

    # TODO: Separate into own function
    # Calculate minimum opacity
    min_step = (1 / 0.05)
    n_clusters = length(unique(filtered_clusters))

    cat_colors = categorical_colors(:seaborn_bright, n_clusters)
    leg_entry = Any[]
    for clst in unique(filtered_clusters)
        n_scens = count(filtered_clusters .== clst)
        color_weight = max(min((1.0 / (n_scens / min_step)), 0.6), 0.1)

        push!(leg_entry, series!(ax, filtered_data[:, filtered_clusters .== clst]', solid_color=(cat_colors[clst], color_weight)))
    end

    Legend(g[1,2], leg_entry, "Cluster " .* string.(1:n_clusters), framevisible=false)

    return g
end
