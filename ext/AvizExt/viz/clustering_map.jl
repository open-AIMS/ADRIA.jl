function ADRIA.viz.spatial_clustering(rs::Union{Domain,ResultSet},
    data::AbstractMatrix,
    clusters::Vector{Int64};
    opts::Dict=Dict(),
    fig_opts::Dict=Dict(),
    axis_opts::Dict=Dict())

    f = Figure(; fig_opts...)
    g = f[1, 1] = GridLayout()
    spatial_clustering!(g, rs, data, clusters; opts=opts, axis_opts=axis_opts)

    return f
end
function spatial_clustering!(g, rs, data, clusters; opts::Dict=Dict(),
    axis_opts::Dict=Dict())

    opts[:highlight] = get(opts, :highlight, cluster_colors(clusters))

    d = collect(dropdims(mean(data, dims=:timesteps), dims=:timesteps))
    ADRIA.viz.map!(g, rs, d; opts=opts)

    return g
end

function cluster_colors(clusters::Vector{Int64})
    n_clusters = length(unique(filter(c -> c != 0, clusters)))
    cat_colors = categorical_colors(:seaborn_bright, n_clusters)
    return [c == 0 ? "black" : cat_colors[c] for c in clusters]
end
