function ADRIA.viz.spatial_clustering(rs::Union{Domain,ResultSet},
    data::AbstractVector,
    clusters::Vector{Int64};
    opts::Dict=Dict(),
    fig_opts::Dict=Dict(),
    axis_opts::Dict=Dict())
    return spatial_clustering!(rs, data, clusters; opts=opts)
end
function spatial_clustering!(rs, data, clusters; opts::Dict=Dict(),
    fig_opts::Dict=Dict(),
    axis_opts::Dict=Dict())

    opts[:highlight] = get(opts, :highlight, cluster_colors(clusters))

    f = Figure(; fig_opts...)
    g = f[1, 1] = GridLayout()
    ADRIA.viz.map!(g, rs, collect(data); opts=opts)

    return f
end

function cluster_colors(clusters::Vector{Int64})
    n_clusters = length(unique(filter(c -> c != 0, clusters)))
    cat_colors = categorical_colors(:seaborn_bright, n_clusters)
    return [c == 0 ? "black" : cat_colors[c] for c in clusters]
end
