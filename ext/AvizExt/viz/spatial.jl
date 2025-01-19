import ArchGDAL as AG
using Graphs, GraphMakie, SimpleWeightedGraphs

using ADRIA: _get_geom_col

"""
    _get_geoms(gdf::DataFrame, geom_col::Symbol)

Retrieve the vector of geometries from a specified column.
"""
function _get_geoms(gdf::DataFrame, geom_col::Symbol)
    return GeoMakie.to_multipoly(gdf[:, geom_col])
end

"""
    _get_geoms(gdf::DataFrame)

Retrieve the vector of geometries from a GeoDataFrame.
"""
function _get_geoms(gdf::DataFrame)
    return _get_geoms(gdf, _get_geom_col(gdf))
end

function set_figure_defaults(fig_opts::OPT_TYPE)::OPT_TYPE
    fig_opts[:size] = get(fig_opts, :size, (600, 900))
    fig_opts[:xticklabelsize] = get(fig_opts, :xticklabelsize, 14)
    fig_opts[:yticklabelsize] = get(fig_opts, :yticklabelsize, 14)

    return fig_opts
end

function set_axis_defaults(axis_opts::OPT_TYPE)::OPT_TYPE
    axis_opts[:title] = get(axis_opts, :title, "Study Area")
    axis_opts[:xlabel] = get(axis_opts, :xlabel, "Longitude")
    axis_opts[:ylabel] = get(axis_opts, :ylabel, "Latitude")
    axis_opts[:xgridwidth] = get(axis_opts, :xgridwidth, 0.5)
    axis_opts[:ygridwidth] = get(axis_opts, :ygridwidth, 0.5)
    axis_opts[:dest] = get(axis_opts, :dest, "+proj=latlong +datum=WGS84")
    axis_opts[:xgridvisible] = get(axis_opts, :xgridvisible, false)
    axis_opts[:ygridvisible] = get(axis_opts, :ygridvisible, false)

    return axis_opts
end

"""
    create_map!(
        f::Union{GridLayout,GridPosition},
        geodata::Vector{<:GeoMakie.GeometryBasics.MultiPolygon},
        data::Observable,
        highlight::Union{Vector,Tuple,Nothing},
        show_colorbar::Bool=true,
        colorbar_label::String="",
        color_map::$COLORMAP_TYPE_DOCSTRING,
        legend_params::Union{Tuple,Nothing}=nothing,
        axis_opts::Dict{Symbol, <:Any}=set_axis_defaults(Dict{Symbol,Any}())
    )

Create a spatial choropleth figure.

# Arguments
- `f` : Makie figure to create plot in
- `geodata` : MultiPolygon features to display
- `data` : Values to use for choropleth
- `highlight` : Stroke colors for each location
- `show_colorbar` : Whether to show a colorbar (true) or not (false)
- `colorbar_label` : Label to use for color bar
- `color_map` : Type of colormap to use,
    See: https://docs.makie.org/stable/documentation/colors/#colormaps
- `legend_params` : Legend parameters
- `axis_opts` : Additional options to pass to adjust Axis attributes
  See: https://docs.makie.org/v0.19/api/index.html#Axis
"""
function create_map!(
    f::Union{GridLayout,GridPosition},
    geodata::Vector{<:GeoMakie.MultiPolygon},
    data::Observable,
    highlight::Union{Vector,Tuple,Nothing};
    show_colorbar::Bool=true,
    colorbar_label::String="",
    color_map::COLORMAP_TYPE(T)=:grayC,
    legend_params::Union{Tuple,Nothing}=nothing,
    axis_opts::OPT_TYPE=set_axis_defaults(DEFAULT_OPT_TYPE())
)::Union{GridLayout,GridPosition} where {T<:Real}
    spatial = GeoAxis(
        f[1, 1];
        axis_opts...
    )

    # spatial.xticklabelsize = 14
    # spatial.yticklabelsize = 14

    min_val = @lift(minimum($data))
    max_val = @lift(maximum($data))

    # Plot geodata polygons using data as internal color
    color_range = min_val[] < 0 ? (min_val[], max_val[]) : (0, max_val[])

    poly!(
        spatial,
        geodata;
        color=data,
        colormap=color_map,
        colorrange=color_range,
        strokecolor=(:black, 0.05),
        strokewidth=1.0
    )

    if show_colorbar
        Colorbar(
            f[1, 2];
            colorrange=color_range,
            colormap=color_map,
            label=colorbar_label,
            height=Relative(0.70)
        )
    end

    # Overlay locations to be highlighted
    # `poly!()` cannot handle multiple strokecolors being specified at the moment
    # so we instead overlay each cluster.
    if !isnothing(highlight)
        if highlight isa Tuple
            poly!(
                spatial,
                geodata;
                color="transparent",
                strokecolor=highlight,
                strokewidth=0.5,
                linestyle=:solid,
                overdraw=true
            )
        else
            hl_groups = unique(highlight)

            for color in hl_groups
                m = findall(highlight .== [color])

                poly!(
                    spatial,
                    geodata[m];
                    color="transparent",
                    strokecolor=color,
                    strokewidth=0.5,
                    linestyle=:solid,
                    overdraw=true
                )
            end
        end

        if !isnothing(legend_params)
            # Plot Legend only if highlight colors are present
            Legend(f[1, 3], legend_params...; framevisible=false)
        end
    end

    if typeof(f) == GridLayout
        trim!(f)
    end

    return f
end

"""
    ADRIA.viz.map(rs::Union{Domain,ResultSet}, y::Union{YAXArray,AbstractVector{<:Real}}; diverging::Bool=false, opts::OPT_TYPE=DEFAULT_OPT_TYPE(), fig_opts::OPT_TYPE=set_figure_defaults(DEFAULT_OPT_TYPE()), axis_opts::OPT_TYPE=set_axis_defaults(DEFAULT_OPT_TYPE()))
    ADRIA.viz.map(rs::Union{Domain,ResultSet}; opts::OPT_TYPE=DEFAULT_OPT_TYPE(), fig_opts::OPT_TYPE=set_figure_defaults(DEFAULT_OPT_TYPE()), axis_opts::OPT_TYPE=set_axis_defaults(DEFAULT_OPT_TYPE()))
    ADRIA.viz.map!(g::Union{GridLayout,GridPosition}, rs::Union{Domain,ResultSet}, y::AbstractVector{<:Real}; opts::OPT_TYPE=DEFAULT_OPT_TYPE(), axis_opts::OPT_TYPE=set_axis_defaults(DEFAULT_OPT_TYPE()))

Plot spatial choropleth of outcomes.

# Arguments
- `rs` : ResultSet
- `y` : results of scenario metric
- `diverging` : If true, uses a diverging color map.
- `opts` : Aviz options
    - `colorbar_label`, label for colorbar. Defaults to "Relative Cover"
    - `color_map`, preferred colormap for plotting heatmaps
- `axis_opts` : Additional options to pass to adjust Axis attributes
  See: https://docs.makie.org/v0.19/api/index.html#Axis

# Returns
GridPosition
"""
function ADRIA.viz.map(
    rs::Union{Domain,ResultSet},
    y::Union{YAXArray,AbstractVector{<:Real}};
    diverging::Bool=false,
    opts::OPT_TYPE=DEFAULT_OPT_TYPE(),
    fig_opts::OPT_TYPE=set_figure_defaults(DEFAULT_OPT_TYPE()),
    axis_opts::OPT_TYPE=set_axis_defaults(DEFAULT_OPT_TYPE())
)
    set_plot_opts!(y, opts, :colorbar_label; metadata_key=:metric_name)

    if diverging
        opts[:color_map] = _diverging_cmap(y)
    end

    f = Figure(; fig_opts...)
    g = f[1, 1] = GridLayout()

    ADRIA.viz.map!(g, rs, collect(y); opts, axis_opts)

    return f
end
function ADRIA.viz.map(
    rs::Union{Domain,ResultSet};
    opts::OPT_TYPE=DEFAULT_OPT_TYPE(),
    fig_opts::OPT_TYPE=set_figure_defaults(DEFAULT_OPT_TYPE()),
    axis_opts::OPT_TYPE=set_axis_defaults(DEFAULT_OPT_TYPE())
)
    f = Figure(; fig_opts...)
    g = f[1, 1] = GridLayout()

    opts[:colorbar_label] = get(opts, :colorbar_label, "Coral Real Estate [%]")

    opts[:show_management_zones] = get(opts, :show_management_zones, false)
    if opts[:show_management_zones]
        local highlight
        try
            highlight = Symbol.(lowercase.(rs.loc_data.zone_type))
        catch
            # Annoyingly, the case of the name may have changed...
            highlight = Symbol.(lowercase.(rs.loc_data.ZONE_TYPE))
        end
        opts[:highlight] = highlight
    end

    ADRIA.viz.map!(g, rs, rs.loc_data.k * 100.0; opts, axis_opts)

    return f
end
function ADRIA.viz.map!(
    g::Union{GridLayout,GridPosition},
    rs::Union{Domain,ResultSet},
    y::AbstractVector{<:Real};
    opts::OPT_TYPE=DEFAULT_OPT_TYPE(),
    axis_opts::OPT_TYPE=set_axis_defaults(DEFAULT_OPT_TYPE())
)
    geodata = _get_geoms(rs.loc_data)
    data = Observable(collect(y))

    highlight = get(opts, :highlight, nothing)
    c_label = get(opts, :colorbar_label, "")
    legend_params = get(opts, :legend_params, nothing)
    show_colorbar = get(opts, :show_colorbar, true)
    color_map = get(opts, :color_map, :grayC)
    if !(:dest in keys(axis_opts))
        col = _get_geom_col(rs.loc_data)
        axis_opts[:dest] = ADRIA.AG.toPROJ4(ADRIA.AG.getspatialref(rs.loc_data[1, col]))
    end

    return create_map!(
        g,
        geodata,
        data,
        highlight;
        show_colorbar=show_colorbar,
        colorbar_label=c_label,
        color_map=color_map,
        legend_params=legend_params,
        axis_opts=axis_opts
    )
end

"""
    ADRIA.viz.map(rs::Union{Domain,ResultSet}, outputs_matrix::Matrix, map_titles::Vector{String}; opts::OPT_TYPE=DEFAULT_OPT_TYPE(), fig_opts::OPT_TYPE=set_figure_defaults(DEFAULT_OPT_TYPE()), axis_opts::OPT_TYPE=set_axis_defaults(DEFAULT_OPT_TYPE()))
    ADRIA.viz.map!(g::Union{GridLayout,GridPosition}, rs::Union{Domain,ResultSet}, outputs_matrix::Matrix, map_titles::Vector{String}; opts::OPT_TYPE=DEFAULT_OPT_TYPE(), axis_opts::OPT_TYPE=set_axis_defaults(DEFAULT_OPT_TYPE()))

Plot a series of maps from an arbitrary (n_locs*n_maps) matrix of outputs.

# Arguments
- `rs` : ResultSet
- `g` : Figure GridPosition or GridLayout.
- `outputs_matrix` : Matrix of outputs where n_locs is the numberof locations and n_maps is the number of different
    maps to plot.
- `map_titles` : Titles for each map to be plotted.
- `opts` : Aviz options
    - `colorbar_label`, label for colorbar. Defaults to "Relative Cover"
    - `color_map`, preferred colormap for plotting heatmaps
- `axis_opts` : Additional options to pass to adjust Axis attributes
  See: https://docs.makie.org/v0.19/api/index.html#Axis
"""
function ADRIA.viz.map(
    rs::Union{Domain,ResultSet},
    outputs_matrix::Matrix,
    map_titles::Vector{String};
    criteria::Vector{Symbol}=Array(M.criteria),
    opts::OPT_TYPE=DEFAULT_OPT_TYPE(),
    fig_opts::OPT_TYPE=set_figure_defaults(DEFAULT_OPT_TYPE()),
    axis_opts::OPT_TYPE=set_axis_defaults(DEFAULT_OPT_TYPE())
)
    f = Figure(; fig_opts...)
    g = f[1, 1] = GridLayout()
    ADRIA.viz.map!(
        g, rs, outputs_matrix, map_titles; criteria=criteria, opts=opts, axis_opts=axis_opts
    )
    return f
end
function ADRIA.viz.map!(
    g::Union{GridLayout,GridPosition},
    rs::Union{Domain,ResultSet},
    outputs_matrix::Matrix,
    map_titles::Vector{String};
    opts::OPT_TYPE=DEFAULT_OPT_TYPE(),
    axis_opts::OPT_TYPE=set_axis_defaults(DEFAULT_OPT_TYPE())
)
    if length(rs.loc_data.site_id) != size(outputs_matrix, 1)
        error("Only unfiltered decision matrices can be plotted.")
    end

    opts[:color_map] = get(opts, :color_map, :viridis)
    opts[:colorbar_limits] = get(opts, :colorbar_limits, (0.0, 1.0))

    n_plots::Int64 = length(map_titles)
    n_rows, n_cols = _calc_gridsize(n_plots)
    step::Int64 = 1

    for row in 1:n_rows, col in 1:n_cols
        if step > length(map_titles)
            break
        end
        axis_opts[:title] = map_titles[step]
        ADRIA.viz.map!(
            g[row, col],
            rs,
            vec(outputs_matrix[:, step]);
            opts=opts,
            axis_opts=axis_opts
        )

        step += 1
    end
    return g
end

function _diverging_cmap(outcomes::YAXArray)::Vector{RGB{Float64}}
    min_val, max_val = extrema(outcomes)

    # Hande only positive or only negative value cases
    min_val = min_val > 0 ? 0 : min_val
    max_val = max_val < 0 ? 0 : max_val

    mid_val = -min_val / (max_val - min_val)
    return diverging_palette(10, 200; mid=mid_val)
end

"""
    ADRIA.viz.map(gdf::DataFrame; geom_col=:geometry, color=nothing)

Plot an arbitrary GeoDataFrame, optionally specifying a color for each feature.

# Arguments
- `gdf` : GeoDataFrame to plot
- `geom_col` : Column in GeoDataFrame that holds feature data
- `color` : Colors to use for each feature
"""
function ADRIA.viz.map(gdf::DataFrame; geom_col=:geometry, color=nothing)
    f = Figure(; size=(600, 900), figure_padding=0.1)
    ga = GeoAxis(
        f[1, 1];
        dest="+proj=latlong +datum=WGS84",
        xlabel="Longitude",
        ylabel="Latitude",
        xticklabelpad=15,
        yticklabelpad=10,
        xticklabelsize=10,
        yticklabelsize=10,
        aspect=DataAspect(),
        xgridwidth=0.5,
        ygridwidth=0.5
    )

    ADRIA.viz.map!(ga, gdf; geom_col=geom_col, color=color)

    display(f)

    return f
end

"""
    ADRIA.viz.map!(gdf::DataFrame; geom_col=:geometry, color=nothing)::Nothing
"""
function ADRIA.viz.map!(gdf::DataFrame; geom_col=:geometry, color=nothing)::Nothing
    ga = current_axis()
    ADRIA.viz.map!(ga, gdf; geom_col=geom_col, color=color)

    return nothing
end
function ADRIA.viz.map!(
    ga::GeoAxis, gdf::DataFrame; geom_col=:geometry, color=nothing
)::Nothing
    plottable = _get_geoms(gdf, geom_col)

    if !isnothing(color)
        poly!(ga, plottable; color=color)
    else
        poly!(ga, plottable)
    end

    return nothing
end

"""
    ADRIA.viz.connectivity(dom::Domain; in_method=nothing, out_method=eigenvector_centrality, opts::Dict{Symbol, <:Any}=Dict{Symbol, <:Any}(), fig_opts::Dict{Symbol, <:Any}=set_figure_defaults(Dict{Symbol,Any}()), axis_opts::Dict{Symbol, <:Any}=set_axis_defaults(Dict{Symbol,Any}()))
    ADRIA.viz.connectivity(dom::Domain, conn::AbstractMatrix; in_method=nothing, out_method=eigenvector_centrality, opts::Dict{Symbol, <:Any}=Dict{Symbol, <:Any}(), fig_opts::Dict{Symbol, <:Any}=set_figure_defaults(Dict{Symbol,Any}()), axis_opts::Dict{Symbol, <:Any}=set_axis_defaults(Dict{Symbol,Any}()))
    ADRIA.viz.connectivity(dom::Domain, network::SimpleWeightedDiGraph, conn_weights::AbstractVector{<:Real}; opts::Dict{Symbol, <:Any}=Dict{Symbol, <:Any}(), fig_opts::Dict{Symbol, <:Any}=set_figure_defaults(Dict{Symbol,Any}()), axis_opts::Dict{Symbol, <:Any}=set_axis_defaults(Dict{Symbol,Any}()))
    ADRIA.viz.connectivity!(g::Union{GridLayout, GridPosition}, dom::Domain,  network::SimpleWeightedDiGraph, conn_weights::AbstractVector{<:Real}; opts::Dict{Symbol, <:Any}=Dict{Symbol, <:Any}(), axis_opts::Dict{Symbol, <:Any}=set_axis_defaults(Dict{Symbol,Any}()))

Produce visualization of connectivity between reef sites with node size and edge visibility
weighted by the connectivity values and node weights.

# Examples

Basic visualization, plotting out-centralities by default:

```julia
dom = ADRIA.load_domain("<Path to Domain>")
ADRIA.viz.connectivity(dom)

# Plot indegree centrality instead
ADRIA.viz.connectivity(dom; in_method=indegree_centrality, out_method=nothing)
```

Finer grain control:

```julia
dom = ADRIA.load_domain("<Path to Domain>")

in_conn, out_conn, network = ADRIA.connectivity_strength(dom.conn; out_method=eigenvector_centrality)

# Plot in centrality
ADRIA.viz.connectivity(
    dom,
    network,
    in_conn;
    opts=opts,
    fig_opts=fig_opts,
    axis_opts=axis_opts
)

# Plot out centrality
ADRIA.viz.connectivity(
    dom,
    network,
    out_conn;
    opts=opts,
    fig_opts=fig_opts,
    axis_opts=axis_opts
)
```

# Arguments
- `dom` : Domain
- `network` : SimpleWeightedDiGraph calculated from the connectivity matrix
- `conn_weights` : Connectivity weighted for each node
- `opts` : AvizOpts
    - `edge_color`, vector of colours for edges. Defaults to reasonable weighting
    - `node_color`, vector of colours for node. Defaults to `conn_weights`
    - `node_size`, size of nodes in the graph
- `fig_opts` : Figure options
- `axis_opts` : Axis options
"""
function ADRIA.viz.connectivity(
    dom::Domain;
    in_method=nothing,
    out_method=eigenvector_centrality,
    opts::OPT_TYPE=DEFAULT_OPT_TYPE(),
    fig_opts::OPT_TYPE=set_figure_defaults(DEFAULT_OPT_TYPE()),
    axis_opts::OPT_TYPE=set_axis_defaults(DEFAULT_OPT_TYPE())
)
    return ADRIA.viz.connectivity(
        dom, dom.conn; in_method, out_method, opts, fig_opts, axis_opts
    )
end
function ADRIA.viz.connectivity(
    dom::Domain,
    conn::AbstractMatrix;
    in_method=nothing,
    out_method=eigenvector_centrality,
    opts::OPT_TYPE=DEFAULT_OPT_TYPE(),
    fig_opts::OPT_TYPE=set_figure_defaults(DEFAULT_OPT_TYPE()),
    axis_opts::OPT_TYPE=set_axis_defaults(DEFAULT_OPT_TYPE())
)
    if !isnothing(in_method) && !isnothing(out_method)
        @warn "Both in and out centrality measures provided. Plotting out centralities."
        _, conn_weight, network = ADRIA.connectivity_strength(conn; in_method, out_method)
    elseif !isnothing(in_method) && isnothing(out_method)
        conn_weight, _, network = ADRIA.connectivity_strength(
            conn; in_method, out_method=outdegree_centrality
        )
    elseif isnothing(in_method) && isnothing(out_method)
        error("Measure for in or out centralities needs to be provided.")
    else
        if isnothing(in_method)
            in_method = indegree_centrality
        end

        _, conn_weight, network = ADRIA.connectivity_strength(
            conn; in_method=in_method, out_method=outdegree_centrality
        )
    end

    return ADRIA.viz.connectivity(dom, network, conn_weight; opts, fig_opts, axis_opts)
end
function ADRIA.viz.connectivity(
    dom::Domain,
    network::SimpleWeightedDiGraph,
    conn_weights::AbstractVector{<:Real};
    opts::OPT_TYPE=DEFAULT_OPT_TYPE(),
    fig_opts::OPT_TYPE=set_figure_defaults(DEFAULT_OPT_TYPE()),
    axis_opts::OPT_TYPE=set_axis_defaults(DEFAULT_OPT_TYPE())
)
    f = Figure(; fig_opts...)
    g = f[1, 1] = GridLayout()

    ADRIA.viz.connectivity!(g, dom, network, conn_weights; opts, axis_opts)

    return f
end
function ADRIA.viz.connectivity!(
    g::Union{GridLayout,GridPosition},
    dom::Domain,
    network::SimpleWeightedDiGraph,
    conn_weights::AbstractVector{<:Real};
    opts::OPT_TYPE=DEFAULT_OPT_TYPE(),
    axis_opts::OPT_TYPE=DEFAULT_OPT_TYPE(Dict{Symbol,Any}())
)
    axis_opts[:title] = get(axis_opts, :title, "Study Area")
    axis_opts[:xlabel] = get(axis_opts, :xlabel, "Longitude")
    axis_opts[:ylabel] = get(axis_opts, :ylabel, "Latitude")

    spatial = GeoAxis(
        g[1, 1];
        dest="+proj=latlong +datum=WGS84",
        axis_opts...
    )

    geodata = _get_geoms(dom.loc_data)

    spatial.yticklabelpad = 50
    spatial.ytickalign = 10

    # Calculate alpha values for edges based on connectivity strength and weighting
    edge_col = Vector{Tuple{Symbol,Float64}}(undef, ne(network))
    norm_coef = maximum(conn_weights)
    for (ind, e) in enumerate(edges(network))
        if (e.src == e.dst)
            edge_col[ind] = (:black, 0.0)
            continue
        end
        edge_col[ind] = (:black, conn_weights[e.src] * e.weight / norm_coef)
    end

    # Rescale node size to be visible
    node_size = conn_weights ./ maximum(conn_weights) .* 10.0

    # Extract graph kwargs and set defaults
    edge_col = get(opts, :edge_color, edge_col)
    node_size = get(opts, :node_size, node_size)
    node_color = get(opts, :node_color, node_size)

    # Plot geodata polygons
    poly!(
        spatial,
        geodata;
        color=:white,
        strokecolor=(:black, 0.25),
        strokewidth=1.0
    )

    # Plot the connectivity graph
    graphplot!(
        spatial,
        network;
        layout=ADRIA.centroids(dom),
        edge_color=edge_col,
        node_size=node_size,
        node_color=node_color,
        edge_plottype=:linesegments
    )

    return g
end
