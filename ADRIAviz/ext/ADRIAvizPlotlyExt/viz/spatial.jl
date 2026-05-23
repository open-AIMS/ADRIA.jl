import GeoInterface as GI
import ArchGDAL as AG
using DataFrames: DataFrame, nrow, hasproperty
using Graphs: edges, ne, eigenvector_centrality, indegree_centrality, outdegree_centrality
using SimpleWeightedGraphs: SimpleWeightedDiGraph
using ADRIA: Domain, ResultSet, _get_geom_col, centroids

# =============================================================================
# GeoJSON conversion helpers
#
# GeoPackage reef-site layers from ADRIA are stored in WGS84 (EPSG:4326).
# GeoInterface.coordinates(geom) returns the same nested-array structure that
# GeoJSON expects, so no JSON parsing is required.
# =============================================================================

"""
    _get_site_ids(gdf::DataFrame) -> Vector{String}

Return string site identifiers. Tries `:site_id`, then `:reef_siteid`, then
falls back to row indices. These values become both the GeoJSON feature `id`
fields and the `locations` array in a Plotly choropleth trace — they must match.
"""
function _get_site_ids(gdf::DataFrame)::Vector{String}
    if hasproperty(gdf, :site_id)
        return string.(gdf.site_id)
    elseif hasproperty(gdf, :reef_siteid)
        return string.(gdf.reef_siteid)
    else
        return string.(1:nrow(gdf))
    end
end

"""
    _geom_to_geojson_dict(geom) -> Dict{String,Any}

Convert an ArchGDAL Polygon or MultiPolygon to a GeoJSON geometry Dict.
Uses `GeoInterface.coordinates`, which returns the same nested-array structure
that GeoJSON expects. No JSON parsing dependency is introduced.
"""
function _geom_to_geojson_dict(geom)::Dict{String,Any}
    trait = GI.geomtrait(geom)
    type_str = if trait isa GI.MultiPolygonTrait
        "MultiPolygon"
    elseif trait isa GI.PolygonTrait
        "Polygon"
    else
        @warn "Unexpected geometry trait $(trait); treating as Polygon"
        "Polygon"
    end
    return Dict{String,Any}("type" => type_str, "coordinates" => GI.coordinates(geom))
end

"""
    _gdf_to_geojson(gdf::DataFrame; geom_col=nothing) -> Dict{String,Any}

Build a Plotly-compatible GeoJSON FeatureCollection Dict from a GeoDataFrame
(backed by a GeoPackage or any ArchGDAL source). Each row becomes a Feature
whose `id` matches the site identifier string used in the `locations` array
of the corresponding choropleth trace.
"""
function _gdf_to_geojson(
    gdf::DataFrame;
    geom_col::Union{Symbol,Nothing}=nothing
)::Dict{String,Any}
    col = isnothing(geom_col) ? ADRIA._get_geom_col(gdf) : geom_col
    col == false && throw(ArgumentError("No geometry column found in GeoDataFrame"))
    ids = _get_site_ids(gdf)
    features = [
        Dict{String,Any}(
            "type" => "Feature",
            "id" => ids[i],
            "geometry" => _geom_to_geojson_dict(gdf[i, col]),
            "properties" => Dict{String,Any}()
        )
        for i in 1:nrow(gdf)
    ]
    return Dict{String,Any}("type" => "FeatureCollection", "features" => features)
end

# =============================================================================
# Layout helpers
# =============================================================================

"""
    _map_geo_layout(; title="", width=700, height=900) -> Layout

Base Layout for all geographic maps. `fitbounds="locations"` auto-zooms to
the GeoJSON feature extent. No Mapbox token is required.
"""
function _map_geo_layout(; title::String="", width::Int=700, height::Int=900)
    return Layout(;
        font=attr(; family="Open Sans, sans-serif", size=12),
        paper_bgcolor="white",
        title_text=title,
        width=width,
        height=height,
        geo=attr(;
            showframe=false,
            showcoastlines=true,
            coastlinecolor="#aaaaaa",
            showland=true,
            landcolor="#f5f5f5",
            showocean=true,
            oceancolor="#e0eef5",
            fitbounds="locations",
            projection_type="mercator"
        ),
        margin=attr(; l=0, r=0, t=40, b=0)
    )
end

"""
    _geo_domain(row, col, n_rows, n_cols; margin=0.02) -> (x0, x1, y0, y1)

Compute the [0,1]-normalised domain rectangle for one panel in a multi-geo
subplot grid. `row` and `col` are 1-indexed.
"""
function _geo_domain(
    row::Int, col::Int, n_rows::Int, n_cols::Int; margin::Float64=0.02
)
    w = (1.0 - (n_cols + 1) * margin) / n_cols
    h = (1.0 - (n_rows + 1) * margin) / n_rows
    x0 = margin + (col - 1) * (w + margin)
    x1 = x0 + w
    y1 = 1.0 - margin - (row - 1) * (h + margin)
    y0 = y1 - h
    return x0, x1, y0, y1
end

# =============================================================================
# ADRIA.viz.map — plain GeoDataFrame
# =============================================================================

"""
    ADRIA.viz.map(gdf::DataFrame; color=nothing, colorbar_label="",
                  colorscale="Viridis", title="", width=700, height=900)

Plot a GeoDataFrame (backed by a GeoPackage). When `color` is `nothing`, all
features are rendered with a uniform grey fill showing reef outlines. Pass a
numeric vector (`length == nrow(gdf)`) to produce a choropleth coloured by
value.

No Mapbox token is required — Plotly's `choropleth` trace is used with a
custom GeoJSON FeatureCollection derived from the GeoPackage geometries.
`fitbounds="locations"` auto-zooms to the reef-site extent.

# Arguments
- `color`          : Per-site numeric values, or `nothing` for outline-only
- `colorbar_label` : Colorbar title text
- `colorscale`     : Plotly colorscale name (default `"Viridis"`)
- `title`          : Figure title
- `width`, `height`: Figure dimensions in pixels
"""
function ADRIA.viz.map(
    gdf::DataFrame;
    color::Union{AbstractVector{<:Real},Nothing}=nothing,
    colorbar_label::String="",
    colorscale::String="Viridis",
    title::String="",
    width::Int=700,
    height::Int=900,
    kwargs...
)
    col = ADRIA._get_geom_col(gdf)
    col == false && throw(ArgumentError("No geometry column found in GeoDataFrame"))
    id_vec = _get_site_ids(gdf)
    geojson = _gdf_to_geojson(gdf; geom_col=col)
    n = nrow(gdf)

    trace = if isnothing(color)
        choropleth(;
            geojson=geojson,
            locations=id_vec,
            z=ones(n),
            colorscale=[[0, "#cccccc"], [1, "#cccccc"]],
            showscale=false,
            marker_line_color="#555555",
            marker_line_width=0.5,
            type="choropleth",
            name=""
        )
    else
        choropleth(;
            geojson=geojson,
            locations=id_vec,
            z=collect(Float64, color),
            colorscale=colorscale,
            colorbar=attr(; title_text=colorbar_label, thickness=15),
            marker_line_color="#555555",
            marker_line_width=0.5,
            type="choropleth"
        )
    end

    return Plot(trace, _map_geo_layout(; title=title, width=width, height=height))
end

# =============================================================================
# ADRIA.viz.map — ResultSet / Domain with outcome values
# =============================================================================

"""
    ADRIA.viz.map(rs::Union{Domain,ResultSet},
                  y::Union{YAXArray,AbstractVector{<:Real}};
                  diverging=false, colorbar_label="", colorscale=nothing,
                  title="Study Area", width=700, height=900)

Choropleth map of a per-location outcome vector `y` over the reef sites in
`rs`. `y` must have `length == nrow(rs.loc_data)`.

Pass `diverging=true` to use the `"RdBu"` colorscale (suitable for
intervention-vs-counterfactual difference metrics).
"""
function ADRIA.viz.map(
    rs::Union{Domain,ResultSet},
    y::Union{YAXArray,AbstractVector{<:Real}};
    diverging::Bool=false,
    colorbar_label::String="",
    colorscale::Union{String,Nothing}=nothing,
    title::String="Study Area",
    width::Int=700,
    height::Int=900,
    kwargs...
)
    if isempty(colorbar_label) && y isa YAXArray
        colorbar_label = ADRIAviz.outcome_label(y; metadata_key=:metric_name)
    end
    cs = if !isnothing(colorscale)
        colorscale
    elseif diverging
        "RdBu"
    else
        "Viridis"
    end
    return ADRIA.viz.map(
        rs.loc_data;
        color=collect(Float64, y),
        colorbar_label=colorbar_label,
        colorscale=cs,
        title=title,
        width=width,
        height=height
    )
end

"""
    ADRIA.viz.map(rs::Union{Domain,ResultSet}; ...)

Plot the reef-site carrying capacity (`k * 100 %`) as the default metric.
"""
function ADRIA.viz.map(
    rs::Union{Domain,ResultSet};
    title::String="Study Area",
    width::Int=700,
    height::Int=900,
    kwargs...
)
    return ADRIA.viz.map(
        rs,
        rs.loc_data.k * 100.0;
        colorbar_label="Coral Real Estate [%]",
        title=title,
        width=width,
        height=height
    )
end

# =============================================================================
# ADRIA.viz.map — multi-panel matrix form
# =============================================================================

"""
    ADRIA.viz.map(rs::Union{Domain,ResultSet}, outputs_matrix::AbstractMatrix,
                  map_titles::Vector{String}; ...)

Render a grid of choropleth maps, one per column of `outputs_matrix`.
`outputs_matrix` must have `nrow(rs.loc_data)` rows and
`length(map_titles)` columns. Color range is shared across all panels.

# Arguments
- `colorscale`     : Plotly colorscale name (default `"Viridis"`)
- `colorbar_label` : Colorbar title text
- `title`          : Overall figure title
- `width`, `height`: Figure dimensions in pixels
"""
function ADRIA.viz.map(
    rs::Union{Domain,ResultSet},
    outputs_matrix::AbstractMatrix,
    map_titles::Vector{String};
    colorscale::String="Viridis",
    colorbar_label::String="",
    title::String="",
    width::Int=1200,
    height::Int=600,
    kwargs...
)
    n_plots = length(map_titles)
    n_rows, n_cols = _calc_gridsize(n_plots)
    id_vec = _get_site_ids(rs.loc_data)
    geojson = _gdf_to_geojson(rs.loc_data)

    # Shared color range across all panels
    cmin, cmax = extrema(filter(isfinite, vec(Float64.(outputs_matrix))))

    traces = AbstractTrace[]
    layout_attrs = Dict{Symbol,Any}(
        :title_text => title,
        :width => width,
        :height => height,
        :font => attr(; family="Open Sans, sans-serif", size=12),
        :paper_bgcolor => "white",
        :annotations => PlotlyBase.PlotlyAttribute[]
    )

    for k in 1:n_plots
        geo_ref = k == 1 ? "geo" : "geo$(k)"
        row = div(k - 1, n_cols) + 1
        col = mod(k - 1, n_cols) + 1
        x0, x1, y0, y1 = _geo_domain(row, col, n_rows, n_cols)

        t = choropleth(;
            geojson=geojson,
            locations=id_vec,
            z=vec(Float64.(outputs_matrix[:, k])),
            zmin=cmin,
            zmax=cmax,
            colorscale=colorscale,
            colorbar=attr(; title_text=colorbar_label, thickness=12, x=1.0),
            showscale=(k == n_plots),
            marker_line_color="#555555",
            marker_line_width=0.3,
            geo=geo_ref,
            type="choropleth",
            name=map_titles[k]
        )
        push!(traces, t)

        geo_attr_key = k == 1 ? :geo : Symbol("geo$(k)")
        layout_attrs[geo_attr_key] = attr(;
            domain=attr(; x=[x0, x1], y=[y0, y1]),
            showframe=false,
            showcoastlines=true,
            coastlinecolor="#aaaaaa",
            showland=true,
            landcolor="#f5f5f5",
            fitbounds="locations",
            projection_type="mercator"
        )

        # Panel title as annotation
        push!(
            layout_attrs[:annotations],
            attr(;
                text=map_titles[k],
                x=(x0 + x1) / 2,
                y=y1 + 0.01,
                xref="paper",
                yref="paper",
                xanchor="center",
                yanchor="bottom",
                showarrow=false,
                font=attr(; size=11)
            )
        )
    end

    return Plot(traces, Layout(; layout_attrs...))
end

# =============================================================================
# ADRIA.viz.connectivity
# =============================================================================

"""
    ADRIA.viz.connectivity(dom::Domain;
                           in_method=nothing,
                           out_method=eigenvector_centrality,
                           title="Connectivity", width=700, height=900)

Connectivity graph of reef sites. Reef polygons are drawn as a grey base map;
edges are semi-transparent lines between site centroids; nodes are sized and
coloured by connectivity weight. No Mapbox token is required.

See also `ADRIA.connectivity_strength` for the underlying graph computation.
"""
function ADRIA.viz.connectivity(
    dom::Domain;
    in_method=nothing,
    out_method=eigenvector_centrality,
    title::String="Connectivity",
    width::Int=700,
    height::Int=900,
    kwargs...
)
    return ADRIA.viz.connectivity(
        dom, dom.conn; in_method, out_method, title, width, height
    )
end

function ADRIA.viz.connectivity(
    dom::Domain,
    conn::AbstractMatrix;
    in_method=nothing,
    out_method=eigenvector_centrality,
    title::String="Connectivity",
    width::Int=700,
    height::Int=900,
    kwargs...
)
    if !isnothing(in_method) && !isnothing(out_method)
        @warn "Both in and out centrality measures provided. Plotting out centralities."
        cs = ADRIA.connectivity_strength(conn; in_method, out_method)
    elseif !isnothing(in_method) && isnothing(out_method)
        cs = ADRIA.connectivity_strength(
            conn; in_method, out_method=outdegree_centrality
        )
    elseif isnothing(in_method) && isnothing(out_method)
        error("At least one of in_method or out_method must be provided.")
    else
        cs = ADRIA.connectivity_strength(
            conn;
            in_method=isnothing(in_method) ? indegree_centrality : in_method,
            out_method=isnothing(out_method) ? outdegree_centrality : out_method
        )
    end
    return ADRIA.viz.connectivity(dom, cs.network, cs.out_conn; title, width, height)
end

"""
    ADRIA.viz.connectivity(dom::Domain, network::SimpleWeightedDiGraph,
                           conn_weights::AbstractVector{<:Real}; ...)

Low-level connectivity plot. `network` and `conn_weights` are as returned by
`ADRIA.connectivity_strength(dom.conn; ...)`.

Traces produced:
1. `choropleth` — reef polygon base map (grey, no colorbar)
2. `scattergeo` — edges as None-separated line segments (single trace)
3. `scattergeo` — nodes sized and coloured by centrality weight
"""
function ADRIA.viz.connectivity(
    dom::Domain,
    network::SimpleWeightedDiGraph,
    conn_weights::AbstractVector{<:Real};
    title::String="Connectivity",
    width::Int=700,
    height::Int=900,
    kwargs...
)
    cents = ADRIA.centroids(dom)
    n_locs = length(cents)
    id_vec = _get_site_ids(dom.loc_data)
    geojson = _gdf_to_geojson(dom.loc_data)

    # ── Base map (reef polygons, uniform grey) ────────────────────────────────
    base_trace = choropleth(;
        geojson=geojson,
        locations=id_vec,
        z=ones(n_locs),
        colorscale=[[0, "#dddddd"], [1, "#dddddd"]],
        showscale=false,
        marker_line_color="#888888",
        marker_line_width=0.5,
        type="choropleth",
        showlegend=false,
        name=""
    )

    # ── Edge trace (None-separated segments in a single scattergeo) ───────────
    norm_coef = max(maximum(conn_weights), 1e-9)
    edge_lons = Union{Float64,Nothing}[]
    edge_lats = Union{Float64,Nothing}[]
    for e in edges(network)
        e.src == e.dst && continue
        conn_weights[e.src] * e.weight / norm_coef < 0.01 && continue
        lon_src, lat_src = cents[e.src]
        lon_dst, lat_dst = cents[e.dst]
        push!(edge_lons, lon_src, lon_dst, nothing)
        push!(edge_lats, lat_src, lat_dst, nothing)
    end
    edge_trace = scattergeo(;
        lon=edge_lons,
        lat=edge_lats,
        mode="lines",
        line=attr(; color="rgba(30,30,30,0.25)", width=0.6),
        type="scattergeo",
        showlegend=false,
        name="Connectivity"
    )

    # ── Node trace (size + colour by centrality weight) ───────────────────────
    node_sizes = max.(3.0, collect(conn_weights) ./ norm_coef .* 14.0)
    node_lons = [c[1] for c in cents]
    node_lats = [c[2] for c in cents]
    node_trace = scattergeo(;
        lon=node_lons,
        lat=node_lats,
        mode="markers",
        marker=attr(;
            size=node_sizes,
            color=collect(conn_weights),
            colorscale="Viridis",
            colorbar=attr(; title_text="Centrality", thickness=12),
            line=attr(; color="white", width=0.5)
        ),
        type="scattergeo",
        name="Centrality",
        text=["Site $(i)" for i in 1:n_locs],
        hoverinfo="text"
    )

    layout = Layout(;
        font=attr(; family="Open Sans, sans-serif", size=12),
        paper_bgcolor="white",
        title_text=title,
        width=width,
        height=height,
        geo=attr(;
            showframe=false,
            showcoastlines=true,
            coastlinecolor="#aaaaaa",
            showland=true,
            landcolor="#f5f5f5",
            showocean=true,
            oceancolor="#e0eef5",
            fitbounds="locations",
            projection_type="mercator"
        ),
        margin=attr(; l=0, r=0, t=40, b=0),
        showlegend=false
    )

    return Plot([base_trace, edge_trace, node_trace], layout)
end
