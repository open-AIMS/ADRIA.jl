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

# Candidate location-id columns, only consulted when no domain context is
# available (see `_get_site_ids`). The authoritative source is always the
# domain's configured `loc_id_col` (see `_loc_id_col` / `_site_ids`).
const _SITE_ID_CANDIDATES = (:reef_siteid, :site_id, :UNIQUE_ID, :reef_id)

"""
    _loc_id_col(x::Union{Domain,ResultSet}) -> Symbol

The canonical per-location identifier column for the modelled spatial scale.
This is the column ADRIA uses to order `conn`/`centroids`, so it is the only
correct choice for matching choropleth `locations` to GeoJSON feature ids and
for labelling reef sites. Never guess column names when this is available.
"""
_loc_id_col(dom::Domain)::Symbol = Symbol(dom.loc_id_col)
_loc_id_col(rs::ResultSet)::Symbol = Symbol(rs.env_layer_md.loc_id_col)

"""
    _get_site_ids(gdf::DataFrame) -> Vector{String}

Fallback site-id resolver for when no `Domain`/`ResultSet` is available to
provide the canonical `loc_id_col`. Tries `_SITE_ID_CANDIDATES` in order, then
falls back to row indices. These values become both the GeoJSON feature `id`
fields and the `locations` array in a Plotly choropleth trace — they must match.

Because several id columns can coexist (e.g. a site-scale run carrying both
`site_id` and a reef-level `UNIQUE_ID`), this warns when the choice is
ambiguous: pass a `Domain`/`ResultSet` to resolve it deterministically.
"""
function _get_site_ids(gdf::DataFrame)::Vector{String}
    present = [c for c in _SITE_ID_CANDIDATES if hasproperty(gdf, c)]
    if isempty(present)
        return string.(1:nrow(gdf))
    end
    if length(present) > 1
        @warn "Multiple candidate location-id columns present $(present); using \
               '$(first(present))'. Pass a Domain/ResultSet so the canonical \
               loc_id_col is used instead."
    end
    return string.(gdf[:, first(present)])
end

"""
    _site_ids(gdf, id_col) -> Vector{String}

Resolve site ids using an explicit column when provided (the domain's
`loc_id_col`), otherwise fall back to the column-guessing `_get_site_ids`.
"""
function _site_ids(gdf::DataFrame, id_col::Union{Symbol,Nothing})::Vector{String}
    if !isnothing(id_col)
        if hasproperty(gdf, id_col)
            return string.(gdf[:, id_col])
        end
        @warn "Configured location-id column '$(id_col)' not found in \
               GeoDataFrame; falling back to heuristic resolution."
    end
    return _get_site_ids(gdf)
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
    geom_col::Union{Symbol,Nothing}=nothing,
    id_col::Union{Symbol,Nothing}=nothing
)::Dict{String,Any}
    col = isnothing(geom_col) ? ADRIA._get_geom_col(gdf) : geom_col
    col == false && throw(ArgumentError("No geometry column found in GeoDataFrame"))
    ids = _site_ids(gdf, id_col)
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
    _map_geo_layout(; title="", width=700, height=900,
                    lon_range=nothing, lat_range=nothing, annotations=nothing) -> Layout

Base Layout for all geographic maps. Uses Plotly's 50 m native coastline
(`resolution=50`, the highest the vector basemap supports) so near-shore reefs
are not drawn over a coarsely-generalised coastline. No Mapbox token is required.

When `lon_range`/`lat_range` are supplied the map is fixed to that extent
(used by the decorated single-reef maps so scale bar and place labels sit
sensibly); otherwise `fitbounds="locations"` auto-zooms to the feature extent.
"""
function _map_geo_layout(;
    title::String="",
    width::Int=700,
    height::Int=900,
    lon_range::Union{AbstractVector,Nothing}=nothing,
    lat_range::Union{AbstractVector,Nothing}=nothing,
    annotations::Union{AbstractVector,Nothing}=nothing
)
    geo = if isnothing(lon_range) || isnothing(lat_range)
        attr(;
            showframe=false, showcoastlines=true, coastlinecolor="#888888",
            resolution=50, showland=true, landcolor="#f5f5f5",
            showocean=true, oceancolor="#e0eef5",
            fitbounds="locations", projection_type="mercator"
        )
    else
        attr(;
            showframe=false, showcoastlines=true, coastlinecolor="#888888",
            resolution=50, showland=true, landcolor="#f5f5f5",
            showocean=true, oceancolor="#e0eef5",
            lonaxis=attr(; range=collect(Float64, lon_range)),
            lataxis=attr(; range=collect(Float64, lat_range)),
            projection_type="mercator"
        )
    end

    layout = Layout(;
        font=attr(; family="Open Sans, sans-serif", size=12),
        paper_bgcolor="white",
        title_text=title,
        width=width,
        height=height,
        geo=geo,
        margin=attr(; l=0, r=0, t=40, b=0)
    )
    isnothing(annotations) || (layout[:annotations] = annotations)
    return layout
end

# =============================================================================
# Map decorations: high-res context (Issue 4)
#   - scale bar sized to the reef extent
#   - paper-anchored north arrow
#   - labelled markers for major coastal centres within `max_km` of the reefs
# =============================================================================

# Major Queensland coastal population centres (name, lon, lat). Used to give
# spatial context; filtered to those within `max_km` of the modelled reefs, so
# non-GBR domains simply show none.
const GBR_COASTAL_PLACES = [
    ("Cooktown", 145.251, -15.468), ("Port Douglas", 145.465, -16.483),
    ("Cairns", 145.770, -16.920), ("Innisfail", 146.030, -17.524),
    ("Cardwell", 146.020, -18.267), ("Ingham", 146.158, -18.650),
    ("Townsville", 146.816, -19.258), ("Ayr", 147.405, -19.574),
    ("Bowen", 148.246, -20.013), ("Airlie Beach", 148.718, -20.267),
    ("Mackay", 149.186, -21.144), ("Yeppoon", 150.742, -23.130),
    ("Rockhampton", 150.510, -23.378), ("Gladstone", 151.256, -23.843),
    ("Bundaberg", 152.349, -24.866), ("Hervey Bay", 152.844, -25.290),
    ("Brisbane", 153.025, -27.470)
]

"""Great-circle distance between two lon/lat points, in kilometres."""
function _haversine_km(lon1, lat1, lon2, lat2)::Float64
    R = 6371.0
    dlat = deg2rad(lat2 - lat1)
    dlon = deg2rad(lon2 - lon1)
    a = sin(dlat / 2)^2 + cosd(lat1) * cosd(lat2) * sin(dlon / 2)^2
    return 2R * asin(min(1.0, sqrt(a)))
end

"""Round a distance down to a cartographically 'nice' scale-bar length (km)."""
function _nice_length(x::Real)::Int
    options = [1, 2, 5, 10, 20, 25, 50, 100, 200, 500, 1000]
    idx = findlast(<=(x), options)
    return options[isnothing(idx) ? 1 : idx]
end

"""
    _map_decorations(gdf; max_km=100.0)
        -> (traces, annotations, lon_range, lat_range)

Build scale-bar / place-name scattergeo traces and a north-arrow annotation for
a reef map, plus the lon/lat display range expanded to include any qualifying
coastal centres. Returns empty decorations if centroids cannot be computed.
"""
function _map_decorations(gdf::DataFrame; max_km::Float64=100.0)
    cents = ADRIA.centroids(gdf)
    isempty(cents) &&
        return (AbstractTrace[], PlotlyBase.PlotlyAttribute[], nothing, nothing)

    lons = Float64[c[1] for c in cents]
    lats = Float64[c[2] for c in cents]
    lon_min, lon_max = extrema(lons)
    lat_min, lat_max = extrema(lats)

    # Coastal centres within max_km of any reef
    places = NamedTuple{(:name, :lon, :lat),Tuple{String,Float64,Float64}}[]
    for (name, plon, plat) in GBR_COASTAL_PLACES
        d = minimum(_haversine_km(plon, plat, lons[i], lats[i]) for i in eachindex(lons))
        d <= max_km && push!(places, (name=name, lon=plon, lat=plat))
    end

    traces = AbstractTrace[]
    if !isempty(places)
        push!(
            traces,
            scattergeo(;
                lon=[p.lon for p in places], lat=[p.lat for p in places],
                text=[p.name for p in places], mode="markers+text",
                marker=attr(; size=6, color="#333333", symbol="circle"),
                # Label sits to the upper-left of the marker.
                textposition="top left", textfont=attr(; size=11, color="#333333"),
                hoverinfo="text", showlegend=false, name="Places"
            )
        )
    end

    # Display range = reefs + qualifying places, with a small pad
    all_lons = vcat(lons, [p.lon for p in places])
    all_lats = vcat(lats, [p.lat for p in places])
    lo_min, lo_max = extrema(all_lons)
    la_min, la_max = extrema(all_lats)
    lon_pad = max(0.05, 0.08 * (lo_max - lo_min))
    lat_pad = max(0.05, 0.08 * (la_max - la_min))
    lon_range = [lo_min - lon_pad, lo_max + lon_pad]
    lat_range = [la_min - lat_pad, la_max + lat_pad]

    # Scale bar (~1/4 of the reef extent width, rounded to a nice length)
    lat_mid = (lat_min + lat_max) / 2
    width_km = _haversine_km(lon_min, lat_mid, lon_max, lat_mid)
    bar_km = _nice_length(max(width_km / 4, 1.0))
    bar_deg = bar_km / (111.32 * cosd(lat_mid))
    bx0 = lon_range[1] + 0.04 * (lon_range[2] - lon_range[1])
    by = lat_range[1] + 0.06 * (lat_range[2] - lat_range[1])
    push!(
        traces,
        scattergeo(;
            lon=[bx0, bx0 + bar_deg], lat=[by, by], mode="lines",
            line=attr(; color="black", width=3), hoverinfo="skip",
            showlegend=false, name="scale"
        )
    )
    push!(
        traces,
        scattergeo(;
            lon=[bx0 + bar_deg / 2], lat=[by], text=["$(bar_km) km"], mode="text",
            textposition="top center", textfont=attr(; size=11, color="black"),
            hoverinfo="skip", showlegend=false, name="scale_label"
        )
    )

    # North arrow (paper-anchored; north is up in mercator)
    north = attr(;
        text="N", x=0.07, y=0.92, xref="paper", yref="paper",
        showarrow=true, ax=0, ay=28, arrowhead=2, arrowsize=1.3, arrowwidth=2,
        arrowcolor="black", font=attr(; size=15, color="black"),
        xanchor="center", yanchor="bottom"
    )

    return traces, PlotlyBase.PlotlyAttribute[north], lon_range, lat_range
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
    id_col::Union{Symbol,Nothing}=nothing,
    decorate::Bool=false,
    kwargs...
)
    col = ADRIA._get_geom_col(gdf)
    col == false && throw(ArgumentError("No geometry column found in GeoDataFrame"))
    id_vec = _site_ids(gdf, id_col)
    geojson = _gdf_to_geojson(gdf; geom_col=col, id_col=id_col)
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

    # Context decorations (scale bar, north arrow, coastal place names).
    # Guarded: a decoration failure must never break the underlying map.
    deco_traces, annotations, lon_range, lat_range = if decorate
        try
            _map_decorations(gdf)
        catch err
            @warn "Map decoration failed; rendering plain map." exception = err
            (AbstractTrace[], nothing, nothing, nothing)
        end
    else
        (AbstractTrace[], nothing, nothing, nothing)
    end

    layout = _map_geo_layout(;
        title=title, width=width, height=height,
        lon_range=lon_range, lat_range=lat_range,
        annotations=(
            isnothing(annotations) || isempty(annotations) ? nothing : annotations
        )
    )

    return Plot(AbstractTrace[trace; deco_traces...], layout)
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
        height=height,
        id_col=_loc_id_col(rs),
        decorate=true
    )
end

"""
    ADRIA.viz.map(rs::Union{Domain,ResultSet},
                  y::Union{YAXArray,AbstractVector{<:Real}},
                  clusters::AbstractVector{<:Integer};
                  colorbar_label="", title="Study Area", width=1200, height=900)

Two-panel choropleth map: left panel shows per-location metric `y` coloured by
value; right panel shows cluster membership as discrete categorical colours.
Plotly choropleth does not support per-location stroke colours, so cluster
membership is shown as a separate panel rather than as outline highlights.
"""
function ADRIA.viz.map(
    rs::Union{Domain,ResultSet},
    y::Union{YAXArray,AbstractVector{<:Real}},
    clusters::AbstractVector{<:Integer};
    colorbar_label::String="",
    title::String="Study Area",
    width::Int=1200,
    height::Int=900,
    kwargs...
)
    if isempty(colorbar_label) && y isa YAXArray
        colorbar_label = ADRIAviz.outcome_label(y; metadata_key=:metric_name)
    end

    id_col = _loc_id_col(rs)
    id_vec = _site_ids(rs.loc_data, id_col)
    geojson = _gdf_to_geojson(rs.loc_data; id_col=id_col)

    metric_vals = collect(Float64, y)
    cluster_ids = collect(Int, clusters)
    cluster_levels = sort(unique(cluster_ids))
    n_clusters = length(cluster_levels)

    # Distinct qualitative colors for up to ~10 clusters.
    _qual = [
        "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
        "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf"
    ]
    cluster_colorscale = [
        [n_clusters == 1 ? 0.0 : (i - 1) / (n_clusters - 1),
            _qual[mod1(i, length(_qual))]]
        for i in eachindex(cluster_levels)
    ]

    trace_metric = choropleth(;
        geojson=geojson,
        locations=id_vec,
        z=metric_vals,
        colorscale="Viridis",
        colorbar=attr(; title_text=colorbar_label, thickness=12, x=0.46),
        marker_line_color="#555555",
        marker_line_width=0.3,
        geo="geo",
        type="choropleth",
        name="Metric"
    )

    trace_cluster = choropleth(;
        geojson=geojson,
        locations=id_vec,
        z=Float64.(cluster_ids),
        zmin=Float64(minimum(cluster_ids)),
        zmax=Float64(maximum(cluster_ids)),
        colorscale=cluster_colorscale,
        colorbar=attr(; title_text="Cluster", thickness=12, x=1.0,
            tickvals=Float64.(cluster_levels),
            ticktext=string.(cluster_levels)),
        marker_line_color="#555555",
        marker_line_width=0.3,
        geo="geo2",
        type="choropleth",
        name="Cluster"
    )

    geo_common = attr(;
        showframe=false,
        showcoastlines=true,
        coastlinecolor="#888888",
        resolution=50,
        showland=true,
        landcolor="#f5f5f5",
        fitbounds="locations",
        projection_type="mercator"
    )

    layout = Layout(;
        title_text=title,
        width=width,
        height=height,
        font=attr(; family="Open Sans, sans-serif", size=12),
        paper_bgcolor="white",
        geo=merge(geo_common, attr(; domain=attr(; x=[0.0, 0.45], y=[0.0, 1.0]))),
        geo2=merge(geo_common, attr(; domain=attr(; x=[0.55, 1.0], y=[0.0, 1.0]))),
        annotations=[
            attr(; text="Metric", x=0.225, y=1.02, xref="paper", yref="paper",
                showarrow=false, font=attr(; size=13)),
            attr(; text="Cluster", x=0.775, y=1.02, xref="paper", yref="paper",
                showarrow=false, font=attr(; size=13))
        ]
    )

    return PlotlyBase.Plot([trace_metric, trace_cluster], layout)
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
    id_col = _loc_id_col(rs)
    id_vec = _site_ids(rs.loc_data, id_col)
    geojson = _gdf_to_geojson(rs.loc_data; id_col=id_col)

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
            coastlinecolor="#888888",
            resolution=50,
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
    min_edge_weight::Real=0.01,
    kwargs...
)
    cents = ADRIA.centroids(dom)
    n_locs = length(cents)
    id_col = _loc_id_col(dom)
    id_vec = _site_ids(dom.loc_data, id_col)
    geojson = _gdf_to_geojson(dom.loc_data; id_col=id_col)

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

    # ── Edge traces ───────────────────────────────────────────────────────────
    # Normalise against the maximum *edge* weight (not node centrality, which is
    # a different scale — mixing them previously filtered out nearly every edge).
    # Edges are grouped into weight classes so line width/opacity can scale with
    # connection strength (a single None-separated trace cannot vary per-segment).
    edge_list = collect(edges(network))
    max_w = max(maximum((e.weight for e in edge_list); init=0.0), 1e-9)

    n_bins = 3
    bin_lons = [Union{Float64,Nothing}[] for _ in 1:n_bins]
    bin_lats = [Union{Float64,Nothing}[] for _ in 1:n_bins]
    for e in edge_list
        e.src == e.dst && continue
        w = e.weight / max_w
        w < min_edge_weight && continue
        b = clamp(ceil(Int, w * n_bins), 1, n_bins)
        lon_src, lat_src = cents[e.src]
        lon_dst, lat_dst = cents[e.dst]
        push!(bin_lons[b], lon_src, lon_dst, nothing)
        push!(bin_lats[b], lat_src, lat_dst, nothing)
    end

    edge_widths = (0.5, 1.2, 2.2)
    edge_opacities = (0.25, 0.45, 0.7)
    edge_traces = AbstractTrace[]
    for b in 1:n_bins
        isempty(bin_lons[b]) && continue
        push!(
            edge_traces,
            scattergeo(;
                lon=bin_lons[b],
                lat=bin_lats[b],
                mode="lines",
                line=attr(;
                    color="rgba(30,30,30,$(edge_opacities[b]))", width=edge_widths[b]
                ),
                type="scattergeo",
                showlegend=false,
                hoverinfo="skip",
                name="Connectivity"
            )
        )
    end

    # ── Node trace (size + colour by centrality weight) ───────────────────────
    norm_coef = max(maximum(conn_weights), 1e-9)
    node_sizes = max.(4.0, collect(conn_weights) ./ norm_coef .* 16.0)
    node_lons = [c[1] for c in cents]
    node_lats = [c[2] for c in cents]
    node_text = [
        "$(id_vec[i])<br>Centrality: $(round(conn_weights[i]; digits=3))"
        for i in 1:n_locs
    ]
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
        text=node_text,
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
            coastlinecolor="#888888",
            resolution=50,
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

    return Plot(AbstractTrace[base_trace; edge_traces...; node_trace], layout)
end
