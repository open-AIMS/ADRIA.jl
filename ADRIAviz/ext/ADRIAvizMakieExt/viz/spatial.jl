import ArchGDAL as AG
using Graphs, GraphMakie, SimpleWeightedGraphs
using DataFrames: DataFrame

using ADRIA: _get_geom_col, Domain, ResultSet

# Import shared spatial utilities
using ADRIAviz: _loc_id_col, compute_map_decorations, MapDecorationData, validate_extent

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

# =============================================================================
# Map decorations: scale bar, north arrow, coastal places
# =============================================================================

"""
    _render_coastal_places!(ax::GeoAxis, deco::MapDecorationData)

Render coastal place labels on a Makie GeoAxis from decoration data.
"""
function _render_coastal_places!(
    ax::GeoAxis,
    deco::MapDecorationData;
    decoration_fontsize::Int=9
)
    if isempty(deco.places)
        return nothing
    end
    scatter!(
        ax,
        [p.lon for p in deco.places],
        [p.lat for p in deco.places];
        color=:black,
        markersize=4,
        label="Places"
    )
    for place in deco.places
        text!(
            ax,
            place.lon,
            place.lat;
            text=place.name,
            fontsize=decoration_fontsize,
            align=(:left, :top),
            offset=(-5, 5)
        )
    end
end

"""
    _render_scale_bar!(ax::GeoAxis, deco::MapDecorationData)

Render a scale bar on a Makie GeoAxis from decoration data.
"""
function _render_scale_bar!(
    ax::GeoAxis,
    deco::MapDecorationData;
    decoration_fontsize::Int=9
)
    linesegments!(
        ax,
        [deco.scale_bar_x0, deco.scale_bar_x0 + deco.scale_bar_deg],
        [deco.scale_bar_y, deco.scale_bar_y];
        color=:black,
        linewidth=2
    )
    text!(
        ax,
        deco.scale_bar_x0 + deco.scale_bar_deg / 2,
        deco.scale_bar_y + 0.5;
        text="$(deco.scale_bar_km) km",
        fontsize=decoration_fontsize,
        align=(:center, :bottom)
    )
end

"""
    _render_coastlines!(ax::GeoAxis, resolution::Int=50)

Render Natural Earth coastlines onto a GeoAxis. `resolution` is one of 110, 50, or 10.
50m and 10m data require NaturalEarth.jl (`pkg> add NaturalEarth`).
110m data is bundled with GeoMakie and requires no download.
"""
function _render_coastlines!(ax::GeoAxis, resolution::Int=50)
    resolution in (110, 50, 10) || error(
        "Coastline resolution must be 10, 50, or 110 (got $resolution)"
    )
    ne_pkgid = Base.PkgId(
        Base.UUID("436b0209-26ab-4e65-94a9-6526d86fea76"), "NaturalEarth"
    )
    isnothing(Base.locate_package(ne_pkgid)) && error(
        "NaturalEarth.jl is required for $(resolution)m coastlines.\n" *
        "Install it with:  pkg> add NaturalEarth"
    )
    coastline_data = GeoMakie.coastlines(resolution)

    lines!(ax, coastline_data; color=(:gray40, 0.6), linewidth=0.5)
    return nothing
end

"""
    _render_north_arrow!(ax::GeoAxis, deco::MapDecorationData)

Render a minimal north arrow in the upper-right corner of the map extent.
"""
function _render_north_arrow!(
    ax::GeoAxis,
    deco::MapDecorationData;
    north_arrow_fontsize::Int=20
)
    lon_span = deco.lon_range[2] - deco.lon_range[1]
    lat_span = deco.lat_range[2] - deco.lat_range[1]
    x = deco.scale_bar_x0 + deco.scale_bar_deg + 0.025 * lon_span
    y = deco.scale_bar_y
    text!(ax, x, y; text="N", fontsize=north_arrow_fontsize, align=(:center, :bottom))
    text!(
        ax,
        x,
        y;
        text="\u2191",
        fontsize=north_arrow_fontsize,
        align=(:center, :top)
    )
    return nothing
end

"""
    _render_map_decorations!(ax::GeoAxis, gdf::DataFrame; kwargs...)

Render map decorations (places, scale bar, north arrow) on a Makie GeoAxis.
Fails silently if decorations cannot be computed.
"""
function _render_map_decorations!(
    ax::GeoAxis,
    gdf::DataFrame;
    min_population::Int=50000,
    max_km::Float64=150.0,
    show_coastal_places::Bool=true,
    show_scale_bar::Bool=true,
    show_north_arrow::Bool=true,
    decoration_fontsize::Int=10,
    north_arrow_fontsize::Int=20
)
    deco_data = compute_map_decorations(gdf; min_population=min_population, max_km=max_km)
    isnothing(deco_data) && return nothing

    show_coastal_places &&
        _render_coastal_places!(ax, deco_data; decoration_fontsize=decoration_fontsize)
    show_scale_bar &&
        _render_scale_bar!(ax, deco_data; decoration_fontsize=decoration_fontsize)
    show_north_arrow &&
        _render_north_arrow!(ax, deco_data; north_arrow_fontsize=north_arrow_fontsize)
    return nothing
end

"""
    _figure_size(n_rows::Int, n_cols::Int)::Tuple{Int,Int}

Compute a landscape-oriented figure size in pixels for a grid of `n_rows` x `n_cols` panels.
Single panel: 1000x800. Multi-panel: 600 wide x 500 tall per panel with 20 px gaps.
"""
function _figure_size(n_rows::Int, n_cols::Int)::Tuple{Int,Int}
    if n_rows == 1 && n_cols == 1
        return (1000, 800)
    end
    panel_w, panel_h, gap = 600, 400, 20
    width = n_cols * panel_w + (n_cols - 1) * gap
    height = n_rows * panel_h + (n_rows - 1) * gap
    return (width, height)
end

"""
    _adaptive_gap(n_panels::Int)::Tuple{Int,Int}

Compute adaptive row and column gaps (in pixels) for a multi-panel grid.
Gaps shrink as panel count increases to reclaim space without crowding.
"""
function _adaptive_gap(n_panels::Int)::Tuple{Int,Int}
    rowgap = max(10, 25 - div(n_panels, 4) * 5)
    colgap = max(8, 20 - div(n_panels, 3) * 4)
    return (rowgap, colgap)
end

function set_figure_defaults(fig_opts::OPT_TYPE; n_rows::Int=1, n_cols::Int=1)::OPT_TYPE
    w, h = _figure_size(n_rows, n_cols)
    w, h = max(w, 400), max(h, 300)
    fig_opts[:size] = get(fig_opts, :size, (w, h))

    actual_size = fig_opts[:size]
    default_padding = max(12, Int(round(0.03 * maximum(actual_size))))
    fig_opts[:figure_padding] = get(fig_opts, :figure_padding, default_padding)

    return fig_opts
end

function set_axis_defaults(axis_opts::OPT_TYPE; n_panels::Int=1)::OPT_TYPE
    axis_opts[:title] = get(axis_opts, :title, "Study Area")
    axis_opts[:xlabel] = get(axis_opts, :xlabel, "Longitude")
    axis_opts[:ylabel] = get(axis_opts, :ylabel, "Latitude")
    axis_opts[:xgridwidth] = get(axis_opts, :xgridwidth, 0.5)
    axis_opts[:ygridwidth] = get(axis_opts, :ygridwidth, 0.5)
    axis_opts[:dest] = get(axis_opts, :dest, "+proj=latlong +datum=WGS84")
    axis_opts[:xgridvisible] = get(axis_opts, :xgridvisible, false)
    axis_opts[:ygridvisible] = get(axis_opts, :ygridvisible, false)

    # Intentional baseline alignment with shared typography tiers: single-panel title 14->16 and 5+-panel tick 7->8.
    set_typography_defaults!(axis_opts; n_panels=n_panels)

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
        color_map::$COLORMAP_TYPE=:grayC,
        legend_params::Union{Tuple,Nothing}=nothing,
        axis_opts::Dict{Symbol, <:Any}=set_axis_defaults(Dict{Symbol,Any}());
        gdf=nothing,
        show_coastlines=true,
        coastline_resolution=50,
        show_coastal_places=true,
        show_scale_bar=true,
        show_north_arrow=true
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
- `gdf` : GeoDataFrame (required for coastal decorations)
- `show_coastlines` : Overlay Natural Earth coastlines; default `true` (110m bundled; 50/10m require NaturalEarth.jl)
- `coastline_resolution` : Coastline resolution in metres: 110, 50, or 10 (default 50)
- `show_coastal_places` : Overlay nearby Queensland coastal population centres (default `true`)
- `show_scale_bar` : Overlay a scale bar in the lower-left corner (default `true`)
- `show_north_arrow` : Overlay a north arrow in the upper-right corner (default `true`)
"""
function create_map!(
    f::Union{GridLayout,GridPosition},
    geodata::Vector{<:GeoMakie.MultiPolygon},
    data::Observable,
    highlight::Union{Vector,Tuple,Nothing};
    show_colorbar::Bool=true,
    colorbar_label::String="",
    color_map::COLORMAP_TYPE=:grayC,
    legend_params::Union{Tuple,Nothing}=nothing,
    axis_opts::OPT_TYPE=set_axis_defaults(DEFAULT_OPT_TYPE()),
    gdf::Union{DataFrame,Nothing}=nothing,
    show_coastlines::Bool=true,
    coastline_resolution::Int=110,
    show_coastal_places::Bool=true,
    show_scale_bar::Bool=true,
    show_north_arrow::Bool=true,
    min_population::Int=50000,
    max_km::Float64=150.0,
    decoration_fontsize::Int=9,
    north_arrow_fontsize::Int=12
)::Union{GridLayout,GridPosition}
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

    # Coastlines rendered before reefs so reefs appear on top
    show_coastlines && _render_coastlines!(spatial, coastline_resolution)

    poly!(
        spatial,
        geodata;
        color=data,
        colormap=color_map,
        colorrange=color_range,
        strokecolor=(:black, 0.05),
        strokewidth=1.0
    )

    # Overlays (places, scale bar, north arrow) rendered after reefs
    if !isnothing(gdf) && (show_coastal_places || show_scale_bar || show_north_arrow)
        deco_data = compute_map_decorations(
            gdf; min_population=min_population, max_km=max_km
        )
        if !isnothing(deco_data)
            # Set axis limits to the computed extent of reefs and nearby places
            limits!(
                spatial,
                deco_data.lon_range,
                deco_data.lat_range
            )

            _render_map_decorations!(
                spatial,
                gdf;
                min_population=min_population,
                max_km=max_km,
                show_coastal_places=show_coastal_places,
                show_scale_bar=show_scale_bar,
                show_north_arrow=show_north_arrow,
                decoration_fontsize=decoration_fontsize,
                north_arrow_fontsize=north_arrow_fontsize
            )
        end
    else
        # Even without decorations, we need to set axis limits to the data extent
        if !isnothing(gdf)
            deco_data = compute_map_decorations(
                gdf; min_population=min_population, max_km=max_km
            )
            if !isnothing(deco_data)
                limits!(
                    spatial,
                    deco_data.lon_range,
                    deco_data.lat_range
                )
            end
        end
    end

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
                isempty(m) && continue

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
- `rs` : ResultSet or Domain
- `y` : results of scenario metric
- `diverging` : If true, uses a diverging color map.
- `opts` : Aviz options
    - `colorbar_label` : label for colorbar. Defaults to "Relative Cover"
    - `color_map` : preferred colormap for plotting heatmaps
    - `show_coastlines` : overlay Natural Earth coastlines (default `true`).
      110m data is bundled with GeoMakie; 50m/10m require `pkg> add NaturalEarth`.
    - `coastline_resolution` : 110, 50, or 10 (default 50)
    - `show_coastal_places` : overlay nearby Queensland coastal centres (default `true`)
    - `show_scale_bar` : overlay a scale bar in the lower-left corner (default `true`)
    - `show_north_arrow` : overlay a north arrow in the upper-right corner (default `true`)
- `axis_opts` : Additional options to pass to adjust Axis attributes
  See: https://docs.makie.org/v0.19/api/index.html#Axis

# Returns
Figure
"""
function ADRIA.viz.map(
    rs::Union{Domain,ResultSet},
    y::Union{YAXArray,AbstractVector{<:Real}},
    highlight::AbstractVector;
    diverging::Bool=false,
    opts::OPT_TYPE=DEFAULT_OPT_TYPE(),
    fig_opts::OPT_TYPE=set_figure_defaults(DEFAULT_OPT_TYPE()),
    axis_opts::OPT_TYPE=set_axis_defaults(DEFAULT_OPT_TYPE())
)
    set_plot_opts!(y, opts, :colorbar_label; metadata_key=:metric_name)

    if diverging
        opts[:color_map] = _diverging_cmap(y)
    end

    opts[:highlight] = highlight

    f = Figure(; fig_opts...)
    g = f[1, 1] = GridLayout()

    ADRIA.viz.map!(g, rs, collect(y); opts, axis_opts)

    return f
end
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
    set_axis_defaults(axis_opts)

    highlight = get(opts, :highlight, nothing)
    c_label = get(opts, :colorbar_label, "")
    legend_params = get(opts, :legend_params, nothing)
    show_colorbar = get(opts, :show_colorbar, true)
    color_map = get(opts, :color_map, :grayC)
    show_coastlines = get(opts, :show_coastlines, true)
    coastline_resolution = get(opts, :coastline_resolution, 50)
    show_coastal_places = get(opts, :show_coastal_places, true)
    show_scale_bar = get(opts, :show_scale_bar, true)
    show_north_arrow = get(opts, :show_north_arrow, true)
    tick_sz = get(axis_opts, :yticklabelsize, 10)
    label_sz = get(axis_opts, :ylabelsize, 12)
    deco_fontsize = max(7, tick_sz - 1)
    north_arrow_fontsize = label_sz
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
        axis_opts=axis_opts,
        gdf=rs.loc_data,
        show_coastlines=show_coastlines,
        coastline_resolution=coastline_resolution,
        show_coastal_places=show_coastal_places,
        show_scale_bar=show_scale_bar,
        show_north_arrow=show_north_arrow,
        decoration_fontsize=deco_fontsize,
        north_arrow_fontsize=north_arrow_fontsize
    )
end

"""
    ADRIA.viz.map(rs::Union{Domain,ResultSet}, outputs_matrix::Matrix, map_titles::Vector{String}; opts::OPT_TYPE=DEFAULT_OPT_TYPE(), fig_opts::OPT_TYPE=DEFAULT_OPT_TYPE(), axis_opts::OPT_TYPE=DEFAULT_OPT_TYPE())
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
    opts::OPT_TYPE=DEFAULT_OPT_TYPE(),
    fig_opts::OPT_TYPE=DEFAULT_OPT_TYPE(),
    axis_opts::OPT_TYPE=DEFAULT_OPT_TYPE(),
    colorscale::Union{Symbol,String,Nothing}=nothing,
    colorbar_label::String="",
    title::String="",
    kwargs...
)
    n_rows, n_cols = _calc_gridsize(length(map_titles))
    set_figure_defaults(fig_opts; n_rows=n_rows, n_cols=n_cols)
    set_axis_defaults(axis_opts; n_panels=n_rows * n_cols)

    # Map Plotly-style parameters to Makie opts
    !isnothing(colorscale) && (opts[:color_map] = colorscale)
    !isempty(colorbar_label) && (opts[:colorbar_label] = colorbar_label)

    f = Figure(; fig_opts...)
    g = f[1, 1] = GridLayout()
    ADRIA.viz.map!(
        g, rs, outputs_matrix, map_titles; opts=opts, axis_opts=axis_opts
    )
    resize_to_layout!(f)
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
    if nrow(rs.loc_data) != size(outputs_matrix, 1)
        error("Only unfiltered decision matrices can be plotted.")
    end

    opts[:color_map] = get(opts, :color_map, :viridis)
    opts[:colorbar_limits] = get(opts, :colorbar_limits, (0.0, 1.0))

    n_plots::Int64 = length(map_titles)
    n_rows, n_cols = _calc_gridsize(n_plots)
    set_axis_defaults(axis_opts; n_panels=n_plots)

    # Suppress per-panel colorbars; a shared one is added after the loop
    opts[:show_colorbar] = false
    opts[:show_scale_bar] = get(opts, :show_scale_bar, true)
    opts[:show_north_arrow] = get(opts, :show_north_arrow, true)
    opts[:show_coastal_places] = get(opts, :show_coastal_places, true)

    step::Int64 = 1
    for row = 1:n_rows, col = 1:n_cols
        if step > length(map_titles)
            break
        end
        axis_opts[:title] = map_titles[step]
        # Hide axis labels on interior panels to reduce clutter
        axis_opts[:ylabelvisible] = col == 1
        axis_opts[:xlabelvisible] = row == n_rows
        axis_opts[:yticklabelsvisible] = col == 1
        axis_opts[:xticklabelsvisible] = row == n_rows

        # Show scale bar on all panels, but north arrow only on the first panel
        panel_opts = copy(opts)
        if step > 1
            panel_opts[:show_north_arrow] = false
        end

        ADRIA.viz.map!(
            g[row, col],
            rs,
            vec(outputs_matrix[:, step]);
            opts=panel_opts,
            axis_opts=axis_opts
        )

        step += 1
    end

    # Shared colorbar spanning full grid height
    Colorbar(
        g[1:n_rows, n_cols + 1];
        colorrange=opts[:colorbar_limits],
        colormap=opts[:color_map],
        label=get(opts, :colorbar_label, "")
    )

    # Adaptive gaps: shrink as panel count increases
    rg, cg = _adaptive_gap(n_plots)
    if g isa GridLayout
        rowgap!(g, rg)
        colgap!(g, cg)
    end

    return g
end

function _diverging_cmap(outcomes::YAXArray)::Vector{RGB}
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
function ADRIA.viz.map(
    gdf::DataFrame;
    geom_col=:geometry,
    color=nothing,
    title="",
    fig_opts::OPT_TYPE=DEFAULT_OPT_TYPE(),
    axis_opts::OPT_TYPE=DEFAULT_OPT_TYPE()
)
    fig_opts[:size] = get(fig_opts, :size, (600, 900))
    set_figure_defaults(fig_opts)
    set_axis_defaults(axis_opts)
    !isempty(title) && (axis_opts[:title] = title)

    f = Figure(; fig_opts...)
    ga = GeoAxis(
        f[1, 1];
        xticklabelpad=15,
        yticklabelpad=10,
        aspect=DataAspect(),
        axis_opts...
    )

    ADRIA.viz.map!(ga, gdf; geom_col=geom_col, color=color)

    if !isnothing(color)
        finite_vals = filter(isfinite, color)  # clear out any infs
        Colorbar(
            f[1, 2];
            colormap=:viridis,
            colorrange=(minimum(finite_vals), maximum(finite_vals)),
            vertical=true,
            height=Relative(0.8),
            tellwidth=true
        )
    end

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
dom = ADRIA.load_domain("<Path to Domain>", "<RCP>")
ADRIA.viz.connectivity(dom)

# Plot indegree centrality instead
ADRIA.viz.connectivity(dom; in_method=indegree_centrality, out_method=nothing)
```

Finer grain control:

```julia
dom = ADRIA.load_domain("<Path to Domain>", "<RCP>")

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
    - `edge_quantile`, quantile threshold for culling weak edges (default 0.2). Edges with alpha below this quantile are hidden to improve rendering performance. Set to 0.0 to show all edges.
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
    fig_opts::OPT_TYPE=DEFAULT_OPT_TYPE(),
    axis_opts::OPT_TYPE=DEFAULT_OPT_TYPE()
)
    set_figure_defaults(fig_opts)
    set_axis_defaults(axis_opts)
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
    fig_opts::OPT_TYPE=DEFAULT_OPT_TYPE(),
    axis_opts::OPT_TYPE=DEFAULT_OPT_TYPE()
)
    set_figure_defaults(fig_opts)
    set_axis_defaults(axis_opts)
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
    fig_opts::OPT_TYPE=DEFAULT_OPT_TYPE(),
    axis_opts::OPT_TYPE=DEFAULT_OPT_TYPE()
)
    set_figure_defaults(fig_opts)
    set_axis_defaults(axis_opts)
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
    axis_opts::OPT_TYPE=DEFAULT_OPT_TYPE()
)
    set_axis_defaults(axis_opts)

    show_coastlines = get(opts, :show_coastlines, true)
    coastline_resolution = get(opts, :coastline_resolution, 50)
    show_coastal_places = get(opts, :show_coastal_places, true)
    show_scale_bar = get(opts, :show_scale_bar, true)
    show_north_arrow = get(opts, :show_north_arrow, true)
    min_population = get(opts, :min_population, 50000)
    max_km = get(opts, :max_km, 150.0)
    tick_sz = get(axis_opts, :yticklabelsize, 10)
    label_sz = get(axis_opts, :ylabelsize, 12)
    decoration_fontsize = max(7, tick_sz - 1)
    north_arrow_fontsize = label_sz

    if !(:dest in keys(axis_opts))
        col = _get_geom_col(dom.loc_data)
        axis_opts[:dest] = ADRIA.AG.toPROJ4(ADRIA.AG.getspatialref(dom.loc_data[1, col]))
    end

    spatial = GeoAxis(
        g[1, 1];
        axis_opts...
    )

    spatial.yticklabelpad = 50
    spatial.ytickalign = 10

    # Cache the normalization coefficient once
    norm_coef = maximum(conn_weights)

    # Drop weak edges from the graph before rendering to reduce GPU load.
    # Filtering colors after the fact still uploads all edge geometry to the GPU.
    edge_quantile = get(opts, :edge_quantile, 0.3)
    if edge_quantile > 0.0 && !haskey(opts, :edge_color)
        w = network.weights  # sparse matrix (src → dst)
        nz = SparseArrays.nonzeros(w)
        if !isempty(nz)
            threshold = quantile(nz, edge_quantile)
            rows, cols, vals = SparseArrays.findnz(w)
            mask = vals .>= threshold
            filtered_w = SparseArrays.sparse(
                rows[mask], cols[mask], vals[mask], size(w)...
            )
            network = SimpleWeightedDiGraph(filtered_w)
        end
    end

    # Calculate alpha values for edges based on connectivity strength and weighting (lazy, only if needed)
    edge_col = get(opts, :edge_color) do
        RGBAf.(
            0,
            0,
            0,
            [
                (e.src == e.dst) ? 0.0f0 :
                Float32(conn_weights[e.src] * e.weight / norm_coef)
                for e in edges(network)
            ]
        )
    end

    # Rescale node size to be visible
    node_size = get(opts, :node_size, conn_weights ./ norm_coef .* 10.0)
    node_color = get(opts, :node_color, node_size)

    # Coastlines rendered before graph so nodes/edges appear on top
    show_coastlines && _render_coastlines!(spatial, coastline_resolution)

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

    # Map decorations (places, scale bar, north arrow) rendered last
    if show_coastal_places || show_scale_bar || show_north_arrow
        deco_data = compute_map_decorations(
            dom.loc_data; min_population=min_population, max_km=max_km
        )
        if !isnothing(deco_data)
            limits!(spatial, deco_data.lon_range, deco_data.lat_range)
            _render_map_decorations!(
                spatial,
                dom.loc_data;
                min_population=min_population,
                max_km=max_km,
                show_coastal_places=show_coastal_places,
                show_scale_bar=show_scale_bar,
                show_north_arrow=show_north_arrow,
                decoration_fontsize=decoration_fontsize,
                north_arrow_fontsize=north_arrow_fontsize
            )
        end
    end

    return g
end
