import GeoMakie.GeoJSON: AbstractFeatureCollection, features, bbox

# Temporary monkey-patch to support retrieval of multiple features
Base.getindex(fc::AbstractFeatureCollection, i::UnitRange) = features(fc)[i]
Base.getindex(fc::AbstractFeatureCollection, i::Vector) = features(fc)[i]

"""
    create_map!(f::Union{GridLayout,GridPosition}, geodata::GeoMakie.GeoJSON.FeatureCollection,
        data::Observable, highlight::Union{Vector,Tuple,Nothing},
        centroids::Vector, show_colorbar::Bool=true, colorbar_label::String="",
        legend_params::Union{Tuple,Nothing}=nothing, axis_opts::Dict=Dict())

Create a spatial choropleth figure.

# Arguments
- `f` : Makie figure to create plot in
- `geodata` : FeatureCollection, Geospatial data to display
- `data` : Values to use for choropleth
- `highlight` : Stroke colors for each location
- `centroids` : Vector{Tuple}, of lon and lats
- `show_colorbar` : Whether to show a colorbar (true) or not (false)
- `colorbar_label` : Label to use for color bar
- `colorbar_limits` : Upper and lower limits displayed on colorbar,
    (default is (0.0, maximum(data)))
- `color_map` : Type of colormap to use,
    See: https://docs.makie.org/stable/documentation/colors/#colormaps
- `legend_params` : Legend parameters
- `axis_opts` : Additional options to pass to adjust Axis attributes
  See: https://docs.makie.org/v0.19/api/index.html#Axis
"""
function create_map!(
    f::Union{GridLayout, GridPosition},
    geodata::GeoMakie.GeoJSON.FeatureCollection,
    data::Observable,
    highlight::Union{Vector, Tuple, Nothing},
    centroids::Vector,
    show_colorbar::Bool = true,
    colorbar_label::String = "",
    colorbar_limits::Tuple{Float64, Float64} = (0.0, maximum(data)),
    color_map::Union{Symbol, Vector{Symbol}, RGBA{Float32}, Vector{RGBA{Float32}}} = :grayC,
    legend_params::Union{Tuple, Nothing} = nothing,
    axis_opts::Dict = Dict(),
)
    axis_opts[:title] = get(axis_opts, :title, "Study Area")
    axis_opts[:xlabel] = get(axis_opts, :xlabel, "Longitude")
    axis_opts[:ylabel] = get(axis_opts, :ylabel, "Latitude")

    spatial = GeoAxis(
        f[1, 1];
        dest = "+proj=latlong +datum=WGS84",
        axis_opts...,
    )
    # lon = first.(centroids)
    # lat = last.(centroids)
    # map_buffer = 0.025
    # xlims!(spatial, minimum(lon) - map_buffer, maximum(lon) + map_buffer)
    # ylims!(spatial, minimum(lat) - map_buffer, maximum(lat) + map_buffer)

    spatial.xticklabelsize = 14
    spatial.yticklabelsize = 14

    spatial.yticklabelpad = 50
    spatial.ytickalign = 10


    poly!(
        spatial,
        geodata;
        color = data,
        colormap = color_map,
        strokecolor = (:black, 0.05),
        strokewidth = 1.0,
    )

    if show_colorbar
        Colorbar(
            f[1, 2];
            colormap = color_map,
            label = colorbar_label,
            height = Relative(0.65),
            limits = colorbar_limits,
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
                color = "transparent",
                strokecolor = highlight,
                strokewidth = 0.5,
                linestyle = :solid,
                overdraw = true,
            )
        else
            hl_groups = unique(highlight)

            for color in hl_groups
                m = findall(highlight .== [color])
                subset_feat = FC(; features = geodata[m])

                poly!(
                    spatial,
                    subset_feat;
                    color = "transparent",
                    strokecolor = color,
                    strokewidth = 0.5,
                    linestyle = :solid,
                    overdraw = true,
                )
            end
        end

        if !isnothing(legend_params)
            # Plot Legend only if highlight colors are present
            Legend(f[1, 3], legend_params...; framevisible = false)
        end
    end

    return f
end

"""
    ADRIA.viz.map(rs::Union{Domain,ResultSet}; opts=Dict(by_RCP => false), fig_opts=Dict(), 
        axis_opts=Dict(), series_opts=Dict())
    ADRIA.viz.map(rs::ResultSet, y::NamedDimsArray; opts=Dict(by_RCP => false), fig_opts=Dict(), 
        axis_opts=Dict(), series_opts=Dict())
    ADRIA.viz.map(rs::ResultSet, S::NamedDimsArray, scores::Vector{Float64};
        criteria::Vector{Symbol} = S.criteria, opts::Dict = Dict(), axis_opts::Dict = Dict(),
        fig_opts::Dict = Dict())
    ADRIA.viz.map!(f::Union{GridLayout,GridPosition}, rs::ADRIA.ResultSet, y::NamedDimsArray; 
        opts=Dict(by_RCP => false), axis_opts=Dict(), series_opts=Dict())
    ADRIA.viz.map!(g::Union{GridLayout,GridPosition},rs::ResultSet, S::NamedDimsArray, 
        scores::Vector{Float64}; criteria::Vector{Symbol} = S.criteria, opts::Dict = Dict(), 
        axis_opts::Dict = Dict(), fig_opts::Dict = Dict())

Plot spatial choropleth of outcomes.

# Arguments
- `rs` : ResultSet
- `y` : results of scenario metric
- `S` : A normalised decision matrix calculated using decison.decision_matrices
- `scores` : Aggregated criteria scores.
- `criteria` : Names of criteria to be plotted, if not specified all criteria in 
    S will be plotted.
- `opts` : Aviz options
    - `colorbar_label`, label for colorbar. Defaults to "Relative Cover"
    -`colorbar_limits`, min and max values to be shown on the colorbar. 
        Defaults to (0.0,maximum(y)).
    - `color_map`, preferred colormap for plotting heatmaps
- `axis_opts` : Additional options to pass to adjust Axis attributes
  See: https://docs.makie.org/v0.19/api/index.html#Axis
- `series_opts` : Additional options to pass to adjust Series attributes
  See: https://docs.makie.org/v0.19/api/index.html#series!

# Returns
GridPosition
"""
function ADRIA.viz.map(
    rs::Union{Domain, ResultSet},
    y::NamedDimsArray;
    opts::Dict=Dict(),
    fig_opts::Dict=Dict(),
    axis_opts::Dict=Dict(),
)
    f = Figure(; fig_opts...)
    g = f[1, 1] = GridLayout()

    ADRIA.viz.map!(g, rs, collect(y); opts, axis_opts)

    return f
end
function ADRIA.viz.map(
    rs::Union{Domain, ResultSet};
    opts::Dict{Symbol, <:Any}=Dict{Symbol, Any}(),
    fig_opts::Dict{Symbol, <:Any}=Dict{Symbol, Any}(),
    axis_opts::Dict{Symbol, <:Any}=Dict{Symbol, Any}(),
)
    f = Figure(; fig_opts...)
    g = f[1, 1] = GridLayout()

    opts[:colorbar_label] = get(opts, :colorbar_label, "Coral Real Estate [%]")

    opts[:show_management_zones] = get(opts, :show_management_zones, false)
    if opts[:show_management_zones]
        highlight = Symbol.(lowercase.(rs.site_data.zone_type))
        opts[:highlight] = highlight
    end

    ADRIA.viz.map!(g, rs, rs.site_data.k; opts, axis_opts)

    return f
end
function ADRIA.viz.map!(
    g::Union{GridLayout, GridPosition},
    rs::Union{Domain, ResultSet},
    y::AbstractVector{<:Real};
    opts::Dict=Dict(),
    axis_opts::Dict=Dict(),
)
    geodata = get_geojson_copy(rs)
    data = Observable(collect(y))

    highlight = get(opts, :highlight, nothing)
    c_label = get(opts, :colorbar_label, "")
    legend_params = get(opts, :legend_params, nothing)
    show_colorbar = get(opts, :show_colorbar, true)
    color_map = get(opts, :color_map, :grayC)
    colorbar_limits = get(opts, :colorbar_limits, (0.0, maximum(y)))
    return create_map!(
        g,
        geodata,
        data,
        highlight,
        ADRIA.centroids(rs),
        show_colorbar,
        c_label,
        colorbar_limits,
        color_map,
        legend_params,
        axis_opts,
    )
end
function ADRIA.viz.map(
    rs::ResultSet,
    S::NamedDimsArray,
    scores::Vector{Float64};
    criteria::Vector{Symbol} = S.criteria,
    opts::Dict = Dict(),
    axis_opts::Dict = Dict(),
    fig_opts::Dict = Dict(),
)
    f = Figure(; fig_opts...)
    g = f[1, 1] = GridLayout()
    ADRIA.viz.map!(
        g, rs, S, scores; criteria = criteria, opts = opts, axis_opts = axis_opts
    )
    return f
end
function ADRIA.viz.map!(
    g::Union{GridLayout, GridPosition},
    rs::ResultSet,
    S::NamedDimsArray,
    scores::Vector{Float64};
    criteria::Vector{Symbol} = S.criteria,
    opts::Dict = Dict(),
    axis_opts::Dict = Dict(),
)
    if length(rs.site_data.site_id) != size(S, 1)
        error("Only unfiltered decision matrices can be plotted.")
    end

    opts[:color_map] = get(opts, :color_map, :viridis)
    opts[:colorbar_limits] = get(opts, :colorbar_limits, (0.0, 1.0))

    m_spec = model_spec(rs)
    criteria_names::Vector{String} = m_spec[
        dropdims(
            any(
                reshape(criteria, 1, length(criteria)) .== m_spec[:, "fieldname"]; dims = 2
            );
            dims = 2,
        ), "name"]
    n_criteria::Int64 = length(criteria)
    n_rows, n_cols = _calc_gridsize(n_criteria + 1)
    step::Int64 = 1

    for row in 1:n_rows, col in 1:n_cols
        if step > length(criteria_names)
            ADRIA.viz.map!(
                g[row, col],
                rs,
                vec(scores);
                opts = opts,
                axis_opts = Dict(:title => "Aggregate criteria score"; axis_opts...),
            )
            break
        end
        axis_opts_temp = Dict(:title => criteria_names[step]; axis_opts...)
        ADRIA.viz.map!(
            g[row, col],
            rs,
            vec(S(criteria[step]));
            opts = opts,
            axis_opts = axis_opts_temp
        )

        step += 1
    end

    # Clear empty figures
    return trim!(g)
end

"""
    make_geojson_copy(ds::Union{ResultSet,Domain})::String

Make a temporary copy of GeoPackage as GeoJSON.

# Arguments
`ds` : Domain or ResultSet containing spatial data

# Returns
Path to temporary copy of GeoJSON file.
"""
function make_geojson_copy(ds::Union{ResultSet, Domain})::String
    tmpdir = ADRIA.viz.tmpdir
    local geo_fn = joinpath(tmpdir, "Aviz_$(ds.name).geojson")
    if !isfile(geo_fn)
        try
            GDF.write(geo_fn, ds.site_data; driver="geojson")
        catch
            GDF.write(geo_fn, ds.site_data; geom_columns=(:geom,), driver="geojson")
        end
    end

    return geo_fn
end

"""
    get_geojson(ds::Union{ResultSet,Domain})::FC

Retrieves a temporary copy of spatial data associated with the given Domain or ResultSet as
a FeatureCollection.

# Arguments
- `ds` : The dataset with which the spatial data is associated with

# Returns
FeatureCollection of polygons
"""
function get_geojson_copy(ds::Union{ResultSet, Domain})::FC
    fn = make_geojson_copy(ds)

    # Only return the set Features if filepaths match
    if isdefined(ADRIA.viz, :tmp_geojson)
        if ADRIA.viz.tmpdir == dirname(fn)
            return ADRIA.viz.tmp_geojson
        end
    end

    ADRIA.viz.tmp_geojson = GeoMakie.GeoJSON.read(read(fn))
    return ADRIA.viz.tmp_geojson
end
