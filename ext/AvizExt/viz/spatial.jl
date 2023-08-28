import GeoMakie.GeoJSON: AbstractFeatureCollection, features, bbox

# Temporary monkey-patch to support retrieval of multiple features
Base.getindex(fc::AbstractFeatureCollection, i::UnitRange) = features(fc)[i]
Base.getindex(fc::AbstractFeatureCollection, i::Vector) = features(fc)[i]


"""
    create_map!(f::Union{GridLayout,GridPosition}, geodata::GeoMakie.GeoJSON.FeatureCollection,
        data::Observable, highlight::Union{Vector,Tuple,Nothing},
        centroids::Vector, colorbar_label::String="",
        legend_params::Union{Tuple,Nothing}=nothing, axis_opts::Dict=Dict())

Create a spatial choropleth figure.

# Arguments
- `f` : GLMakie figure to create plot in
- `geodata` : FeatureCollection, Geospatial data to display
- `data` : Values to use for choropleth
- `highlight` : Stroke colors for each location
- `centroids` : Vector{Tuple}, of lon and lats
- `colorbar_label` : Label to use for color bar
- `legend_params` : Legend parameters
- `axis_opts` : Additional options to pass to adjust Axis attributes
  See: https://docs.makie.org/v0.19/api/index.html#Axis
"""
function create_map!(f::Union{GridLayout,GridPosition}, geodata::GeoMakie.GeoJSON.FeatureCollection,
    data::Observable, highlight::Union{Vector,Tuple,Nothing},
    centroids::Vector, colorbar_label::String="",
    legend_params::Union{Tuple,Nothing}=nothing, axis_opts::Dict=Dict())
    lon = first.(centroids)
    lat = last.(centroids)

    axis_opts[:title] = get(axis_opts, :title, "Study Area")
    axis_opts[:xlabel] = get(axis_opts, :xlabel, "Longitude")
    axis_opts[:ylabel] = get(axis_opts, :ylabel, "Latitude")

    map_buffer = 0.025
    spatial = GeoAxis(
        f[1, 1];
        lonlims=(minimum(lon) - map_buffer, maximum(lon) + map_buffer),
        latlims=(minimum(lat) - map_buffer, maximum(lat) + map_buffer),
        dest="+proj=latlong +datum=WGS84",
        axis_opts...
    )
    spatial.xticklabelsize = 12
    spatial.yticklabelsize = 12

    # spatial.xticklabelsvisible = false
    # spatial.yticklabelsvisible = false

    spatial.yticklabelpad = 50
    spatial.ytickalign = 10
    max_val = @lift(maximum($data))

    # Plot geodata polygons using data as internal color
    color_range = (0.0, max_val[])
    color_map = :grayC
    poly!(
        spatial,
        geodata,
        color=data,
        colormap=color_map,
        colorrange=color_range,
        strokecolor=(:black, 0.05),
        strokewidth=1.0
    )
    Colorbar(f[1, 2]; colorrange=color_range, colormap=color_map, label=colorbar_label,
        height=Relative(0.65))

    # Overlay locations to be highlighted
    # `poly!()` cannot handle multiple strokecolors being specified at the moment
    # so we instead overlay each cluster.
    if !isnothing(highlight)
        if highlight isa Tuple
            poly!(
                spatial,
                geodata,
                color="transparent",
                strokecolor=highlight,
                strokewidth=0.5,
                linestyle=:solid,
                overdraw=true,
            )
        else
            hl_groups = unique(highlight)

            for color in hl_groups
                m = findall(highlight .== [color])
                subset_feat = FC(; features=geodata[m])

                poly!(
                    spatial,
                    subset_feat,
                    color="transparent",
                    strokecolor=color,
                    strokewidth=0.5,
                    linestyle=:solid,
                    overdraw=true,
                )
            end
        end

        if !isnothing(legend_params)
            # Plot Legend only if highlight colors are present
            Legend(f[1, 3], legend_params..., framevisible=false)
        end
    end
    # datalims!(spatial)  # auto-adjust limits (doesn't work if there are Infs...)

    return f
end


"""
	ADRIA.viz.map(rs::Union{Domain,ResultSet}; opts=Dict(by_RCP => false), fig_opts=Dict(), axis_opts=Dict(), series_opts=Dict())
	ADRIA.viz.map(rs::ResultSet, y::NamedDimsArray; opts=Dict(by_RCP => false), fig_opts=Dict(), axis_opts=Dict(), series_opts=Dict())
	ADRIA.viz.map!(f::Union{GridLayout,GridPosition}, rs::ADRIA.ResultSet, y::NamedDimsArray; opts=Dict(by_RCP => false), axis_opts=Dict(), series_opts=Dict())

Plot spatial choropleth of outcomes.

# Arguments
- `rs` : ResultSet
- `y` : results of scenario metric
- `opts` : Aviz options
	- `colorbar_label`, label for colorbar. Defaults to "Relative Cover".
- `axis_opts` : Additional options to pass to adjust Axis attributes
  See: https://docs.makie.org/v0.19/api/index.html#Axis
- `series_opts` : Additional options to pass to adjust Series attributes
  See: https://docs.makie.org/v0.19/api/index.html#series!

# Returns
GridPosition
"""
function ADRIA.viz.map(rs::Union{Domain,ResultSet}, y::NamedDimsArray; opts::Dict=Dict(),
    fig_opts::Dict=Dict(), axis_opts::Dict=Dict())
    f = Figure(; fig_opts...)
    g = f[1, 1] = GridLayout()

    ADRIA.viz.map!(g, rs, collect(y); opts, axis_opts)

    return f
end
function ADRIA.viz.map(rs::Union{Domain,ResultSet}; opts::Dict=Dict(), fig_opts::Dict=Dict(),
    axis_opts::Dict=Dict())
    f = Figure(; fig_opts...)
    g = f[1, 1] = GridLayout()

    opts[:colorbar_label] = get(opts, :colorbar_label, "Coral Real Estate [%]")
    ADRIA.viz.map!(g, rs, rs.site_data.k; opts, axis_opts)

    return f
end
function ADRIA.viz.map!(g::Union{GridLayout,GridPosition}, rs::Union{Domain,ResultSet},
    y::Vector; opts::Dict=Dict(), axis_opts::Dict=Dict())

    geodata = get_geojson_copy(rs)
    data = Observable(y)

    highlight = get(opts, :highlight, nothing)
    c_label = get(opts, :colorbar_label, "")
    legend_params = get(opts, :legend_params, nothing)

    return create_map!(g, geodata, data, highlight, ADRIA.centroids(rs), c_label, legend_params, axis_opts)
end


"""
	make_geojson_copy(ds::Union{ResultSet,Domain})::String

Make a temporary copy of GeoPackage as GeoJSON.

# Arguments
`ds` : Domain or ResultSet containing spatial data

# Returns
Path to temporary copy of GeoJSON file.
"""
function make_geojson_copy(ds::Union{ResultSet,Domain})::String
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
function get_geojson_copy(ds::Union{ResultSet,Domain})::FC
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
