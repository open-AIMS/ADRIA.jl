import GeoMakie.GeoJSON: AbstractFeatureCollection, features, bbox

# Temporary monkey-patch to support retrieval of multiple features
Base.getindex(fc::AbstractFeatureCollection, i::UnitRange) = features(fc)[i]
Base.getindex(fc::AbstractFeatureCollection, i::Vector) = features(fc)[i]


"""
    create_map!(f, geodata, data, highlight, centroids, c_label)

Create a spatial choropleth figure.

# Arguments
- `f` : GLMakie figure to create plot in
- `geodata` : FeatureCollection, Geospatial data to display
- `data` : Vector, values to use for choropleth
- `highlight` : Vector, stroke colors for each location
- `centroids` : Vector{Tuple}, of lon and lats
- `c_label` : String, label to use for color bar
"""
function create_map!(f, geodata, data, highlight, centroids, c_label)
    lon = first.(centroids)
    lat = last.(centroids)

    map_buffer = 0.025
    spatial = GeoAxis(
        f[1, 1];
        lonlims=(minimum(lon) - map_buffer, maximum(lon) + map_buffer),
        latlims=(minimum(lat) - map_buffer, maximum(lat) + map_buffer),
        title="Study Area",
        xlabel="Longitude",
        ylabel="Latitude",
        dest="+proj=latlong +datum=WGS84"
    )
    spatial.xticklabelsize = 12
    spatial.yticklabelsize = 12

    # spatial.xticklabelsvisible = false
    # spatial.yticklabelsvisible = false

    spatial.yticklabelpad = 50
    spatial.ytickalign = 10
    m_b = @lift(maximum($data))

    poly!(spatial, geodata, color=data, colormap=:plasma, colorrange=(0.0, m_b[]), strokecolor=(:black, 0.05))

    # Overlay locations to be highlighted
    # `poly!()` cannot handle multiple strokecolors being specified at the moment
    # so we instead overlay each cluster.
    if !isnothing(highlight)
        hl_groups = unique(highlight)
        for clr in hl_groups
            m = findall(highlight .== [clr])
            subset_feat = FC(; features=geodata[m])
            poly!(spatial, subset_feat, color=data[][m], colormap=:plasma,
                colorrange=(0.0, m_b[]), strokecolor=clr, strokewidth=4.0,
                linestyle=:dashdot, overdraw=true)
        end
    end
    # datalims!(spatial)  # auto-adjust limits (doesn't work if there are Infs...)

    Colorbar(f[1, 2]; colorrange=(0.0, m_b[]),
        colormap=:plasma, label=c_label, height=Relative(0.65))

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
function ADRIA.viz.map(rs::Union{Domain,ResultSet}, y::NamedDimsArray; opts::Dict=Dict(), fig_opts::Dict=Dict(), axis_opts::Dict=Dict(), series_opts::Dict=Dict())
    f = Figure(; fig_opts...)
    g = f[1, 1] = GridLayout()

    ADRIA.viz.map!(g, rs, collect(y); opts, axis_opts, series_opts)

    return f
end
function ADRIA.viz.map(rs::Union{Domain,ResultSet}; opts::Dict=Dict(), fig_opts::Dict=Dict(), axis_opts::Dict=Dict(), series_opts::Dict=Dict())
    f = Figure(; fig_opts...)

    ADRIA.viz.map!(f, rs, rs.site_data.k; opts, axis_opts, series_opts)

    return f
end
function ADRIA.viz.map!(g::Union{GridLayout,GridPosition}, rs::Union{Domain,ResultSet}, y::Vector;
    opts::Dict=Dict(), axis_opts::Dict=Dict(), series_opts::Dict=Dict())

    geo_fn = make_geojson_copy(rs)
    geodata = GeoMakie.GeoJSON.read(read(geo_fn))
    data = Observable(y)

    c_label = get(opts, :colorbar_label, "Relative Cover")
    highlight = get(opts, :highlight, nothing)
    return create_map!(g, geodata, data, highlight, ADRIA.centroids(rs), c_label)
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
