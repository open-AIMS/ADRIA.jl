"""
    ADRIA.viz.map(rs::Union{Domain,ResultSet}; opts=Dict(by_RCP => false), fig_opts=Dict(), axis_opts=Dict(), series_opts=Dict())
    ADRIA.viz.map(rs::ResultSet, y::NamedDimsArray; opts=Dict(by_RCP => false), fig_opts=Dict(), axis_opts=Dict(), series_opts=Dict())
    ADRIA.viz.map!(f::Union{GridLayout,GridPosition}, rs::ADRIA.ResultSet, y::NamedDimsArray; opts=Dict(by_RCP => false), axis_opts=Dict(), series_opts=Dict())

Plot spatial outcomes.

# Arguments
- `rs` : ResultSet
- `y` : results of scenario metric
- `opts` : Aviz options 
    - `by_RCP`, color by RCP otherwise color by scenario type. Defaults to false.
- `axis_opts` : Additional options to pass to adjust Axis attributes  
  See: https://docs.makie.org/v0.19/api/index.html#Axis
- `series_opts` : Additional options to pass to adjust Series attributes  
  See: https://docs.makie.org/v0.19/api/index.html#series!

# Returns
GridPosition
"""
function ADRIA.viz.map!(g::Union{GridLayout,GridPosition}, rs::Union{Domain,ResultSet}, y::Vector;
    opts::Dict=Dict(:by_RCP => false), axis_opts::Dict=Dict(), series_opts::Dict=Dict())

    geo_fn = make_geojson_copy(rs)
    geodata = GeoMakie.GeoJSON.read(read(geo_fn))
    data = Observable(y)

    create_map!(g, geodata, data, nothing, ADRIA.centroids(rs))
end
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
