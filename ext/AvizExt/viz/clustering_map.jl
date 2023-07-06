function ADRIA.viz.spatial_clustering(rs::Union{Domain,ResultSet},
                                      data::AbstractMatrix,
                                      clusters::Vector{Int64}; 
                                      opts::Dict=Dict(), 
                                      fig_opts::Dict=Dict(), 
                                      axis_opts::Dict=Dict(), 
                                      series_opts::Dict=Dict())
    f = Figure(; fig_opts...)

    ADRIA.viz.tmpdir = mktempdir()
    geo_fn = make_geojson_copy(rs)
    geodata = GeoMakie.GeoJSON.read(read(geo_fn))

    # Create GeoAxis
    centroids = rs.site_centroids
    lon = first.(centroids)
    lat = last.(centroids)

    map_buffer = 0.01
    spatial = GeoAxis(
        f[1, 1];
        lonlims=(minimum(lon) - map_buffer, maximum(lon) + map_buffer),
        latlims=(minimum(lat) - map_buffer, maximum(lat) + map_buffer),
        title="Study Area",
        xlabel="Longitude",
        ylabel="Latitude",
        dest="+proj=latlong +datum=WGS84"
    )

    stroke_colors = cluster_colors(clusters)
    
    data = Observable(rs.site_data.k)
    m_b = @lift(maximum($data[:]))
    
    poly!(spatial, geodata, color=stroke_colors)#, strokecolor=:red, strokewidth=0.4)
    return f
end

function cluster_colors(clusters::Vector{Int64})
    n_clusters = length(unique(filter(c -> c != 0, clusters)))
    cat_colors = categorical_colors(:seaborn_bright, n_clusters)
    return [c == 0 ? "black" : cat_colors[c] for c in clusters]
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