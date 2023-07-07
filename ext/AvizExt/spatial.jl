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

    return spatial
end
