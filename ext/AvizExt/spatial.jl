function create_map!(f, geodata, data, highlight, centroids)
    lon = first.(centroids)
    lat = last.(centroids)

    map_buffer = 0.1
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

    spatial.xticklabelsvisible = false
    spatial.yticklabelsvisible = false

    spatial.yticklabelpad = -40
    spatial.ytickalign = 10
    m_b = @lift(maximum($data[:]))

    poly!(spatial, geodata, color=data, colormap=:plasma, colorrange=(0.0, m_b[]), strokecolor=highlight, strokewidth=2.0)
    datalims!(spatial)  # auto-adjust limits (doesn't work if there are Infs...)
    Colorbar(f[1, 2]; colorrange=(0.0, m_b[]),
        colormap=:plasma, label="Relative Cover", height=Relative(0.65))

    return spatial
end
