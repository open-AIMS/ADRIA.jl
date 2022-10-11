function create_map!(f, geodata, data, highlight, centroids)
    lon = [c[1] for c in centroids]
    lat = [c[2] for c in centroids]

    map_buffer = 0.005
    spatial = GeoAxis(
        f[1,2];
        lonlims=(minimum(lon) - map_buffer, maximum(lon) + map_buffer),
        latlims=(minimum(lat) - map_buffer, maximum(lat) + map_buffer),
        xlabel="Long",
        ylabel="Lat"
    )
    # datalims!(spatial)  # auto-adjust limits (doesn't work if there are Infs...)
    m_b = @lift(maximum($data[:]))
    
    poly!(spatial, geodata, color=data, colormap=:plasma, colorrange=(0.0, m_b[]), strokecolor=highlight, strokewidth=2.0)
    Colorbar(f[1,1]; colorrange=(0.0, m_b[]),
             colormap=:plasma, label="Relative Cover", height=Relative(0.65))

    return spatial
end
