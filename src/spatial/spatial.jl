"""Functions to interact with spatial datasets."""

"""
    _get_geom_col(gdf::DataFrame)::Union{Symbol, Bool}

Retrieve first column found to hold plottable geometries.

# Returns
Symbol, indicating column name or `false` if no geometries found.
"""
function _get_geom_col(gdf::DataFrame)::Union{Symbol,Bool}
    col = findall(typeof.(eachcol(gdf)) .<: Vector{<:AG.IGeometry})
    if !isempty(col)
        return Symbol(names(gdf)[col[1]])
    end

    return false
end

function get_geometry(df::DataFrame)
    if columnindex(df, :geometry) > 0
        return df.geometry
    elseif columnindex(df, :geom) > 0
        return df.geom
    end

    return error("No geometry data found")
end

"""
    centroids(df::DataFrame)

Extract and return long/lat from a GeoDataFrame.

# Arguments
- `df` : GeoDataFrame

# Returns
Array of tuples (x, y), where x and y relate to long and lat respectively.
"""
function centroids(df::DataFrame)::Vector{Tuple{Float64,Float64}}
    site_centroids::Vector = AG.centroid.(get_geometry(df))
    return collect(zip(AG.getx.(site_centroids, 0), AG.gety.(site_centroids, 0)))
end
function centroids(ds::Union{Domain,ResultSet})::Vector{Tuple{Float64,Float64}}
    return centroids(ds.site_data)
end
