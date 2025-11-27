"""Functions to interact with spatial datasets."""

"""
    _get_geom_col(gdf::DataFrame)::Union{Symbol, Bool}

Retrieve first column found to hold plottable geometries.

# Returns
Symbol, indicating column name or `false` if no geometries found.
"""
function _get_geom_col(gdf::DataFrame)::Union{Symbol,Bool}
    possible_cols = GeoInterface.geometrycolumns(gdf)
    if length(possible_cols) > 0
        return first(possible_cols)
    end

    return false
end

function get_geometry(df::DataFrame)
    col = _get_geom_col(df)
    if col == false
        throw(ArgumentError("No geometry data found"))
    end

    return df[:, col]
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
    loc_centroids::Vector = AG.centroid.(get_geometry(df))
    return collect(zip(AG.getx.(loc_centroids, 0), AG.gety.(loc_centroids, 0)))
end
function centroids(ds::Union{Domain,ResultSet})::Vector{Tuple{Float64,Float64}}
    return centroids(ds.loc_data)
end
