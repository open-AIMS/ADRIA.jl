"""
Shared spatial visualization utilities for all backends.

Provides location ID resolution, cartographic helpers, and decoration data
structures used by both Plotly and Makie spatial viz backends.
"""

using DataFrames: DataFrame, nrow, hasproperty
using Distances: haversine
using ADRIA: Domain, ResultSet, centroids

# =============================================================================
# Location ID Resolution
# =============================================================================

# Candidate location-id columns, consulted when no domain context is available.
const _SITE_ID_CANDIDATES = (:reef_siteid, :site_id, :UNIQUE_ID, :reef_id)

"""
    _loc_id_col(x::Union{Domain,ResultSet}) -> Symbol

The canonical per-location identifier column for the modelled spatial scale.
This is the column ADRIA uses to order `conn`/`centroids`, so it is the only
correct choice for matching spatial features to data and for labelling sites.
Never guess column names when this is available.
"""
_loc_id_col(dom::Domain)::Symbol = Symbol(dom.loc_id_col)
_loc_id_col(rs::ResultSet)::Symbol = Symbol(rs.env_layer_md.loc_id_col)

"""
    _get_site_ids(gdf::DataFrame) -> Vector{String}

Fallback site-id resolver when no Domain/ResultSet provides the canonical
`loc_id_col`. Tries `_SITE_ID_CANDIDATES` in order, then falls back to
row indices. Warns when the choice is ambiguous.
"""
function _get_site_ids(gdf::DataFrame)::Vector{String}
    present = [c for c in _SITE_ID_CANDIDATES if hasproperty(gdf, c)]
    if isempty(present)
        return string.(1:nrow(gdf))
    end
    if length(present) > 1
        @warn "Multiple candidate location-id columns present $(present); using \
               '$(first(present))'. Pass a Domain/ResultSet so the canonical \
               loc_id_col is used instead."
    end
    return string.(gdf[:, first(present)])
end

"""
    _site_ids(gdf, id_col) -> Vector{String}

Resolve site ids using an explicit column when provided (domain's `loc_id_col`),
otherwise fall back to column-guessing via `_get_site_ids`.
"""
function _site_ids(gdf::DataFrame, id_col::Union{Symbol,Nothing})::Vector{String}
    if !isnothing(id_col)
        if hasproperty(gdf, id_col)
            return string.(gdf[:, id_col])
        end
        @warn "Configured location-id column '$(id_col)' not found in \
               GeoDataFrame; falling back to heuristic resolution."
    end
    return _get_site_ids(gdf)
end

# =============================================================================
# Cartographic Helpers
# =============================================================================

"""
    _haversine_km(lon1::Real, lat1::Real, lon2::Real, lat2::Real)::Float64

Great-circle distance between two lon/lat points (degrees) in kilometres.
Uses Distances.jl haversine with WGS84 Earth radius (6371 km).
"""
function _haversine_km(lon1::Real, lat1::Real, lon2::Real, lat2::Real)::Float64
    6371.0 * haversine([deg2rad(lat1), deg2rad(lon1)], [deg2rad(lat2), deg2rad(lon2)])
end

"""
    _nice_length(x::Real)::Int

Round a distance down to a cartographically 'nice' scale-bar length (km).
"""
function _nice_length(x::Real)::Int
    options = [1, 2, 5, 10, 20, 25, 50, 100, 200, 500, 1000]
    idx = findlast(<=(x), options)
    return options[isnothing(idx) ? 1 : idx]
end

# =============================================================================
# Coastal Context Data
# =============================================================================

"""Major Queensland coastal population centres (name, lon, lat)."""
const GBR_COASTAL_PLACES = [
    ("Cooktown", 145.251, -15.468), ("Port Douglas", 145.465, -16.483),
    ("Cairns", 145.770, -16.920), ("Innisfail", 146.030, -17.524),
    ("Cardwell", 146.020, -18.267), ("Ingham", 146.158, -18.650),
    ("Townsville", 146.816, -19.258), ("Ayr", 147.405, -19.574),
    ("Bowen", 148.246, -20.013), ("Airlie Beach", 148.718, -20.267),
    ("Mackay", 149.186, -21.144), ("Yeppoon", 150.742, -23.130),
    ("Rockhampton", 150.510, -23.378), ("Gladstone", 151.256, -23.843),
    ("Bundaberg", 152.349, -24.866), ("Hervey Bay", 152.844, -25.290),
    ("Brisbane", 153.025, -27.470)
]

"""
    CoastalPlace

Named tuple for a coastal place: (name, lon, lat).
"""
const CoastalPlace = NamedTuple{(:name, :lon, :lat),Tuple{String,Float64,Float64}}

"""
    MapDecorationData

Computed data for map decorations: nearby places, scale bar specs, and display extent.
"""
struct MapDecorationData
    places::Vector{CoastalPlace}
    scale_bar_km::Int
    scale_bar_deg::Float64
    scale_bar_x0::Float64
    scale_bar_y::Float64
    lon_range::Vector{Float64}
    lat_range::Vector{Float64}
end

"""
    compute_map_decorations(gdf::DataFrame; max_km=100.0) -> Union{MapDecorationData,Nothing}

Compute decoration data (nearby coastal places, scale bar) for a reef map.

# Arguments
- `gdf::DataFrame`: GeoDataFrame with geometries in WGS84 (EPSG:4326)
- `max_km::Float64`: Maximum distance in km to include coastal places

# Returns
- `MapDecorationData` with computed places, scale bar specs, and display extent
- `Nothing` if centroids cannot be computed

# Assumptions
- Geometries are projected to WGS84 (EPSG:4326, lon/lat coordinates)
- `ADRIA.centroids(gdf)` returns (lon, lat) tuples in degrees
- Display is within ±75° latitude (near-pole scale bar rendering not supported)
"""
function compute_map_decorations(gdf::DataFrame; max_km::Float64=100.0)
    cents = centroids(gdf)
    isempty(cents) && return nothing

    lons = Float64[c[1] for c in cents]
    lats = Float64[c[2] for c in cents]
    lon_min, lon_max = extrema(lons)
    lat_min, lat_max = extrema(lats)

    # Coastal centres within max_km of any reef
    places = CoastalPlace[]
    for (name, plon, plat) in GBR_COASTAL_PLACES
        d = minimum(_haversine_km(plon, plat, lons[i], lats[i]) for i in eachindex(lons))
        d <= max_km && push!(places, (name=name, lon=plon, lat=plat))
    end

    # Display range = reefs + qualifying places, with small pad
    all_lons = vcat(lons, [p.lon for p in places])
    all_lats = vcat(lats, [p.lat for p in places])
    lo_min, lo_max = extrema(all_lons)
    la_min, la_max = extrema(all_lats)
    lon_pad = max(0.05, 0.08 * (lo_max - lo_min))
    lat_pad = max(0.05, 0.08 * (la_max - la_min))
    lon_range = [lo_min - lon_pad, lo_max + lon_pad]
    lat_range = [la_min - lat_pad, la_max + lat_pad]

    # Scale bar (~1/4 of reef extent width, rounded to nice length)
    lat_mid = (lat_min + lat_max) / 2

    # Guard against polar regions (cosd near poles → 0, causing division issues)
    if abs(lat_mid) > 75.0
        @warn "Map decorations requested at high latitude (lat=$lat_mid); \
               skipping scale bar computation for near-polar regions."
        return nothing
    end

    width_km = _haversine_km(lon_min, lat_mid, lon_max, lat_mid)
    bar_km = _nice_length(max(width_km / 4, 1.0))

    cos_lat = cosd(lat_mid)
    if cos_lat < 1e-10  # Shouldn't happen with guard above, but be defensive
        @warn "Scale bar computation failed: degenerate latitude ($lat_mid)"
        return nothing
    end
    bar_deg = bar_km / (111.32 * cos_lat)

    bx0 = lon_range[1] + 0.04 * (lon_range[2] - lon_range[1])
    by = lat_range[1] + 0.06 * (lat_range[2] - lat_range[1])

    return MapDecorationData(
        places, bar_km, bar_deg, bx0, by, lon_range, lat_range
    )
end
