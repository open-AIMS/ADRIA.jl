using Pkg; Pkg.activate(".")
using GeoDataFrames
using ArchGDAL
using Proj

# Load existing domain spatial data
gpkg_path = "sandbox/data/Lizard_Island_Cluster_v0.1/spatial/lizard_cluster.gpkg"
lizard_sites = GeoDataFrames.read(gpkg_path)

# Load raster
raster_path = "sandbox/data/COTS_prob_0.02_cpue_year2025_clean.tif"
cots_densities = zeros(Float64, nrow(lizard_sites))

ArchGDAL.read(raster_path) do dataset
    band = ArchGDAL.getband(dataset, 1)
    geotransform = ArchGDAL.getgeotransform(dataset)
    # Proj transformer from WGS84 (EPSG:4326) to GDA2020/MGA zone 55 (EPSG:7855)
    trans = Proj.Transformation("EPSG:4326", "EPSG:7855", always_xy=true)

    # Simple nearest-neighbor extraction based on centroids
    for i in 1:nrow(lizard_sites)
        lon, lat = lizard_sites.x_coord[i], lizard_sites.y_coord[i]
        
        # Transform coords to match raster
        x, y = trans(lon, lat)
        
        # Calculate pixel indices
        pixel = round(Int, (x - geotransform[1]) / geotransform[2])
        line = round(Int, (y - geotransform[4]) / geotransform[6])
        
        # Check bounds
        if 0 <= pixel < ArchGDAL.width(dataset) && 0 <= line < ArchGDAL.height(dataset)
            val = ArchGDAL.read(dataset, 1, pixel, line, 1, 1)[1]
            cots_densities[i] = Float64(val)
        end
    end
end

lizard_sites.cots_density = cots_densities

# Delete existing layer and write new one
GeoDataFrames.write(gpkg_path, lizard_sites; layer_name="lizard_cluster")
println("Successfully added COTS densities to spatial data.")
