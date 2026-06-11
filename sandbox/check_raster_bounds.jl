using Pkg; Pkg.activate(".")
using ArchGDAL
using GeoDataFrames

ArchGDAL.read("sandbox/data/COTS_prob_0.02_cpue_year2025_clean.tif") do ds
    println("Raster Geotransform: ", ArchGDAL.getgeotransform(ds))
    println("Raster bounds: ", ArchGDAL.width(ds), "x", ArchGDAL.height(ds))
    band = ArchGDAL.getband(ds, 1)
    println("Raster max value: ", maximum(ArchGDAL.read(band)))
    println("Raster min value: ", minimum(ArchGDAL.read(band)))
end

df = GeoDataFrames.read("sandbox/data/Lizard_Island_Cluster_v0.1/spatial/lizard_cluster.gpkg")
println("Points min x/y: ", minimum(df.x_coord), " ", minimum(df.y_coord))
println("Points max x/y: ", maximum(df.x_coord), " ", maximum(df.y_coord))
