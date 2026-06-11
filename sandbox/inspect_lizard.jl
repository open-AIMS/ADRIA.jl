using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
using GeoDataFrames

gpkg_path = joinpath(@__DIR__, "data", "lizard_sites_june.gpkg")
gdf = GeoDataFrames.read(gpkg_path, layer="sites_master")

println("DataFrame Info:")
println(first(gdf, 5))
println("\nColumns: ", names(gdf))
println("Size: ", size(gdf))


