using Pkg; Pkg.activate(".")
using GeoDataFrames

df = GeoDataFrames.read("sandbox/data/Lizard_Island_Cluster_v0.1/spatial/lizard_cluster.gpkg")
densities = coalesce.(df.cots_density, 0.0)
densities = [isnan(x) ? 0.0 : x for x in densities]
println("Sum > 0: ", sum(densities .> 0.0))
println("Sum > 0.1: ", sum(densities .> 0.1))
println("Sum > 0.3: ", sum(densities .> 0.3))
println("Sum > 0.4: ", sum(densities .> 0.4))
println("Sum > 0.5: ", sum(densities .> 0.5))
println("Sum > 0.6: ", sum(densities .> 0.6))
println("Sum > 0.7: ", sum(densities .> 0.7))
println("Maximum: ", maximum(densities))
