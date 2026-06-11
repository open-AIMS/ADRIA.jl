using Pkg; Pkg.activate(".")
using GeoDataFrames
df = GeoDataFrames.read("sandbox/data/Lizard_Island_Cluster_v0.1/spatial/lizard_cluster.gpkg")
println(names(df))
