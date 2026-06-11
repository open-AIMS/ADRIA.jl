using Pkg; Pkg.activate(".")
using GeoDataFrames
df = GeoDataFrames.read("sandbox/data/Lizard_Island_Cluster_v0.1/spatial/lizard_cluster.gpkg")

println("site_id[1]: ", df.site_id[1], " type: ", typeof(df.site_id[1]))
println("UNIQUE_ID[1]: ", df.UNIQUE_ID[1], " type: ", typeof(df.UNIQUE_ID[1]))
