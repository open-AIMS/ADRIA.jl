using Pkg; Pkg.activate(".")
using NetCDF
println("== Projection ==")
ncinfo("sandbox/data/Lizard_Island_Cluster_v0.1/DHWs/dhw_RCP45.nc")
println("\n== Historical ==")
ncinfo("sandbox/data/Lizard_Historical_v0.1/DHWs/dhw_RCPhistorical.nc")
