using Pkg; Pkg.activate(".")
using YAXArrays, NetCDF
println("Loading projection DHW:")
try
    dhw = YAXArrays.Cube("sandbox/data/Lizard_Island_Cluster_v0.1/DHWs/dhw_RCP45.nc")
    println(dhw)
catch e
    println("Failed: ", e)
end

println("\nLoading historical DHW:")
try
    dhw2 = YAXArrays.Cube("sandbox/data/Lizard_Historical_v0.1/DHWs/dhw_RCPhistorical.nc")
    println(dhw2)
catch e
    println("Failed: ", e)
end
