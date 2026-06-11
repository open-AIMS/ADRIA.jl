using Pkg; Pkg.activate(".")
using NetCDF

# Load the single scenario matrix
dhw_matrix = ncread("sandbox/data/dhw_historical.nc", "dhw")
# shape is (31, 2914, 1)

# duplicate it
dhw_matrix_2 = cat(dhw_matrix, dhw_matrix, dims=3)
# shape is (31, 2914, 2)

dhw_file = "sandbox/data/Lizard_Historical_v0.1/DHWs/dhw_RCPhistorical.nc"
rm(dhw_file, force=true)

nccreate(dhw_file, "dhw", 
    "time", collect(1:31), 
    "location", collect(1:2914), 
    "scenario", [1, 2])
ncwrite(dhw_matrix_2, dhw_file, "dhw")

println("Saved dhw_RCPhistorical.nc with shape ", size(dhw_matrix_2))
