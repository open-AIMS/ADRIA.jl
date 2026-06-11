using Pkg; Pkg.activate(".")
using CSV, DataFrames, GeoDataFrames
using NetCDF

# Load sites to get exactly 2914 UNIQUE_IDs in order
lizard_sites = GeoDataFrames.read("sandbox/data/Lizard_Historical_v0.1/spatial/lizard_cluster.gpkg")
unique_ids = lizard_sites.UNIQUE_ID

# Load the filtered CSV
dhw_df = CSV.read("sandbox/data/dhw_historical.csv", DataFrame)

# Dimensions
years = 1985:2024
n_years = length(years)
n_sites = length(unique_ids)

# Create mapping from UNIQUE_ID to site index
id_to_idx = Dict(unique_ids[i] => i for i in 1:n_sites)

# Create DHW array
dhw_matrix = zeros(Float32, n_years, n_sites, 2)

# Fill array
for row in eachrow(dhw_df)
    id = row.RME_UNIQUE_ID
    if haskey(id_to_idx, id)
        site_idx = id_to_idx[id]
        year_idx = row.timestep - 1984 # 1985 -> 1
        
        dhw_matrix[year_idx, site_idx, 1] = row.dhw
        dhw_matrix[year_idx, site_idx, 2] = row.dhw
    end
end

dhw_file = "sandbox/data/Lizard_Historical_v0.1/DHWs/dhw_RCPhistorical.nc"
isfile(dhw_file) && rm(dhw_file)

nccreate(dhw_file, "dhw", 
    "time", collect(1:n_years), 
    "location", collect(1:n_sites), 
    "scenario", [1, 2])
ncwrite(dhw_matrix, dhw_file, "dhw")

println("Created dhw_historical.nc with size ", size(dhw_matrix))
