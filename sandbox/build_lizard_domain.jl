using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
using ADRIA
using GeoDataFrames
using DataFrames
using YAXArrays
using NetCDF
using Statistics
import ArchGDAL as AG
using JSON

# 1. Paths
lizard_gpkg_path = joinpath(@__DIR__, "data", "lizard_sites_june.gpkg")
rme_domain_path = joinpath(@__DIR__, "data", "rme_ml_2025_06_05")
out_dir = joinpath(@__DIR__, "data", "Lizard_Island_Cluster_v0.1")

# Ensure output directories exist
mkpath(out_dir)
mkpath(joinpath(out_dir, "spatial"))
mkpath(joinpath(out_dir, "connectivity"))
mkpath(joinpath(out_dir, "DHWs"))
mkpath(joinpath(out_dir, "waves"))
mkpath(joinpath(out_dir, "cyclones"))

# 2. Load Lizard Island Sites
lizard_sites = GeoDataFrames.read(lizard_gpkg_path, layer="sites_master")

# Calculate centroids if needed
function get_centroid(geom)
    try
        centroid = AG.centroid(geom)
        return (AG.getx(centroid, 0), AG.gety(centroid, 0))
    catch
        return (0.0, 0.0)
    end
end

centroids = get_centroid.(lizard_sites.geometry)
lizard_sites.x_coord = [c[1] for c in centroids]
lizard_sites.y_coord = [c[2] for c in centroids]

# Add standard columns
lizard_sites.k .= 0.5
lizard_sites.depth_med .= 5.0
lizard_sites.recs_per_m2 .= 10.0 # Default value for recruitment

# Save spatial data
# Convert to simple geometries and set CRS
lizard_gpkg_out = joinpath(out_dir, "spatial", "lizard_cluster.gpkg")
# Just write it out
GeoDataFrames.write(lizard_gpkg_out, lizard_sites; layer_name="lizard_cluster")

# 3. Load RME Domain to extract matching data
dom = ADRIA.load_domain(ADRIA.RMEDomain, rme_domain_path, "45")

# 4. Map Lizard sites to RME unique IDs
lizard_rme_ids = lizard_sites.UNIQUE_ID
unique_rme_ids = unique(lizard_rme_ids)

# Find corresponding indices in RME domain
# Note: dom.loc_data.UNIQUE_ID might be a string
rme_idx_map = Dict{String, Int}()
for (i, rme_id) in enumerate(dom.loc_data.UNIQUE_ID)
    rme_idx_map[rme_id] = i
end

# Identify overlapping reefs
matched_rme_ids = filter(id -> haskey(rme_idx_map, id), unique_rme_ids)
println("Matched $(length(matched_rme_ids)) out of $(length(unique_rme_ids)) RME IDs.")

# Build mapping from site index to RME index
site_to_rme_idx = [haskey(rme_idx_map, id) ? rme_idx_map[id] : -1 for id in lizard_rme_ids]

# Count sites per reef for connectivity scaling
sites_per_reef = Dict{String, Int}()
for id in lizard_rme_ids
    sites_per_reef[id] = get(sites_per_reef, id, 0) + 1
end

# 5. Extract and build Connectivity
# We have 2914 sites.
n_sites = nrow(lizard_sites)
conn = zeros(Float64, n_sites, n_sites)

for sink in 1:n_sites
    sink_rme_id = lizard_rme_ids[sink]
    sink_rme_idx = site_to_rme_idx[sink]
    
    if sink_rme_idx == -1
        continue
    end
    
    n_sink_sites = sites_per_reef[sink_rme_id]
    
    for source in 1:n_sites
        source_rme_id = lizard_rme_ids[source]
        source_rme_idx = site_to_rme_idx[source]
        
        if source_rme_idx == -1
            continue
        end
        
        # Base connectivity from RME
        base_conn = dom.conn[source=source_rme_idx, sink=sink_rme_idx][1]
        
        # Simple stopgap: Distribute uniformly among sink sites
        conn[source, sink] = base_conn / n_sink_sites
    end
end

using CSV
using DataFrames
conn_df = DataFrame(conn, :auto)
rename!(conn_df, string.(lizard_sites.site_id); makeunique=true)
CSV.write(joinpath(out_dir, "connectivity", "Lizard_Connectivity.csv"), conn_df)

# 6. Extract DHW
# DHW shape in RME: (timesteps, locations, scenarios)
dhw_r1 = dom.dhw_scens
dhw_sz = size(dhw_r1)
tf = dhw_sz[1]
rme_nlocs = dhw_sz[2]
nscens = length(dhw_sz) > 2 ? dhw_sz[3] : 1

dhw_lizard = zeros(Float32, tf, n_sites, nscens)
for s in 1:n_sites
    rme_idx = site_to_rme_idx[s]
    if rme_idx != -1
        if length(dhw_sz) > 2
            dhw_lizard[:, s, :] .= dhw_r1[:, rme_idx, :]
        else
            dhw_lizard[:, s, 1] .= dhw_r1[:, rme_idx]
        end
    end
end

function write_nc(filename, data, varname, sites, scens, timesteps)
    isfile(filename) && rm(filename)
    if ndims(data) == 3
        nccreate(filename, varname, 
            "time", collect(timesteps), 
            "location", collect(1:length(sites)), 
            "scenario", collect(scens))
        ncwrite(data, filename, varname)
    else
        nccreate(filename, varname, 
            "time", collect(timesteps), 
            "location", collect(1:length(sites)))
        ncwrite(data, filename, varname)
    end
end

dhw_file = joinpath(out_dir, "DHWs", "dhw_RCP45.nc")
write_nc(dhw_file, dhw_lizard, "dhw", lizard_sites.site_id, 1:nscens, 1:tf)

# 7. Extract Initial Cover
# Initial cover shape in dom.init_coral_cover: (sizes, species, locations)
icc_rme = dom.init_coral_cover
icc_sz = size(icc_rme)
nspecies = icc_sz[1]
rme_nlocs = icc_sz[2]

icc_lizard = zeros(Float64, nspecies, n_sites)
for s in 1:n_sites
    rme_idx = site_to_rme_idx[s]
    if rme_idx != -1
        icc_lizard[:, s] .= icc_rme[:, rme_idx]
    end
end

icc_file = joinpath(out_dir, "initial_cover.nc")
isfile(icc_file) && rm(icc_file)
nccreate(icc_file, "covers", 
    "species", collect(1:nspecies), 
    "location", collect(1:n_sites))
ncwrite(icc_lizard, icc_file, "covers")

# 8. DataPackage JSON
metadata = Dict(
    "name" => "Lizard_Island_Cluster",
    "version" => "0.8.0",
    "description" => "Lizard Island Domain built from sub-sampled RME data",
    "spatial" => "spatial/lizard_cluster.gpkg",
    "connectivity" => "connectivity/Lizard_Connectivity.csv",
    "initial_cover" => "initial_cover.nc",
    "dhw" => "DHWs/dhw_RCP45.nc"
)
open(joinpath(out_dir, "datapackage.json"), "w") do f
    write(f, JSON.json(metadata, 4))
end

println("Successfully built Lizard Island domain at: ", out_dir)
