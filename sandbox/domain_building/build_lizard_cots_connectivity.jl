using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using CSV
using DataFrames
using SHA
using TOML
using ADRIA
import ADRIA.GDF as GDF

const DATA_ROOT = normpath(joinpath(@__DIR__, "..", "data"))
const RME_ROOT = joinpath(DATA_ROOT, "rme_ml_2025_06_05", "data_files")
const LIZARD_ROOT = joinpath(DATA_ROOT, "Lizard_Historical_v0.1")
const ID_PATH = joinpath(RME_ROOT, "id", "id_list_2024_12_01.csv")
const RME_GPKG = joinpath(RME_ROOT, "region", "reefmod_gbr.gpkg")
const LIZARD_GPKG = joinpath(LIZARD_ROOT, "spatial", "lizard_cluster.gpkg")
const COTS_SOURCE_DIR = joinpath(RME_ROOT, "con_csv")
const OUTPUT_DIR = joinpath(LIZARD_ROOT, "cots_connectivity")
const OUTPUT_PATH = joinpath(OUTPUT_DIR, "Lizard_COTS_Connectivity.csv")

const COTS_FILES = sort(filter(
    name -> startswith(name, "CONNECT_COTS_") && endswith(name, ".csv"),
    readdir(COTS_SOURCE_DIR)
))

function sha256_file(path::String)::String
    return bytes2hex(open(SHA.sha256, path))
end

"""
Read only the requested rows and columns from a headerless ReefMod connectivity CSV.

The row and column positions use the order in the ReefMod ID list. Streaming avoids
materialising each 3,806 x 3,806 matrix when only the Lizard reefs are required.
"""
function read_connectivity_subset(
    path::String,
    reef_indices::Vector{Int},
    full_size::Int
)::Matrix{Float64}
    local_row = Dict(global_idx => local_idx for (local_idx, global_idx) in enumerate(reef_indices))
    subset = zeros(Float64, length(reef_indices), length(reef_indices))
    rows_found = 0
    rows_seen = 0

    open(path, "r") do io
        for (row_idx, line) in enumerate(eachline(io))
            rows_seen = row_idx
            local_idx = get(local_row, row_idx, 0)
            local_idx == 0 && continue

            values = split(line, ',')
            length(values) == full_size || error(
                "$(basename(path)) row $row_idx has $(length(values)) columns; expected $full_size"
            )
            for (target_idx, global_col) in enumerate(reef_indices)
                subset[local_idx, target_idx] = parse(Float64, values[global_col])
            end
            rows_found += 1
        end
    end

    rows_seen == full_size || error(
        "$(basename(path)) has $rows_seen rows; expected $full_size"
    )
    rows_found == length(reef_indices) || error(
        "$(basename(path)) supplied $rows_found selected rows; expected $(length(reef_indices))"
    )
    return subset
end

isempty(COTS_FILES) && error("No CONNECT_COTS CSV files found in $COTS_SOURCE_DIR")

id_list = CSV.read(ID_PATH, DataFrame; header=false, comment="#")
full_size = nrow(id_list)
rme_spatial = GDF.read(RME_GPKG)
rme_id_col = "RME_GBRMPA_ID" in names(rme_spatial) ? "RME_GBRMPA_ID" : "LABEL_ID"
rme_order = Dict(string(id) => idx for (idx, id) in enumerate(rme_spatial[!, rme_id_col]))
ordered_spatial_rows = [
    get(rme_order, string(id), 0) for id in id_list[!, 1]
]
any(iszero, ordered_spatial_rows) && error("The ReefMod geopackage does not cover the complete ID list")
rme_spatial = rme_spatial[ordered_spatial_rows, :]

lizard_sites = GDF.read(LIZARD_GPKG)
site_ids = String.(DataFrames.make_unique(Symbol.(string.(lizard_sites.site_id)); makeunique=true))
site_reef_ids = string.(lizard_sites.UNIQUE_ID)
rme_unique_ids = string.(rme_spatial.UNIQUE_ID)
rme_index = Dict(id => idx for (idx, id) in enumerate(rme_unique_ids))
missing_reef_ids = setdiff(unique(site_reef_ids), keys(rme_index))
isempty(missing_reef_ids) || error(
    "Lizard reefs absent from the ReefMod ordering: $(sort(missing_reef_ids))"
)

reef_indices = sort(unique(rme_index[id] for id in site_reef_ids))
reef_local_index = Dict(global_idx => local_idx for (local_idx, global_idx) in enumerate(reef_indices))
site_reef_local = [reef_local_index[rme_index[id]] for id in site_reef_ids]

println("Extracting $(length(reef_indices)) Lizard reefs from $(length(COTS_FILES)) COTS matrices...")
reef_conn = zeros(Float64, length(reef_indices), length(reef_indices))
for file in COTS_FILES
    println("  ", file)
    reef_conn .+= read_connectivity_subset(
        joinpath(COTS_SOURCE_DIR, file), reef_indices, full_size
    )
end
reef_conn ./= length(COTS_FILES)

all(isfinite, reef_conn) || error("COTS connectivity contains non-finite values")
all(0.0 .<= reef_conn .<= 1.0) || error("COTS connectivity falls outside [0, 1]")

n_sites = length(site_ids)
sites_per_reef = zeros(Int, length(reef_indices))
for reef_idx in site_reef_local
    sites_per_reef[reef_idx] += 1
end

site_conn = Matrix{Float64}(undef, n_sites, n_sites)
for sink in 1:n_sites
    sink_reef = site_reef_local[sink]
    divisor = sites_per_reef[sink_reef]
    for source in 1:n_sites
        source_reef = site_reef_local[source]
        site_conn[source, sink] = reef_conn[source_reef, sink_reef] / divisor
    end
end

# Every site-level source must preserve the original reef-to-reef probability after
# summing across the sites belonging to a sink reef.
for source in 1:n_sites
    source_reef = site_reef_local[source]
    for sink_reef in axes(reef_conn, 2)
        sink_sites = findall(==(sink_reef), site_reef_local)
        isapprox(
            sum(site_conn[source, sink_sites]),
            reef_conn[source_reef, sink_reef];
            atol=1e-12,
            rtol=1e-12
        ) || error("Site expansion failed mass preservation validation")
    end
end

mkpath(OUTPUT_DIR)
output = DataFrame(site_conn, Symbol.(site_ids))
insertcols!(output, 1, :Source => site_ids)
CSV.write(OUTPUT_PATH, output)

metadata = Dict(
    "method" => "arithmetic mean across available spawning seasons; each reef-to-reef sink probability distributed uniformly among sink sites",
    "orientation" => "rows are sources; columns are sinks",
    "full_reef_count" => full_size,
    "selected_reef_count" => length(reef_indices),
    "site_count" => n_sites,
    "id_list" => Dict("path" => relpath(ID_PATH, LIZARD_ROOT), "sha256" => sha256_file(ID_PATH)),
    "source_files" => [
        Dict(
            "name" => file,
            "sha256" => sha256_file(joinpath(COTS_SOURCE_DIR, file))
        ) for file in COTS_FILES
    ],
    "output" => Dict(
        "name" => basename(OUTPUT_PATH),
        "sha256" => sha256_file(OUTPUT_PATH)
    )
)
open(joinpath(OUTPUT_DIR, "provenance.toml"), "w") do io
    TOML.print(io, metadata; sorted=true)
end

println("Wrote $OUTPUT_PATH")
println("Dimensions: $(size(site_conn)); min=$(minimum(site_conn)); max=$(maximum(site_conn))")
