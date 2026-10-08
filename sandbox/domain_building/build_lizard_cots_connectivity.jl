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
const LIZARD_DOMAIN_VERSION = get(ENV, "LIZARD_DOMAIN_VERSION", "Lizard_Historical_v0.1")
occursin(r"^[A-Za-z0-9_.-]+$", LIZARD_DOMAIN_VERSION) || error(
    "LIZARD_DOMAIN_VERSION contains unsupported path characters"
)
const LIZARD_ROOT = joinpath(DATA_ROOT, LIZARD_DOMAIN_VERSION)
const ID_PATH = joinpath(RME_ROOT, "id", "id_list_2024_12_01.csv")
const RME_GPKG = joinpath(RME_ROOT, "region", "reefmod_gbr.gpkg")
const LIZARD_GPKG = joinpath(LIZARD_ROOT, "spatial", "lizard_cluster.gpkg")
const COTS_SOURCE_DIR = joinpath(RME_ROOT, "con_csv")
const WATER_SOURCE_DIR = joinpath(RME_ROOT, "water_csv")
const OUTBREAK_PREDICTION_PATH = joinpath(DATA_ROOT, "gbrPredsAdj_20262408.csv")
const OUTBREAK_CPUE_THRESHOLD = 0.22
const OUTPUT_DIR = joinpath(LIZARD_ROOT, "cots_connectivity")
const OUTPUT_PATH = joinpath(OUTPUT_DIR, "Lizard_COTS_Connectivity.csv")
const FORCING_DIR = joinpath(OUTPUT_DIR, "annual_forcing")
const WATER_SCENARIO = get(ENV, "COTS_WATER_QUALITY_SCENARIO", "q3baseline")
const SURVIVAL_PATH = joinpath(WATER_SOURCE_DIR, "COTS_LARVAL_REDUCTION_$(WATER_SCENARIO).csv")

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

"""
Read one full GBR connectivity matrix and return the internal Lizard block plus
the survival-weighted contribution from all external GBR source reefs to each
Lizard sink. Matrix rows are sources and columns are sinks.
"""
function read_annual_forcing(
    path::String,
    reef_indices::Vector{Int},
    full_size::Int,
    source_survival::Vector{Float64},
    outbreak_probability::Matrix{Float64},
    habitable_area_ha::Vector{Float64}
)::Tuple{Matrix{Float64},Vector{Float64},Matrix{Float64},Matrix{Float64}}
    length(source_survival) == full_size || error("Larval survival length does not match GBR reef count")
    size(outbreak_probability, 1) == full_size ||
        error("Outbreak-probability rows do not match GBR reef count")
    length(habitable_area_ha) == full_size ||
        error("Habitable-area vector does not match GBR reef count")
    local_row = Dict(global_idx => local_idx for (local_idx, global_idx) in enumerate(reef_indices))
    internal = zeros(Float64, length(reef_indices), length(reef_indices))
    external = zeros(Float64, length(reef_indices))
    external_outbreak_potential = zeros(
        Float64, length(reef_indices), size(outbreak_probability, 2)
    )
    external_outbreak_pelagic = zeros(
        Float64, length(reef_indices), size(outbreak_probability, 2)
    )
    rows_seen = 0

    open(path, "r") do io
        for (source_idx, line) in enumerate(eachline(io))
            rows_seen = source_idx
            values = split(line, ',')
            length(values) == full_size || error(
                "$(basename(path)) row $source_idx has $(length(values)) columns; expected $full_size"
            )
            local_source = get(local_row, source_idx, 0)
            for (local_sink, global_sink) in enumerate(reef_indices)
                probability = parse(Float64, values[global_sink])
                if local_source == 0
                    external[local_sink] += probability * source_survival[source_idx]
                    area_ratio = habitable_area_ha[source_idx] /
                        habitable_area_ha[global_sink]
                    pressure = probability * area_ratio
                    source_probability = @view outbreak_probability[source_idx, :]
                    external_outbreak_potential[local_sink, :] .+=
                        pressure .* source_probability
                    external_outbreak_pelagic[local_sink, :] .+=
                        pressure * source_survival[source_idx] .* source_probability
                else
                    internal[local_source, local_sink] = probability
                end
            end
        end
    end
    rows_seen == full_size || error("$(basename(path)) has $rows_seen rows; expected $full_size")
    return (
        internal,
        external,
        external_outbreak_potential,
        external_outbreak_pelagic
    )
end

isempty(COTS_FILES) && error("No CONNECT_COTS CSV files found in $COTS_SOURCE_DIR")

id_list = CSV.read(ID_PATH, DataFrame; header=false, comment="#")
full_size = nrow(id_list)
rme_ids = string.(id_list[!, 1])
habitable_area_ha = Float64.(id_list[!, 2]) .* 100.0 .* (1.0 .- Float64.(id_list[!, 3]))
all(isfinite, habitable_area_ha) && all(habitable_area_ha .> 0.0) ||
    error("RME habitable reef areas must be positive and finite")

isfile(OUTBREAK_PREDICTION_PATH) ||
    error("Missing GBR outbreak-probability predictions: $OUTBREAK_PREDICTION_PATH")
prediction_df = CSV.read(OUTBREAK_PREDICTION_PATH, DataFrame)
prediction_years = sort(unique(Int.(prediction_df.year)))
year_lookup = Dict(year => idx for (idx, year) in enumerate(prediction_years))
rme_id_lookup = Dict(lowercase(id) => idx for (idx, id) in enumerate(rme_ids))
prediction_aliases = Dict(
    "11-325" => "10-441",
    "11-244e" => "11-288",
    "11-244f" => "11-303",
    "11-244g" => "11-310",
    "11-244h" => "11-311",
    "20198" => "20-198"
)

function prediction_reef_id(name::AbstractString)::String
    matched = match(r"\(([^()]*)\)\s*$", name)
    isnothing(matched) && error("Cannot extract reef ID from prediction name: $name")
    return lowercase(strip(matched.captures[1]))
end

outbreak_probability = fill(NaN, full_size, length(prediction_years))
prediction_source_id = fill("", full_size)
prediction_mapping_status = fill("missing", full_size)
ignored_prediction_ids = Set{String}()
for row in eachrow(prediction_df)
    raw_id = prediction_reef_id(string(row.reefName))
    normalized_id = get(prediction_aliases, raw_id, raw_id)
    source_idx = get(rme_id_lookup, normalized_id, 0)
    if source_idx == 0
        push!(ignored_prediction_ids, raw_id)
        continue
    end
    year_idx = year_lookup[Int(row.year)]
    isfinite(outbreak_probability[source_idx, year_idx]) &&
        error("Duplicate prediction for RME reef $(rme_ids[source_idx]) in $(row.year)")
    probability = Float64(row.outbrProb)
    0.0 <= probability <= 1.0 ||
        error("Outbreak probability outside [0, 1] for $(row.reefName) in $(row.year)")
    outbreak_probability[source_idx, year_idx] = probability
    prediction_source_id[source_idx] = raw_id
    prediction_mapping_status[source_idx] =
        raw_id == normalized_id ? "exact" : "documented_alias"
end

missing_prediction_rows = findall(vec(any(isnan.(outbreak_probability); dims=2)))
missing_prediction_ids = rme_ids[missing_prediction_rows]
missing_prediction_ids == ["99-407"] || error(
    "Unexpected RME reefs without complete outbreak predictions: $missing_prediction_ids"
)
outbreak_probability[missing_prediction_rows, :] .= 0.0
prediction_mapping_status[missing_prediction_rows] .= "missing_zero"
all(isfinite, outbreak_probability) || error("Outbreak-probability mapping is incomplete")

prediction_mapping_df = DataFrame(
    rme_reef_id=rme_ids,
    prediction_reef_id=prediction_source_id,
    mapping_status=prediction_mapping_status
)

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
reef_labels = rme_unique_ids[reef_indices]

isfile(SURVIVAL_PATH) || error("Missing COTS larval-survival file: $SURVIVAL_PATH")
survival_df = CSV.read(SURVIVAL_PATH, DataFrame)
survival_lookup = Dict(string(id) => row for (row, id) in enumerate(survival_df.id))
ordered_survival_rows = [get(survival_lookup, string(id), 0) for id in id_list[!, 1]]
any(iszero, ordered_survival_rows) && error("Larval-survival file does not cover the complete ReefMod ID list")

function forcing_year(filename::String)::Int
    matched = match(r"CONNECT_COTS_(\d{4})_\d{2}\.csv", filename)
    isnothing(matched) && error("Cannot infer forcing year from $filename")
    return parse(Int, matched.captures[1])
end

println("Extracting $(length(reef_indices)) Lizard reefs from $(length(COTS_FILES)) COTS matrices...")
reef_conn = zeros(Float64, length(reef_indices), length(reef_indices))
annual_years = Int[]
annual_reef_conn = Matrix{Float64}[]
annual_source_survival = Vector{Float64}[]
annual_external_supply = Vector{Float64}[]
annual_external_outbreak_potential = Matrix{Float64}[]
annual_external_outbreak_pelagic = Matrix{Float64}[]
for file in COTS_FILES
    year = forcing_year(file)
    survival_col = string(year)
    survival_col in names(survival_df) || error("No larval-survival values for $year")
    full_survival = Float64.(survival_df[ordered_survival_rows, survival_col])
    println("  ", file, " with ", basename(SURVIVAL_PATH))
    internal, external, outbreak_potential, outbreak_pelagic = read_annual_forcing(
        joinpath(COTS_SOURCE_DIR, file),
        reef_indices,
        full_size,
        full_survival,
        outbreak_probability,
        habitable_area_ha
    )
    reef_conn .+= internal
    push!(annual_years, year)
    push!(annual_reef_conn, internal)
    push!(annual_source_survival, full_survival[reef_indices])
    push!(annual_external_supply, external)
    push!(annual_external_outbreak_potential, outbreak_potential)
    push!(annual_external_outbreak_pelagic, outbreak_pelagic)
end
reef_conn ./= length(COTS_FILES)

all(isfinite, reef_conn) || error("COTS connectivity contains non-finite values")
all(0.0 .<= reef_conn .<= 1.0) || error("COTS connectivity falls outside [0, 1]")
all(all(isfinite, matrix) for matrix in annual_reef_conn) || error("Annual connectivity contains non-finite values")
all(all(0.0 .<= matrix .<= 1.0) for matrix in annual_reef_conn) || error("Annual connectivity falls outside [0, 1]")
all(all(0.0 .<= values .<= 1.0) for values in annual_source_survival) || error("Larval survival falls outside [0, 1]")
all(all(isfinite, values) && all(values .>= 0.0) for values in annual_external_supply) || error("External supply coefficients are invalid")
all(all(isfinite, values) && all(values .>= 0.0) for values in annual_external_outbreak_potential) ||
    error("Potential outbreak-boundary coefficients are invalid")
all(all(isfinite, values) && all(values .>= 0.0) for values in annual_external_outbreak_pelagic) ||
    error("Pelagic outbreak-boundary coefficients are invalid")

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

mkpath(FORCING_DIR)
site_map_path = joinpath(FORCING_DIR, "site_reef_map.csv")
CSV.write(
    site_map_path,
    DataFrame(
        site_id=site_ids,
        reef_index=site_reef_local,
        reef_id=reef_labels[site_reef_local],
        reef_habitable_area_ha=habitable_area_ha[reef_indices][site_reef_local]
    )
)

matrix_paths = String[]
for (year, matrix) in zip(annual_years, annual_reef_conn)
    matrix_path = joinpath(FORCING_DIR, "reef_connectivity_$(year).csv")
    matrix_df = DataFrame(matrix, Symbol.(reef_labels))
    insertcols!(matrix_df, 1, :Source_Reef => reef_labels)
    CSV.write(matrix_path, matrix_df)
    push!(matrix_paths, matrix_path)
end

source_survival_path = joinpath(FORCING_DIR, "source_larval_survival.csv")
source_survival_df = DataFrame(reef_id=reef_labels)
for (year, values) in zip(annual_years, annual_source_survival)
    source_survival_df[!, Symbol(string(year))] = values
end
CSV.write(source_survival_path, source_survival_df)

external_supply_path = joinpath(FORCING_DIR, "external_supply_coefficients.csv")
external_supply_df = DataFrame(reef_id=reef_labels)
for (year, values) in zip(annual_years, annual_external_supply)
    external_supply_df[!, Symbol(string(year))] = values
end
CSV.write(external_supply_path, external_supply_df)

prediction_mapping_path = joinpath(FORCING_DIR, "outbreak_prediction_mapping.csv")
CSV.write(prediction_mapping_path, prediction_mapping_df)

outbreak_boundary_path = joinpath(FORCING_DIR, "external_outbreak_boundary.csv")
outbreak_boundary_df = DataFrame(
    reef_id=String[],
    calendar_year=Int[],
    hydrodynamic_year=Int[],
    probability_area_connectivity=Float64[],
    survival_weighted_probability_area_connectivity=Float64[]
)
for (
    forcing_idx,
    hydrodynamic_year,
    potential,
    pelagic
) in zip(
    eachindex(annual_years),
    annual_years,
    annual_external_outbreak_potential,
    annual_external_outbreak_pelagic
)
    for (year_idx, calendar_year) in enumerate(prediction_years)
        for reef_idx in eachindex(reef_labels)
            push!(
                outbreak_boundary_df,
                (
                    reef_labels[reef_idx],
                    calendar_year,
                    hydrodynamic_year,
                    potential[reef_idx, year_idx],
                    pelagic[reef_idx, year_idx]
                )
            )
        end
    end
end
CSV.write(outbreak_boundary_path, outbreak_boundary_df)

forcing_outputs = vcat(
    matrix_paths,
    [
        site_map_path,
        source_survival_path,
        external_supply_path,
        prediction_mapping_path,
        outbreak_boundary_path
    ]
)
forcing_metadata = Dict(
    "schema_version" => 2,
    "orientation" => "reef connectivity rows are sources and columns are sinks",
    "state_density_units" => "COTS ha^-1",
    "external_source_units" => "pre-dispersal recruits ha^-1 yr^-1",
    "years" => annual_years,
    "water_quality_scenario" => WATER_SCENARIO,
    "water_quality_source" => Dict(
        "path" => relpath(SURVIVAL_PATH, LIZARD_ROOT),
        "sha256" => sha256_file(SURVIVAL_PATH)
    ),
    "external_supply_definition" => "sum of external-GBR source-to-Lizard connectivity multiplied by source larval survival; multiply by absolute external pre-dispersal recruit production in recruits ha^-1 yr^-1",
    "site_map" => basename(site_map_path),
    "source_survival" => basename(source_survival_path),
    "external_supply" => basename(external_supply_path),
    "external_outbreak_boundary" => basename(outbreak_boundary_path),
    "outbreak_prediction_mapping" => basename(prediction_mapping_path),
    "matrix_files" => basename.(matrix_paths),
    "outbreak_probability_source" => Dict(
        "path" => relpath(OUTBREAK_PREDICTION_PATH, LIZARD_ROOT),
        "sha256" => sha256_file(OUTBREAK_PREDICTION_PATH),
        "calendar_years" => prediction_years,
        "event_definition" => "probability COTS CPUE exceeds 0.22 COTS per tow",
        "cpue_threshold_cots_per_tow" => OUTBREAK_CPUE_THRESHOLD,
        "interpretation" => "probability is used continuously as expected upstream outbreak activity, not thresholded as a probability",
        "documented_alias_count" => count(==("documented_alias"), prediction_mapping_status),
        "missing_zero_ids" => missing_prediction_ids,
        "ignored_prediction_id_count" => length(ignored_prediction_ids)
    ),
    "outbreak_boundary_definition" => "sum over external GBR sources of outbreak probability multiplied by source-to-sink connectivity and source-habitable-area/sink-habitable-area; pelagic column additionally applies source larval survival",
    "outbreak_boundary_units" => "dimensionless multiplier of pre-pelagic outbreak production density in COTS ha^-1 yr^-1",
    "outputs" => [
        Dict("name" => basename(path), "sha256" => sha256_file(path)) for path in forcing_outputs
    ]
)
open(joinpath(FORCING_DIR, "provenance.toml"), "w") do io
    TOML.print(io, forcing_metadata; sorted=true)
end

metadata = Dict(
    "method" => "static control: arithmetic mean across available spawning seasons; each reef-to-reef sink probability distributed uniformly among sink sites",
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
println("Wrote annual reef forcing and external-boundary coefficients to $FORCING_DIR")
