"""
Build a site x year background coral mortality table following Mellin et al. (2019) /
COTS_Model (crash_mechanism_protocol.md, 2026-10-13).

Disease and unknown disturbances act on log coral cover: in an event year
log(C) += B_disease (or B_unknown), so mortality m = 1 - exp(B_d * disease + B_u * unknown).
Event records are used where LTMP observed the reef:
  1. CoralBlox-params-calibration `datasets/ltmp_data/disturbances.nc` (2008-2022), then
  2. the COTS_Model LTMP tables exported by export_mellin_disturbance.R (1996-2017).
Other reef-years take the expected loss at the COTS_Model simulated incidence (1% per
cell-year for each of disease and unknown): m = 1 - (1 - p(1 - e^B_d))(1 - p(1 - e^B_u)).

    julia --project=sandbox sandbox/calibration/build_mellin_background_mortality.jl <export dir> <output dir>
"""

using ADRIA, CSV, DataFrames, SHA, TOML, Dates

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const DOMAIN = joinpath(ROOT, "sandbox", "data", "Lizard_Historical_v0.3")
const DISTURBANCES_NC = joinpath(ROOT, "..", "CoralBlox-params-calibration", "datasets",
    "ltmp_data", "disturbances.nc")
const YEARS = 1991:2024
const SIMULATED_INCIDENCE = 0.01

file_sha256(path) = bytes2hex(open(sha256, path))

"""LTMP reef code (e.g. 14114100104) to the GBRMPA reef base (14-114)."""
ltmp_base(code::AbstractString) = string(code[1:2], "-", code[3:5])
gbrmpa_base(id::AbstractString) = replace(id, r"[a-z]$" => "")

function main(export_dir::String, output::String)
    ispath(output) && error("Refusing to overwrite $output")
    coefficients = CSV.read(joinpath(export_dir, "mellin_coefficients.csv"), DataFrame)
    B = Dict(String(r.term) => Float64(r.log_cover_effect_mean) for r in eachrow(coefficients))
    B_d, B_u = B["disease"], B["unknown"]
    expected = 1.0 - (1.0 - SIMULATED_INCIDENCE * (1.0 - exp(B_d))) *
        (1.0 - SIMULATED_INCIDENCE * (1.0 - exp(B_u)))

    # Source 2: COTS_Model LTMP tables (reef_id exact or base), 1996-2017.
    mellin = CSV.read(joinpath(export_dir, "mellin_ltmp_events.csv"), DataFrame;
        types=Dict(:reef_id => String))
    mellin_events = Dict{Tuple{String,Int,String},Int}()
    for r in eachrow(mellin)
        key = (gbrmpa_base(r.reef_id), Int(r.year), String(r.kind))
        mellin_events[key] = max(get(mellin_events, key, 0), Int(r.event))
    end

    # Source 1: CoralBlox LTMP disturbance flags, 2008-2022.
    ds = ADRIA.YAXArrays.open_dataset(DISTURBANCES_NC)
    layer = ds.layer
    DD = ADRIA.YAXArrays.DD
    nc_locs = String.(collect(DD.lookup(layer, :locations)))
    nc_years = Int.(collect(DD.lookup(layer, :timesteps)))
    nc_kinds = String.(collect(DD.lookup(layer, :disturbances)))
    flags = Array(layer)
    nc_events = Dict{Tuple{String,Int,String},Int}()
    for (i, loc) in enumerate(nc_locs), (j, year) in enumerate(nc_years)
        for (kind, label) in (("disease", "Disease"), ("unknown", "Unknown"))
            k = findfirst(==(label), nc_kinds)
            key = (ltmp_base(loc), year, kind)
            nc_events[key] = max(get(nc_events, key, 0), Int(flags[i, j, k] != 0))
        end
    end

    sites = ADRIA.GDF.read(joinpath(DOMAIN, "spatial", "lizard_cluster.gpkg"))
    site_ids = String.(string.(sites.site_id))
    reefs = gbrmpa_base.(String.(sites.GBRMPA_ID))
    mortality = zeros(length(site_ids), length(YEARS))
    source = fill("expected", length(site_ids), length(YEARS))
    for (i, reef) in enumerate(reefs), (j, year) in enumerate(YEARS)
        record = nothing
        for (events, label) in ((nc_events, "ltmp_coralblox_nc"), (mellin_events, "ltmp_cots_model"))
            if haskey(events, (reef, year, "disease")) || haskey(events, (reef, year, "unknown"))
                record = (get(events, (reef, year, "disease"), 0), get(events, (reef, year, "unknown"), 0), label)
                break
            end
        end
        if isnothing(record)
            mortality[i, j] = expected
        else
            mortality[i, j] = 1.0 - exp(B_d * record[1] + B_u * record[2])
            source[i, j] = record[3]
        end
    end
    all(0.0 .<= mortality .< 1.0) || error("Background mortality outside [0, 1)")

    mkpath(output)
    table = DataFrame(site_id=site_ids)
    for (j, year) in enumerate(YEARS)
        table[!, string(year)] = mortality[:, j]
    end
    CSV.write(joinpath(output, "background_mortality.csv"), table)
    event_rows = [(site_ids[i], reefs[i], YEARS[j], mortality[i, j], source[i, j])
        for i in eachindex(site_ids), j in eachindex(YEARS) if source[i, j] != "expected" && mortality[i, j] > 0]
    CSV.write(joinpath(output, "recorded_events.csv"),
        DataFrame(event_rows, [:site_id, :reef, :year, :mortality, :source]) |> unique)
    metadata = Dict(
        "status" => "experimental_input_not_promoted",
        "method" => "Mellin et al. 2019 / COTS_Model: m = 1 - exp(B_d*disease + B_u*unknown); expected value at 1% simulated incidence where LTMP did not observe the reef-year",
        "B_disease_log_cover" => B_d, "B_unknown_log_cover" => B_u,
        "simulated_incidence_per_cell_year" => SIMULATED_INCIDENCE,
        "expected_annual_mortality" => expected,
        "years" => [first(YEARS), last(YEARS)],
        "site_years_recorded" => count(!=("expected"), source),
        "site_years_with_recorded_event" => count(i -> source[i] != "expected" && mortality[i] > 0, eachindex(source)),
        "inputs" => Dict(
            "mellin_coefficients.csv" => file_sha256(joinpath(export_dir, "mellin_coefficients.csv")),
            "mellin_ltmp_events.csv" => file_sha256(joinpath(export_dir, "mellin_ltmp_events.csv")),
            "mellin_provenance.txt" => read(joinpath(export_dir, "provenance.txt"), String),
            "disturbances_nc" => DISTURBANCES_NC, "disturbances_nc_sha256" => file_sha256(DISTURBANCES_NC),
            "domain_gpkg_sha256" => file_sha256(joinpath(DOMAIN, "spatial", "lizard_cluster.gpkg")),
            "builder_sha256" => file_sha256(@__FILE__)),
        "created" => string(now()),
    )
    open(joinpath(output, "metadata.toml"), "w") do io
        TOML.print(io, metadata)
    end
    println("expected annual mortality ", round(expected; digits=5), "; recorded site-years ",
        metadata["site_years_recorded"], "; with events ", metadata["site_years_with_recorded_event"])
    events = DataFrame(event_rows, [:site_id, :reef, :year, :mortality, :source])
    show(unique(events[:, [:reef, :year, :mortality, :source]]), allrows=true)
end

main(ARGS[1], ARGS[2])
