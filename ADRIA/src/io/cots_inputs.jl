# Optional COTS inputs of a data package; see the `CotsInputs` docstring for the layout.
"""
Name of the data-package folder holding the optional COTS inputs.
"""
const COTS_CONN_DIR = "cots_connectivity"
const COTS_FORCING_DIR = "annual_forcing"
const COTS_INITIAL_STATE_FN = "cots_initial_state.csv"
const COTS_LARVAL_SURVIVAL_FN = "cots_larval_survival.csv"
const COTS_SETTLEMENT_MODIFIER_FN = "cots_settlement_modifier.csv"

"""
    CotsReefForcing

Annual reef-scale COTS larval connectivity, and the larval supply arriving from outbreaks
outside the domain. Required by separated recruitment and by the annual connectivity modes.

# Fields
- `years` : hydrodynamic years, one per reef connectivity matrix
- `reef_connectivity` : reef x reef source-to-sink connectivity, one matrix per year
- `site_to_reef` : reef index of each domain location
- `reef_habitable_area_ha` : habitable area of each reef (ha)
- `outbreak_years` : calendar years with outside outbreak supply (empty if not provided)
- `outbreak_supply` : outside supply coefficients, reef x outbreak year x hydrodynamic year
"""
struct CotsReefForcing
    years::Vector{Int}
    reef_connectivity::Vector{Matrix{Float64}}
    site_to_reef::Vector{Int}
    reef_habitable_area_ha::Vector{Float64}
    outbreak_years::Vector{Int}
    outbreak_supply::Array{Float64,3}
end

"""
    CotsInputs

The optional COTS inputs of a domain. Each field is `nothing` when the data package does
not provide it, and a data package without any of them runs COTS with larvae dispersing
through the coral connectivity and COTSMod initialising the populations.

# Data package layout

```
<data package>/
    cots_connectivity/                    site-to-site COTS larval connectivity (CSV)
        annual_forcing/                   reef-scale connectivity, one matrix per year
            provenance.toml               `years`, `site_map`, `matrix_files` and,
                                          optionally, `external_outbreak_boundary`
            <site_map>.csv                `site_id`, `reef_index`, `reef_habitable_area_ha`
            <matrix_files>.csv            reef x reef, first column holds the reef ids
            <external_outbreak_boundary>.csv
                                          `reef_id`, `calendar_year`, `hydrodynamic_year`,
                                          `probability_area_connectivity`
    cots_initial_state.csv                `site_id`, `juvenile_ha`, `subadult_ha`, `adult_ha`
    cots_larval_survival.csv              `site_id`, `year`, `multiplier`
    cots_settlement_modifier.csv          `site_id`, `multiplier`
```

`site_id` values are the domain's location ids. Densities are COTS ha^-1; multipliers are
dimensionless. Years in `cots_larval_survival.csv` and `calendar_year` in the outbreak
boundary are calendar years, matched against the domain's timeframe.

# Fields
- `conn` : site-to-site COTS larval connectivity
- `reef_forcing` : annual reef-scale connectivity and outside outbreak supply
- `initial_state` : starting COTS density (COTS ha^-1), locations x (juvenile, subadult,
    adult)
- `larval_survival` : larval survival multipliers, timesteps x locations
- `settlement_modifier` : settlement multipliers in [0, 1], one per location
"""
Base.@kwdef struct CotsInputs
    conn::Union{Nothing,YAXArray} = nothing
    reef_forcing::Union{Nothing,CotsReefForcing} = nothing
    initial_state::Union{Nothing,Matrix{Float64}} = nothing
    larval_survival::Union{Nothing,Matrix{Float64}} = nothing
    settlement_modifier::Union{Nothing,Vector{Float64}} = nothing
end

"""
    load_cots_inputs(dpkg_path::String, loc_ids::Vector{String}, timeframe::Vector{Int64})::CotsInputs

Load the optional COTS inputs of the data package at `dpkg_path` for the locations
`loc_ids` (in domain order) and calendar years `timeframe`.
"""
function load_cots_inputs(
    dpkg_path::String, loc_ids::Vector{String}, timeframe::Vector{Int64}
)::CotsInputs
    conn_dir = joinpath(dpkg_path, COTS_CONN_DIR)
    conn = if isdir(conn_dir) && any(endswith(".csv"), readdir(conn_dir))
        loaded = location_connectivity(conn_dir, loc_ids)
        isempty(loaded.truncated) ||
            error("COTS connectivity lacks domain locations: $(loaded.truncated)")
        loaded.conn
    else
        nothing
    end

    forcing_dir = joinpath(conn_dir, COTS_FORCING_DIR)
    reef_forcing = if isfile(joinpath(forcing_dir, "provenance.toml"))
        load_cots_reef_forcing(forcing_dir, loc_ids)
    else
        nothing
    end

    path(fn) = joinpath(dpkg_path, fn)
    return CotsInputs(;
        conn=conn,
        reef_forcing=reef_forcing,
        initial_state=isfile(path(COTS_INITIAL_STATE_FN)) ?
                      load_cots_initial_state(path(COTS_INITIAL_STATE_FN), loc_ids) :
                      nothing,
        larval_survival=isfile(path(COTS_LARVAL_SURVIVAL_FN)) ?
                        load_cots_larval_survival(
            path(COTS_LARVAL_SURVIVAL_FN), loc_ids, timeframe
        ) : nothing,
        settlement_modifier=isfile(path(COTS_SETTLEMENT_MODIFIER_FN)) ?
                            load_cots_settlement_modifier(
            path(COTS_SETTLEMENT_MODIFIER_FN), loc_ids
        ) : nothing
    )
end

"""Row of each domain location in `ids`; errors if any location is missing."""
function _cots_site_rows(
    ids::AbstractVector, loc_ids::Vector{String}, what::String
)::Vector{Int}
    lookup = Dict(string(id) => row for (row, id) in enumerate(ids))
    rows = [get(lookup, id, 0) for id in loc_ids]
    any(iszero, rows) && error("$what does not cover every domain location")
    return rows
end

"""
    load_cots_reef_forcing(forcing_dir::String, loc_ids::Vector{String})::CotsReefForcing

Read annual reef-scale COTS connectivity, as described by `provenance.toml` in
`forcing_dir`.
"""
function load_cots_reef_forcing(
    forcing_dir::String, loc_ids::Vector{String}
)::CotsReefForcing
    metadata = TOML.parsefile(joinpath(forcing_dir, "provenance.toml"))
    years = Int.(metadata["years"])

    site_map = CSV.read(joinpath(forcing_dir, metadata["site_map"]), DataFrame)
    rows = _cots_site_rows(site_map.site_id, loc_ids, "COTS reef site map")
    site_to_reef = Int.(site_map.reef_index[rows])

    reef_connectivity = Matrix{Float64}[]
    reef_order = String[]
    for fn in metadata["matrix_files"]
        matrix_df = CSV.read(joinpath(forcing_dir, fn), DataFrame)
        order = string.(matrix_df[!, 1])
        isempty(reef_order) && (reef_order = order)
        order == reef_order || error("Annual COTS reef matrices use different reef orders")
        push!(reef_connectivity, Matrix{Float64}(matrix_df[:, 2:end]))
    end
    length(reef_connectivity) == length(years) ||
        error("Expected one COTS reef matrix per year in provenance.toml")
    n_reefs = length(reef_order)
    all(size(m) == (n_reefs, n_reefs) for m in reef_connectivity) ||
        error("Annual COTS reef matrices have invalid dimensions")

    reef_habitable_area_ha = [
        only(
            unique(Float64.(site_map.reef_habitable_area_ha[site_map.reef_index .== reef]))
        )
        for reef in 1:n_reefs
    ]
    all(isfinite, reef_habitable_area_ha) && all(reef_habitable_area_ha .> 0.0) ||
        error("COTS reef habitable areas must be positive and finite")

    outbreak_years = Int[]
    outbreak_supply = zeros(Float64, n_reefs, 0, length(years))
    if haskey(metadata, "external_outbreak_boundary")
        boundary = CSV.read(
            joinpath(forcing_dir, metadata["external_outbreak_boundary"]), DataFrame
        )
        sort(unique(Int.(boundary.hydrodynamic_year))) == sort(years) ||
            error("Outbreak boundary hydrodynamic years do not match the reef matrices")
        outbreak_years = sort(unique(Int.(boundary.calendar_year)))
        reef_idx = Dict(id => i for (i, id) in enumerate(reef_order))
        year_idx = Dict(y => i for (i, y) in enumerate(outbreak_years))
        hydro_idx = Dict(y => i for (i, y) in enumerate(years))
        outbreak_supply = fill(NaN, n_reefs, length(outbreak_years), length(years))
        for row in eachrow(boundary)
            r = get(reef_idx, string(row.reef_id), 0)
            r > 0 || error("Outbreak boundary contains an unknown reef: $(row.reef_id)")
            outbreak_supply[r, year_idx[Int(row.calendar_year)], hydro_idx[Int(row.hydrodynamic_year)]] = Float64(
                row.probability_area_connectivity
            )
        end
        all(isfinite, outbreak_supply) || error("Outbreak boundary is incomplete")
    end

    return CotsReefForcing(
        years, reef_connectivity, site_to_reef, reef_habitable_area_ha, outbreak_years,
        outbreak_supply
    )
end

"""
    load_cots_initial_state(path::String, loc_ids::Vector{String})::Matrix{Float64}

Starting COTS density (COTS ha^-1) by location (rows, in domain order) and stage
(juvenile, subadult, adult).
"""
function load_cots_initial_state(path::String, loc_ids::Vector{String})::Matrix{Float64}
    table = CSV.read(path, DataFrame; types=Dict(:site_id => String))
    columns = ["juvenile_ha", "subadult_ha", "adult_ha"]
    all(in(names(table)), ["site_id"; columns]) ||
        error("$path needs columns site_id, $(join(columns, ", "))")
    allunique(table.site_id) || error("Duplicate locations in $path")
    rows = _cots_site_rows(table.site_id, loc_ids, path)
    state = Matrix{Float64}(table[rows, columns])
    all(isfinite, state) && all(>=(0.0), state) ||
        error("COTS starting densities must be finite and non-negative")
    return state
end

"""
    load_cots_larval_survival(path::String, loc_ids::Vector{String}, timeframe::Vector{Int64})::Matrix{Float64}

Larval survival multipliers by timestep and location. Location-years absent from the table
are 1; years outside `timeframe` are ignored.
"""
function load_cots_larval_survival(
    path::String, loc_ids::Vector{String}, timeframe::Vector{Int64}
)::Matrix{Float64}
    table = CSV.read(path, DataFrame; types=Dict(:site_id => String))
    all(in(names(table)), ["site_id", "year", "multiplier"]) ||
        error("$path needs columns site_id, year, multiplier")
    sites = Dict(id => loc for (loc, id) in enumerate(loc_ids))
    steps = Dict(year => t for (t, year) in enumerate(timeframe))
    multipliers = ones(Float64, length(timeframe), length(loc_ids))
    for row in eachrow(table)
        loc = get(sites, row.site_id, 0)
        loc > 0 ||
            error("COTS larval-survival location $(row.site_id) is not in the domain")
        value = Float64(row.multiplier)
        isfinite(value) && value >= 0.0 ||
            error("COTS larval-survival multipliers must be finite and non-negative")
        t = get(steps, Int(row.year), 0)
        t > 0 && (multipliers[t, loc] = value)
    end
    return multipliers
end

"""
    load_cots_settlement_modifier(path::String, loc_ids::Vector{String})::Vector{Float64}

Settlement multiplier in [0, 1] for each location, constant over time. Every location must
appear exactly once.
"""
function load_cots_settlement_modifier(
    path::String, loc_ids::Vector{String}
)::Vector{Float64}
    table = CSV.read(path, DataFrame; types=Dict(:site_id => String))
    all(in(names(table)), ["site_id", "multiplier"]) ||
        error("$path needs columns site_id, multiplier")
    allunique(table.site_id) || error("Duplicate locations in $path")
    unknown = setdiff(table.site_id, loc_ids)
    isempty(unknown) || error("$path has locations not in the domain: $(unknown)")
    rows = _cots_site_rows(table.site_id, loc_ids, path)
    multipliers = Float64.(table.multiplier[rows])
    all(m -> isfinite(m) && 0.0 <= m <= 1.0, multipliers) ||
        error("COTS settlement modifiers must be finite and in [0, 1]")
    return multipliers
end

"""
    larval_conn_index(forcing::CotsReefForcing, timestep::Int, mode::Int, seed::Int)::Int

Index of the annual reef matrix used at `timestep`: in turn (`mode` 1) or drawn with a
generator seeded by `seed + timestep` (`mode` 2), so the draws do not depend on other
random numbers in the scenario.
"""
function larval_conn_index(
    forcing::CotsReefForcing, timestep::Int, mode::Int, seed::Int
)::Int
    mode == 1 && return mod1(timestep, length(forcing.years))
    mode == 2 && return rand(MersenneTwister(seed + timestep), eachindex(forcing.years))
    throw(ArgumentError("Larval connectivity mode must be 1 or 2 for annual matrices"))
end

"""
    outbreak_supply_by_location(forcing, forcing_idx, calendar_year, production)::Vector{Float64}

Larvae arriving at each location from outbreaks outside the domain in `calendar_year`:
`production` times the reef's outbreak-boundary coefficient for that year, under annual
matrix `forcing_idx` (0 = the mean over the annual matrices). Zero in years the boundary
does not cover.
"""
function outbreak_supply_by_location(
    forcing::CotsReefForcing, forcing_idx::Int, calendar_year::Int, production::Float64
)::Vector{Float64}
    year_idx = findfirst(==(calendar_year), forcing.outbreak_years)
    isnothing(year_idx) && return zeros(Float64, length(forcing.site_to_reef))
    coefficients = if forcing_idx == 0
        vec(mean(@view(forcing.outbreak_supply[:, year_idx, :]); dims=2))
    else
        @view(forcing.outbreak_supply[:, year_idx, forcing_idx])
    end
    return [production * coefficients[reef] for reef in forcing.site_to_reef]
end
