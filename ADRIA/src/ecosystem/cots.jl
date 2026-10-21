# ADRIA compatibility adapter for the external COTSMod package.
#
# ADRIA owns parameter-table sampling, domain wiring, coral-cover tensor conversion,
# result logging, and scenario timing. COTSMod owns the COTS state transition,
# predation, larval dispersal, initialization, and external-supply mechanics.

const AbstractCotsModel = COTSMod.AbstractCOTSModel
const CotsPreyMap = COTSMod.COTSPreyMap
const CotsPreyState = COTSMod.COTSPreyState
const CotsHuman = COTSMod.COTSHuman
const CotsState = COTSMod.COTSState
# Last projection-start reset applied by run_scenario (protocol 2026-10-21, M2); nothing when off.
const COTS_RESET_LOG = Ref{Any}(nothing)

struct COTSConnectivityForcing
    years::Vector{Int}
    reef_connectivity::Vector{Matrix{Float64}}
    site_to_reef::Vector{Int}
    reef_habitable_area_ha::Vector{Float64}
    source_survival::Matrix{Float64}
    external_supply::Matrix{Float64}
    outbreak_years::Vector{Int}
    external_outbreak_potential::Array{Float64,3}
    external_outbreak_pelagic::Array{Float64,3}
end

function load_cots_connectivity_forcing(
    forcing_dir::AbstractString,
    loc_ids::AbstractVector{<:AbstractString}
)::COTSConnectivityForcing
    metadata_path = joinpath(forcing_dir, "provenance.toml")
    isfile(metadata_path) || error("Missing annual COTS forcing metadata: $metadata_path")
    metadata = TOML.parsefile(metadata_path)
    years = Int.(metadata["years"])

    site_map = CSV.read(joinpath(forcing_dir, metadata["site_map"]), DataFrame)
    site_lookup = Dict(string(id) => row for (row, id) in enumerate(site_map.site_id))
    ordered_rows = [get(site_lookup, string(id), 0) for id in loc_ids]
    any(iszero, ordered_rows) && error("Annual COTS forcing does not cover every domain location")
    site_to_reef = Int.(site_map.reef_index[ordered_rows])

    reef_connectivity = Matrix{Float64}[]
    reef_order = String[]
    for filename in metadata["matrix_files"]
        matrix_df = CSV.read(joinpath(forcing_dir, filename), DataFrame)
        current_order = string.(matrix_df[!, 1])
        if isempty(reef_order)
            reef_order = current_order
        else
            current_order == reef_order || error("Annual COTS matrices use different reef orders")
        end
        push!(reef_connectivity, Matrix{Float64}(matrix_df[:, 2:end]))
    end

    function ordered_forcing_matrix(filename::AbstractString)::Matrix{Float64}
        table = CSV.read(joinpath(forcing_dir, filename), DataFrame)
        lookup = Dict(string(id) => row for (row, id) in enumerate(table[!, 1]))
        rows = [get(lookup, id, 0) for id in reef_order]
        any(iszero, rows) && error("Annual COTS forcing table has a different reef set")
        return Matrix{Float64}(table[rows, Symbol.(string.(years))])
    end

    source_survival = ordered_forcing_matrix(metadata["source_survival"])
    external_supply = ordered_forcing_matrix(metadata["external_supply"])
    n_reefs = length(reef_order)
    all(size(matrix) == (n_reefs, n_reefs) for matrix in reef_connectivity) ||
        error("Annual COTS reef matrices have invalid dimensions")
    size(source_survival) == (n_reefs, length(years)) ||
        error("COTS larval-survival forcing has invalid dimensions")
    size(external_supply) == (n_reefs, length(years)) ||
        error("COTS external-supply forcing has invalid dimensions")
    "reef_habitable_area_ha" in names(site_map) ||
        error("Annual COTS forcing site map lacks reef_habitable_area_ha; rebuild it")
    reef_habitable_area_ha = [
        only(unique(Float64.(site_map.reef_habitable_area_ha[site_map.reef_index .== reef])))
        for reef in 1:n_reefs
    ]
    all(isfinite, reef_habitable_area_ha) && all(reef_habitable_area_ha .> 0.0) ||
        error("COTS forcing reef areas must be positive and finite")

    haskey(metadata, "external_outbreak_boundary") ||
        error("Annual COTS forcing lacks evidence-derived outbreak boundary; rebuild it")
    boundary = CSV.read(
        joinpath(forcing_dir, metadata["external_outbreak_boundary"]), DataFrame
    )
    outbreak_years = sort(unique(Int.(boundary.calendar_year)))
    boundary_hydrodynamic_years = sort(unique(Int.(boundary.hydrodynamic_year)))
    boundary_hydrodynamic_years == sort(years) ||
        error("Outbreak boundary hydrodynamic years do not match connectivity forcing")
    reef_lookup = Dict(id => idx for (idx, id) in enumerate(reef_order))
    year_lookup = Dict(year => idx for (idx, year) in enumerate(outbreak_years))
    hydro_lookup = Dict(year => idx for (idx, year) in enumerate(years))
    external_outbreak_potential = fill(
        NaN, n_reefs, length(outbreak_years), length(years)
    )
    external_outbreak_pelagic = similar(external_outbreak_potential)
    for row in eachrow(boundary)
        reef_idx = get(reef_lookup, string(row.reef_id), 0)
        reef_idx > 0 || error("Outbreak boundary contains an unknown sink reef")
        year_idx = year_lookup[Int(row.calendar_year)]
        hydro_idx = hydro_lookup[Int(row.hydrodynamic_year)]
        external_outbreak_potential[reef_idx, year_idx, hydro_idx] =
            Float64(row.probability_area_connectivity)
        external_outbreak_pelagic[reef_idx, year_idx, hydro_idx] =
            Float64(row.survival_weighted_probability_area_connectivity)
    end
    all(isfinite, external_outbreak_potential) ||
        error("Potential outbreak boundary is incomplete")
    all(isfinite, external_outbreak_pelagic) ||
        error("Pelagic outbreak boundary is incomplete")
    return COTSConnectivityForcing(
        years,
        reef_connectivity,
        site_to_reef,
        reef_habitable_area_ha,
        source_survival,
        external_supply,
        outbreak_years,
        external_outbreak_potential,
        external_outbreak_pelagic
    )
end

function cots_forcing_index(
    forcing::COTSConnectivityForcing,
    timestep::Int,
    mode::AbstractString,
    seed::Int
)::Int
    normalized = lowercase(mode)
    normalized == "cycle" && return mod1(timestep, length(forcing.years))
    if normalized == "sample"
        return rand(MersenneTwister(seed + timestep), eachindex(forcing.years))
    end
    error("COTS temporal connectivity mode must be mean, cycle, or sample; got $mode")
end

# Diagnostic intervention only: suppress settlement from the modeled reef
# network in selected calendar years. With no interval, the legacy scalar is
# returned unchanged. External boundary supply is applied separately.
function cots_internal_settlement_scalar(
    base_scalar::Float64,
    calendar_year::Int,
    blackout_years::Union{Nothing,UnitRange{Int}}
)::Float64
    return !isnothing(blackout_years) && calendar_year in blackout_years ?
        0.0 : base_scalar
end

function cots_external_outbreak_supply_by_site(
    forcing::COTSConnectivityForcing,
    forcing_index::Int,
    calendar_year::Int,
    outbreak_production::Float64
)::Tuple{Vector{Float64},Vector{Float64}}
    year_index = findfirst(==(calendar_year), forcing.outbreak_years)
    isnothing(year_index) && return (
        zeros(Float64, length(forcing.site_to_reef)),
        zeros(Float64, length(forcing.site_to_reef))
    )
    potential_coefficients = if forcing_index == 0
        vec(mean(
            @view(forcing.external_outbreak_potential[:, year_index, :]);
            dims=2
        ))
    else
        @view(forcing.external_outbreak_potential[:, year_index, forcing_index])
    end
    pelagic_coefficients = if forcing_index == 0
        vec(mean(
            @view(forcing.external_outbreak_pelagic[:, year_index, :]);
            dims=2
        ))
    else
        @view(forcing.external_outbreak_pelagic[:, year_index, forcing_index])
    end
    potential = outbreak_production .* potential_coefficients
    pelagic = outbreak_production .* pelagic_coefficients
    return (
        [potential[reef] for reef in forcing.site_to_reef],
        [pelagic[reef] for reef in forcing.site_to_reef]
    )
end

function cots_external_supply_by_site(
    forcing::COTSConnectivityForcing,
    forcing_index::Int,
    source_recruit_production::Float64
)::Vector{Float64}
    sites_per_reef = [count(==(reef), forcing.site_to_reef) for reef in axes(forcing.external_supply, 1)]
    return [
        source_recruit_production * forcing.external_supply[reef, forcing_index] / sites_per_reef[reef]
        for reef in forcing.site_to_reef
    ]
end

function _cots_runtime_params(params::COTSMod.COTSParams)::COTSMod.COTSParams
    return params
end

function _cots_runtime_params(params::NamedTuple)::COTSMod.COTSParams
    return COTSMod.COTSParams(
        a=Float64(params.a),
        b=Float64(params.b),
        IMM=Float64(params.IMM),
        p_tilde=Float64(params.p_tilde),
        C_max=Float64(params.C_max),
        m1=Float64(params.m1),
        m2=Float64(params.m2),
        m3=Float64(params.m3),
        a_F=Float64(params.a_F),
        a_S=Float64(params.a_S),
        h=Float64(get(params, :h, 0.0)),
        eta_F=Float64(get(params, :eta_F, 1.0)),
        eta_S=Float64(get(params, :eta_S, 1.0)),
        eta_starve=Float64(get(params, :eta_starve, 2.0)),
        eta_imm=Float64(get(params, :eta_imm, 2.0)),
        imm_threshold=Float64(get(params, :imm_threshold, 0.35)),
        fecundity_gate=Bool(get(params, :fecundity_gate, false)),
        a_ricker=Float64(get(params, :a_ricker, 6.0)),
        b_ricker=Float64(get(params, :b_ricker, 0.1)),
        tau_condition=Float64(get(params, :tau_condition, 5.0)),
        allee_threshold=Float64(get(params, :allee_threshold, 1.0)),
        separated_recruitment=Bool(get(params, :separated_recruitment, false)),
        larval_fecundity=Float64(get(params, :larval_fecundity, 6.0)),
        settlement_probability=Float64(get(params, :settlement_probability, 1.0)),
        juvenile_storage=Bool(get(params, :juvenile_storage, false)),
        juvenile_maturation_min=Float64(get(params, :juvenile_maturation_min, 0.05)),
        juvenile_maturation_max=Float64(get(params, :juvenile_maturation_max, 1.0)),
        juvenile_maturation_cover=Float64(get(params, :juvenile_maturation_cover, 0.4)),
        habitat_mediation=Bool(get(params, :habitat_mediation, false)),
        settlement_floor=Float64(get(params, :settlement_floor, 0.05)),
        low_cover_mortality=Bool(get(params, :low_cover_mortality, false)),
        low_cover_threshold=Float64(get(params, :low_cover_threshold, 0.10)),
        low_cover_strength=Float64(get(params, :low_cover_strength, 0.0)),
        size_weighted_fecundity=Bool(get(params, :size_weighted_fecundity, false)),
        size_small_diameter_mm=Float64(get(params, :size_small_diameter_mm, 200.0)),
        size_large_diameter_mm=Float64(get(params, :size_large_diameter_mm, 300.0)),
        size_fecundity_slope_per_mm=Float64(get(params, :size_fecundity_slope_per_mm, 0.0115)),
        size_initial_large_fraction=Float64(get(params, :size_initial_large_fraction, 0.5)),
        size_growth_max=Float64(get(params, :size_growth_max, 0.5)),
        size_growth_cover=Float64(get(params, :size_growth_cover, 0.4)),
        size_growth_food_mediated=Bool(get(params, :size_growth_food_mediated, false)),
        lagged_adult_hazard=Bool(get(params, :lagged_adult_hazard, false)),
        hazard_food_coupled=Bool(get(params, :hazard_food_coupled, false)),
        hazard_rho=Float64(get(params, :hazard_rho, 0.5)),
        hazard_threshold_ha=Float64(get(params, :hazard_threshold_ha, 1.5)),
        hazard_width_ha=Float64(get(params, :hazard_width_ha, 0.5)),
        hazard_coral_threshold=Float64(get(params, :hazard_coral_threshold, 0.15)),
        hazard_coral_width=Float64(get(params, :hazard_coral_width, 0.03)),
        hazard_max=Float64(get(params, :hazard_max, 0.0)),
        preferred_prey_starvation=Bool(get(params, :preferred_prey_starvation, false)),
        starvation_preferred_threshold=Float64(get(params, :starvation_preferred_threshold, 0.05)),
        starvation_memory_years=Float64(get(params, :starvation_memory_years, 0.0)),
        adult_senescence=Bool(get(params, :adult_senescence, false)),
        senescence_age=Int(get(params, :senescence_age, 6)),
        senescence_mortality=Float64(get(params, :senescence_mortality, 0.8)),
        per_capita_consumption=Bool(get(params, :per_capita_consumption, false)),
        consumption_m2_per_adult_year=Float64(get(params, :consumption_m2_per_adult_year, 10.0)),
        prey_switching=Bool(get(params, :prey_switching, false)),
        switch_midpoint=Float64(get(params, :switch_midpoint, 0.10)),
        switch_steepness=Float64(get(params, :switch_steepness, 70.0)),
        preferred_condition=Bool(get(params, :preferred_condition, false)),
        condition_saturation=Float64(get(params, :condition_saturation, 0.10)),
        condition_memory_years=Float64(get(params, :condition_memory_years, 1.0)),
        diet_linked_starvation=Bool(get(params, :diet_linked_starvation, false)),
        nonpreferred_feeding_efficiency=Float64(get(params, :nonpreferred_feeding_efficiency, 1.0)),
        juvenile_predation=Bool(get(params, :juvenile_predation, false)),
        juvenile_predation_capacity=Float64(get(params, :juvenile_predation_capacity, 20.0)),
        preferred_maturation=Bool(get(params, :preferred_maturation, false)),
        allee_stretched=Bool(get(params, :allee_stretched, false)),
        allee_exponent=Float64(get(params, :allee_exponent, 0.5)),
        starvation_post_grazing=Bool(get(params, :starvation_post_grazing, false)),
        starvation_grace_years=Int(get(params, :starvation_grace_years, 0)),
        outbreak_collapse=Bool(get(params, :outbreak_collapse, false)),
        collapse_threshold_ha=Float64(get(params, :collapse_threshold_ha, 67.5)),
        collapse_duration_years=Int(get(params, :collapse_duration_years, 4)),
        collapse_background_ha=Float64(get(params, :collapse_background_ha, 0.25)),
        collapse_juveniles=Bool(get(params, :collapse_juveniles, false))
    )
end

cots_timestep!(model::CotsHuman, F::Float64, S::Float64) = COTSMod.cots_timestep!(model, F, S)

function init_cots_populations(n_locs::Int, params; outbreak_fraction::Float64=0.25, seed_locs::Union{Set{Int},Nothing}=nothing, init_density::Float64=0.1, rng=Random.GLOBAL_RNG)::CotsState
    return COTSMod.init_cots_populations(n_locs, _cots_runtime_params(params); outbreak_fraction=outbreak_fraction, seed_locs=seed_locs, init_density=init_density, rng=rng)
end

function init_cots_from_spatial(n_locs::Int, params, probabilities::AbstractVector{<:Real}; initial_multiplier::Float64=parse(Float64, get(ENV, "COTS_INITIAL_MULTIPLIER", "1.5")))::CotsState
    return COTSMod.init_cots_from_spatial(n_locs, _cots_runtime_params(params), probabilities; initial_multiplier=initial_multiplier)
end

function cots_mortality!(C_cover_t::AbstractArray{Float64,3}, cots_models::CotsState, prey_map::CotsPreyMap)::Nothing
    return COTSMod.cots_mortality!(C_cover_t, cots_models, prey_map)
end

function disperse_cots_larvae!(cots_models::CotsState, conn::SparseMatrixCSC{Float64,Int64}; immigration_scalar::Float64=1.0)::Nothing
    return COTSMod.disperse_cots_larvae!(cots_models, conn; immigration_scalar=immigration_scalar)
end

function disperse_cots_larvae_by_reef!(cots_models::CotsState, reef_conn::AbstractMatrix, site_to_reef::AbstractVector{<:Integer}; immigration_scalar::Float64=1.0, source_survival=nothing, reef_areas=nothing)::Nothing
    return COTSMod.disperse_cots_larvae_by_reef!(cots_models, reef_conn, site_to_reef; immigration_scalar=immigration_scalar, source_survival=source_survival, reef_areas=reef_areas)
end

function inject_upstream_pulse!(cots_models::CotsState, pulse_locs::Set{Int}; pulse_val::Float64=2.0)::Nothing
    return COTSMod.inject_upstream_pulse!(cots_models, pulse_locs; pulse_val=pulse_val)
end

function inject_upstream_pulse!(cots_models::CotsState, pulse_locs::Set{Int}, pulse_vals::AbstractVector{Float64})::Nothing
    return COTSMod.inject_upstream_pulse!(cots_models, pulse_locs, pulse_vals)
end

function initialize_cots(n_locs::Int, params; enabled::Bool=true, spatial_initial_density::Union{AbstractVector{<:Real},Nothing}=nothing, seed_locs::Union{Set{Int},Nothing}=nothing, init_density::Float64=0.1, initial_multiplier::Float64=parse(Float64, get(ENV, "COTS_INITIAL_MULTIPLIER", "1.5")), rng=Random.GLOBAL_RNG)::CotsState
    return COTSMod.initialize_cots(n_locs, _cots_runtime_params(params); enabled=enabled, spatial_initial_density=spatial_initial_density, seed_locs=seed_locs, init_density=init_density, initial_multiplier=initial_multiplier, rng=rng)
end

"""Read an opt-in, site-keyed COTS initial state in COTS ha^-1 by age class."""
function load_cots_initial_state_csv(path::AbstractString, loc_ids::AbstractVector{<:AbstractString})
    isfile(path) || error("Missing COTS initial-state CSV: $path")
    table = CSV.read(path, DataFrame)
    columns = (:site_id, :juvenile_ha, :subadult_ha, :adult_ha)
    all(in(propertynames(table)), columns) || error("COTS initial-state CSV lacks required columns")
    nrow(table) == length(loc_ids) || error("COTS initial-state CSV has the wrong site count")
    ids = String.(table.site_id)
    length(unique(ids)) == length(ids) || error("COTS initial-state CSV has duplicate site IDs")
    lookup = Dict(id => row for (row, id) in enumerate(ids))
    rows = [get(lookup, String(id), 0) for id in loc_ids]
    any(iszero, rows) && error("COTS initial-state CSV does not cover every domain site")
    state = Matrix{Float64}(table[rows, collect(columns[2:end])])
    all(isfinite, state) && all(>=(0.0), state) ||
        error("COTS initial-state densities must be finite and non-negative")
    return state
end

"""
    load_cots_forced_adults(path, loc_ids, calendar_years)

Read an opt-in site x calendar-year table of prescribed adult COTS density
(COTS ha^-1 of site area; columns `site_id` and one per year) into a
`timesteps x locations` matrix for forced-grazing coral calibration.
"""
function load_cots_forced_adults(path::AbstractString, loc_ids::AbstractVector{<:AbstractString},
    calendar_years::AbstractVector{<:Integer})::Matrix{Float64}
    isfile(path) || error("Missing forced COTS adults CSV: $path")
    table = CSV.read(path, DataFrame; types=Dict(:site_id => String))
    ids = String.(table.site_id)
    length(unique(ids)) == length(ids) || error("Forced COTS CSV has duplicate site IDs")
    lookup = Dict(id => row for (row, id) in enumerate(ids))
    rows = [get(lookup, String(id), 0) for id in loc_ids]
    any(iszero, rows) && error("Forced COTS CSV does not cover every domain site")
    # Collect first: calendar years may be a DimensionalData vector, and a
    # comprehension over it would keep the dimension wrapper.
    year_columns = String[string(Int(year)) for year in collect(calendar_years)]
    all(in(names(table)), year_columns) || error("Forced COTS CSV does not cover every timestep")
    adults = permutedims(Matrix{Float64}(table[rows, year_columns]))
    all(isfinite, adults) && all(>=(0.0), adults) ||
        error("Forced COTS adult densities must be finite and non-negative")
    return adults
end

"""
    parse_cots_supply_window(spec) -> Union{Nothing,Tuple{Int,Int,Float64}}

Labelled counterfactual (`COTS_EXTERNAL_SUPPLY_WINDOW = "<start>:<end>:<factor>"`;
crash_mechanism_protocol.md, 2026-10-14): multiply the evidence-derived external outbreak
supply by `factor` in calendar years `start..end`. Empty = established path.
"""
function parse_cots_supply_window(spec::AbstractString)::Union{Nothing,Tuple{Int,Int,Float64}}
    isempty(spec) && return nothing
    parts = split(spec, ":")
    length(parts) == 3 || error("COTS_EXTERNAL_SUPPLY_WINDOW must be <start>:<end>:<factor>")
    start_year, end_year, factor = parse(Int, parts[1]), parse(Int, parts[2]), parse(Float64, parts[3])
    start_year <= end_year || error("COTS_EXTERNAL_SUPPLY_WINDOW years are reversed")
    isfinite(factor) && factor >= 0.0 || error("COTS_EXTERNAL_SUPPLY_WINDOW factor must be finite and non-negative")
    return (start_year, end_year, factor)
end

"""Supply multiplier for `calendar_year` under a parsed window (1 outside it or when unset)."""
cots_supply_window_factor(window::Union{Nothing,Tuple{Int,Int,Float64}}, calendar_year::Integer)::Float64 =
    isnothing(window) || !(window[1] <= calendar_year <= window[2]) ? 1.0 : window[3]

"""
    load_cots_larval_survival(path, loc_ids, calendar_years) -> Matrix{Float64}

Opt-in larval-survival index (`COTS_LARVAL_SURVIVAL_CSV`; crash_mechanism_protocol.md,
2026-10-15): columns `site_id, year, multiplier` (dimensionless, finite, >= 0). Returns
`[timestep, location]` multipliers; site-years absent from the table are 1. Unknown sites or
invalid multipliers are errors.
"""
function load_cots_larval_survival(path::AbstractString, loc_ids::AbstractVector{<:AbstractString},
    calendar_years::AbstractVector{<:Integer})::Matrix{Float64}
    isfile(path) || error("Missing COTS larval-survival CSV: $path")
    table = CSV.read(path, DataFrame; types=Dict(:site_id => String))
    all(in(names(table)), ["site_id", "year", "multiplier"]) ||
        error("COTS larval-survival CSV needs columns site_id, year, multiplier")
    sites = Dict(String(id) => loc for (loc, id) in enumerate(loc_ids))
    years = Dict(Int(year) => t for (t, year) in enumerate(collect(calendar_years)))
    multipliers = ones(Float64, length(years), length(loc_ids))
    for row in eachrow(table)
        loc = get(sites, row.site_id, 0)
        loc == 0 && error("COTS larval-survival site $(row.site_id) is not a domain site")
        value = Float64(row.multiplier)
        isfinite(value) && value >= 0.0 || error("COTS larval-survival multipliers must be finite and non-negative")
        t = get(years, Int(row.year), 0)
        t == 0 && continue
        multipliers[t, loc] = value
    end
    return multipliers
end

"""
    load_cots_external_survival(path, calendar_years) -> Vector{Float64}

Opt-in year-specific survival multiplier on external pelagic supply
(`COTS_EXTERNAL_SURVIVAL_CSV`; crash_mechanism_protocol.md, 2026-10-19, twenty-fourth round):
columns `year, multiplier` (dimensionless, finite, >= 0), one row per year. Returns one value per
timestep; years absent from the table are 1. Duplicate years or invalid multipliers are errors.
"""
function load_cots_external_survival(path::AbstractString,
    calendar_years::AbstractVector{<:Integer})::Vector{Float64}
    isfile(path) || error("Missing COTS external-survival CSV: $path")
    table = CSV.read(path, DataFrame)
    all(in(names(table)), ["year", "multiplier"]) ||
        error("COTS external-survival CSV needs columns year, multiplier")
    allunique(table.year) || error("Duplicate years in COTS external-survival CSV")
    values = Dict(Int(r.year) => Float64(r.multiplier) for r in eachrow(table))
    all(v -> isfinite(v) && v >= 0.0, Base.values(values)) ||
        error("COTS external-survival multipliers must be finite and non-negative")
    return [get(values, Int(y), 1.0) for y in calendar_years]
end

"""Scale each location's pre-pelagic larval production (separated recruitment) by the
timestep's survival multiplier, at the source, before dispersal."""
function apply_cots_larval_survival!(cots_models::CotsState, multipliers::AbstractVector{<:Real})::Nothing
    length(multipliers) == length(cots_models) ||
        throw(DimensionMismatch("larval-survival multipliers must match locations"))
    for (model, m) in zip(cots_models, multipliers)
        model.last_local_fecundity *= m
    end
    return nothing
end

"""
    load_cots_settlement_modifier(path, loc_ids) -> Vector{Float64}

Opt-in reef-habitat settlement modifier (`COTS_SETTLEMENT_MODIFIER_CSV`;
crash_mechanism_protocol.md, 2026-10-19, wave exposure): columns `site_id, multiplier`
(dimensionless, finite, in [0, 1]), time-invariant. Every domain site must appear exactly
once; unknown, duplicate or missing sites and invalid multipliers are errors.
"""
function load_cots_settlement_modifier(path::AbstractString,
    loc_ids::AbstractVector{<:AbstractString})::Vector{Float64}
    isfile(path) || error("Missing COTS settlement-modifier CSV: $path")
    table = CSV.read(path, DataFrame; types=Dict(:site_id => String))
    all(in(names(table)), ["site_id", "multiplier"]) ||
        error("COTS settlement-modifier CSV needs columns site_id, multiplier")
    sites = Dict(String(id) => loc for (loc, id) in enumerate(loc_ids))
    multipliers = fill(NaN, length(loc_ids))
    for row in eachrow(table)
        loc = get(sites, row.site_id, 0)
        loc == 0 && error("COTS settlement-modifier site $(row.site_id) is not a domain site")
        isnan(multipliers[loc]) || error("Duplicate COTS settlement-modifier site $(row.site_id)")
        value = Float64(row.multiplier)
        isfinite(value) && 0.0 <= value <= 1.0 ||
            error("COTS settlement modifiers must be finite and in [0, 1]")
        multipliers[loc] = value
    end
    any(isnan, multipliers) && error("COTS settlement-modifier CSV is missing domain sites")
    return multipliers
end

"""Scale each location's settlement (separated recruitment) by its habitat modifier, after the
COTS transition and before dispersal. The settlement gate scales local retention, internal
immigration, external supply and pulses exactly once; the background immigration already
settled in the transition (`N[1]`, separated recruitment) is scaled here as well."""
function apply_cots_settlement_modifier!(cots_models::CotsState, multipliers::AbstractVector{<:Real})::Nothing
    length(multipliers) == length(cots_models) ||
        throw(DimensionMismatch("settlement modifiers must match locations"))
    for (model, m) in zip(cots_models, multipliers)
        model.params.separated_recruitment ||
            error("COTS settlement modifier requires separated recruitment")
        model.last_settlement_gate *= m
        model.N[1] *= m
        model.last_background_immigration *= m
    end
    return nothing
end

"""
    load_cots_cull_removals(path, loc_ids, calendar_years) -> Array{Float64,3}

Opt-in observed culling (`COTS_CULL_REMOVALS_CSV`; crash_mechanism_protocol.md,
2026-10-14). Columns `site_id, report_year, phase, juveniles_ha, adults_ha` (COTS ha^-1
of site area). Returns requests indexed `[timestep, k, location]` with k = 1 pre-spawning
juveniles, 2 pre-spawning adults, 3 post-spawning juveniles, 4 post-spawning adults; a
row applies at the timestep whose calendar year equals `report_year`. Rows outside the
simulated years are ignored; unknown sites, phases or invalid densities are errors.
"""
function load_cots_cull_removals(path::AbstractString, loc_ids::AbstractVector{<:AbstractString},
    calendar_years::AbstractVector{<:Integer})::Array{Float64,3}
    isfile(path) || error("Missing COTS cull removals CSV: $path")
    table = CSV.read(path, DataFrame; types=Dict(:site_id => String, :phase => String))
    required = ["site_id", "report_year", "phase", "juveniles_ha", "adults_ha"]
    all(in(names(table)), required) || error("COTS cull CSV needs columns $(required)")
    sites = Dict(String(id) => loc for (loc, id) in enumerate(loc_ids))
    years = Dict(Int(year) => t for (t, year) in enumerate(collect(calendar_years)))
    requests = zeros(Float64, length(years), 4, length(loc_ids))
    for row in eachrow(table)
        loc = get(sites, row.site_id, 0)
        loc == 0 && error("COTS cull CSV site $(row.site_id) is not a domain site")
        row.phase in ("pre", "post") || error("COTS cull phase must be pre or post")
        juveniles, adults = Float64(row.juveniles_ha), Float64(row.adults_ha)
        isfinite(juveniles) && isfinite(adults) && juveniles >= 0.0 && adults >= 0.0 ||
            error("COTS cull densities must be finite and non-negative")
        t = get(years, Int(row.report_year), 0)
        t == 0 && continue
        offset = row.phase == "pre" ? 0 : 2
        requests[t, offset + 1, loc] += juveniles
        requests[t, offset + 2, loc] += adults
    end
    return requests
end

"""
    reef_proportional_requests(cots_models, requests, stage, groups, areas) -> Vector{Float64}

Opt-in reef-level cull allocation (`COTS_CULL_ALLOCATION=reef_proportional`;
crash_mechanism_protocol.md, 2026-10-19). Each reef's requested count
`R = sum(requests .* areas)` (COTS) is reassigned to its sites in proportion to the COTS
present in `stage` (2 juveniles, 3 adults) times site area, as divers search the reef and
cull where COTS are. Site requests become `R * N_s / sum(N .* areas)` (COTS ha^-1), so the
reef total is fully removable whenever the reef holds at least `R`; otherwise every site is
cleared. A reef with no COTS keeps zero requests (the total is unmet).
"""
function reef_proportional_requests(cots_models::CotsState, requests::AbstractVector{<:Real},
    stage::Int, groups::AbstractVector{<:AbstractVector{Int}}, areas::AbstractVector{<:Real})::Vector{Float64}
    stage in (2, 3) || throw(ArgumentError("stage must be 2 (juveniles) or 3 (adults)"))
    length(requests) == length(areas) == length(cots_models) ||
        throw(DimensionMismatch("cull requests, areas and locations must match"))
    out = zeros(Float64, length(requests))
    for sites in groups
        total = sum(Float64(requests[s]) * areas[s] for s in sites)
        total > 0.0 || continue
        stock = sum(cots_models[s].N[stage] * areas[s] for s in sites)
        stock > 0.0 || continue
        for s in sites
            out[s] = total * cots_models[s].N[stage] / stock
        end
    end
    return out
end

"""Apply one phase (`:pre` or `:post`) of the timestep's cull requests and log what was
removed into `removed_log[tstep, k, :]` (same k layout as the requests). With `groups`
and `areas`, requests are first reallocated within each reef (`reef_proportional_requests`)."""
function apply_cots_cull!(cots_models::CotsState, requests::Array{Float64,3},
    removed_log::Array{Float64,3}, tstep::Int, phase::Symbol;
    groups=nothing, areas=nothing)::Nothing
    phase in (:pre, :post) || throw(ArgumentError("phase must be :pre or :post"))
    k = phase == :pre ? 1 : 3
    juvenile_requests, adult_requests = if isnothing(groups)
        @view(requests[tstep, k, :]), @view(requests[tstep, k + 1, :])
    else
        reef_proportional_requests(cots_models, @view(requests[tstep, k, :]), 2, groups, areas),
        reef_proportional_requests(cots_models, @view(requests[tstep, k + 1, :]), 3, groups, areas)
    end
    juveniles, adults = COTSMod.apply_removals!(cots_models, juvenile_requests, adult_requests)
    removed_log[tstep, k, :] .= juveniles
    removed_log[tstep, k + 1, :] .= adults
    return nothing
end

"""
    load_cots_cull_effort(path, group_ids) -> Vector{Float64}

Opt-in forward effort culling (`COTS_CULL_EFFORT_CSV`; crash_mechanism_protocol.md, 2026-10-21,
M1). Columns `reef_id, hours_per_year` (diver-hours yr^-1, summed dive bottom time). Returns
hours aligned with `group_ids` (the reef id of each cull group); reefs absent from the table get 0.
Unknown reefs or invalid hours are errors.
"""
function load_cots_cull_effort(path::AbstractString, group_ids::AbstractVector{<:AbstractString})::Vector{Float64}
    isfile(path) || error("Missing COTS cull effort CSV: $path")
    table = CSV.read(path, DataFrame; types=Dict(:reef_id => String))
    all(in(names(table)), ["reef_id", "hours_per_year"]) ||
        error("COTS cull effort CSV needs columns reef_id, hours_per_year")
    lookup = Dict(String(id) => g for (g, id) in enumerate(group_ids))
    hours = zeros(Float64, length(group_ids))
    for row in eachrow(table)
        g = get(lookup, row.reef_id, 0)
        g == 0 && error("COTS cull effort reef $(row.reef_id) is not a domain reef")
        h = Float64(row.hours_per_year)
        isfinite(h) && h >= 0.0 || error("COTS cull effort hours must be finite and non-negative")
        hours[g] += h
    end
    return hours
end

"""Catch per diver-hour at adult density `D` (COTS ha^-1): the SALAD-cull anchor
`D = c * CPUE^k` inverted, capped at `cpue_max` (COTS h^-1)."""
cots_cull_cpue(D::Real; c::Real=7.874, k::Real=0.549, cpue_max::Real=20.0)::Float64 =
    D > 0.0 ? min((Float64(D) / c)^(1.0 / k), Float64(cpue_max)) : 0.0

"""Reef stock (COTS) and mean adult density (COTS ha^-1 of site area) per cull group.
`areas` must be site areas in ha: stock = sum(N_adult [COTS ha^-1] x area [ha])."""
function cots_group_stock(cots_models::CotsState, groups::AbstractVector{<:AbstractVector{Int}},
    areas::AbstractVector{<:Real})::Tuple{Vector{Float64},Vector{Float64}}
    stock = [sum(cots_models[s].N[3] * areas[s] for s in sites; init=0.0) for sites in groups]
    area = [sum(Float64(areas[s]) for s in sites; init=0.0) for sites in groups]
    return stock, [a > 0.0 ? x / a : 0.0 for (x, a) in zip(stock, area)]
end

"""
    targeted_cull_hours(cots_models, total_hours, groups, areas; threshold_ha) -> Vector{Float64}

`targeted` mode: `total_hours` reallocated to groups whose adult density is at least
`threshold_ha`, in proportion to their stock. No group qualifies -> all zeros.
"""
function targeted_cull_hours(cots_models::CotsState, total_hours::Real,
    groups::AbstractVector{<:AbstractVector{Int}}, areas::AbstractVector{<:Real};
    threshold_ha::Real=3.0)::Vector{Float64}
    stock, density = cots_group_stock(cots_models, groups, areas)
    eligible = density .>= threshold_ha
    total = sum(stock[eligible]; init=0.0)
    total > 0.0 || return zeros(Float64, length(groups))
    return [eligible[g] ? Float64(total_hours) * stock[g] / total : 0.0 for g in eachindex(groups)]
end

"""
    effort_cull_requests(cots_models, hours, groups, areas; c, k, cpue_max) -> Vector{Float64}

Adult removal requests (COTS ha^-1 per site) for one phase: each group removes
`min(stock, hours * cpue(D))` COTS, allocated to its sites in proportion to adults x area.
`areas` are site areas in ha (removals are absolute counts of COTS).
"""
function effort_cull_requests(cots_models::CotsState, hours::AbstractVector{<:Real},
    groups::AbstractVector{<:AbstractVector{Int}}, areas::AbstractVector{<:Real};
    c::Real=7.874, k::Real=0.549, cpue_max::Real=20.0)::Vector{Float64}
    length(hours) == length(groups) || throw(DimensionMismatch("effort hours and groups must match"))
    stock, density = cots_group_stock(cots_models, groups, areas)
    out = zeros(Float64, length(cots_models))
    for (g, sites) in enumerate(groups)
        (hours[g] > 0.0 && stock[g] > 0.0) || continue
        removed = min(stock[g], Float64(hours[g]) * cots_cull_cpue(density[g]; c=c, k=k, cpue_max=cpue_max))
        for s in sites
            out[s] = removed * cots_models[s].N[3] / stock[g]
        end
    end
    return out
end

"""Apply one phase of forward effort culling (adults only) and add the removals to
`removed_log[tstep, k + 1, :]` (k = 1 pre, 3 post), after any observed culls."""
function apply_cots_effort_cull!(cots_models::CotsState, hours::AbstractVector{<:Real},
    groups, areas, removed_log::Array{Float64,3}, tstep::Int, phase::Symbol; kwargs...)::Nothing
    phase in (:pre, :post) || throw(ArgumentError("phase must be :pre or :post"))
    requests = effort_cull_requests(cots_models, hours, groups, areas; kwargs...)
    _, adults = COTSMod.apply_removals!(cots_models, zeros(Float64, length(requests)), requests)
    removed_log[tstep, phase == :pre ? 2 : 4, :] .+= adults
    return nothing
end

"""
    load_cots_state_reset(path, loc_ids) -> (year, groups, targets, labels)

Opt-in projection-start reset (`COTS_STATE_RESET_CSV`; protocol 2026-10-21, M2). Columns
`site_id, group, year, target_adults_ha`; one year for the whole table; each group has one target.
"""
function load_cots_state_reset(path::AbstractString, loc_ids::AbstractVector{<:AbstractString})
    isfile(path) || error("Missing COTS state reset CSV: $path")
    table = CSV.read(path, DataFrame; types=Dict(:site_id => String, :group => String))
    all(in(names(table)), ["site_id", "group", "year", "target_adults_ha"]) ||
        error("COTS state reset CSV needs columns site_id, group, year, target_adults_ha")
    years = unique(Int.(table.year))
    length(years) == 1 || error("COTS state reset CSV must hold a single year")
    sites = Dict(String(id) => loc for (loc, id) in enumerate(loc_ids))
    labels = unique(String.(table.group))
    groups, targets = Vector{Int}[], Float64[]
    for label in labels
        rows = table[table.group .== label, :]
        target = unique(Float64.(rows.target_adults_ha))
        length(target) == 1 && isfinite(target[1]) && target[1] >= 0.0 ||
            error("COTS state reset group $label needs one finite non-negative target")
        locs = [get(sites, id, 0) for id in rows.site_id]
        any(iszero, locs) && error("COTS state reset group $label has a site outside the domain")
        push!(groups, locs)
        push!(targets, target[1])
    end
    length(unique(vcat(groups...))) == sum(length, groups) || error("COTS state reset groups overlap")
    return only(years), groups, targets, labels
end

"""Scale one location's COTS state: adults (with age profile and large adults) by `f_adult`,
recruits and juveniles (N[1], N[2]) by `f_juvenile`."""
function scale_cots_state!(model, f_adult::Real, f_juvenile::Real)::Nothing
    model.N[3] *= f_adult
    isempty(model.adult_ages) || (model.adult_ages .*= f_adult)
    model.large_adults *= f_adult
    model.N[1] *= f_juvenile
    model.N[2] *= f_juvenile
    return nothing
end

"""
    apply_cots_state_reset!(cots_models, groups, targets, areas; juveniles=true, others=:median)
        -> (factors, other_factor)

Rescale each group's area-weighted adult density to its target (factor `f = target / current`;
`f = 1` where the group holds no adults). Juveniles take the same factor when `juveniles`.
Locations outside every group take the median group factor (`others = :median`) or 1 (`:none`).
"""
function apply_cots_state_reset!(cots_models::CotsState, groups::AbstractVector{<:AbstractVector{Int}},
    targets::AbstractVector{<:Real}, areas::AbstractVector{<:Real};
    juveniles::Bool=true, others::Symbol=:median)
    others in (:median, :none) || throw(ArgumentError("others must be :median or :none"))
    _, density = cots_group_stock(cots_models, groups, areas)
    factors = [d > 0.0 ? Float64(t) / d : 1.0 for (t, d) in zip(targets, density)]
    for (sites, f) in zip(groups, factors), s in sites
        scale_cots_state!(cots_models[s], f, juveniles ? f : 1.0)
    end
    other = others == :median && !isempty(factors) ? median(factors) : 1.0
    grouped = Set(vcat(groups...))
    for loc in eachindex(cots_models)
        loc in grouped || scale_cots_state!(cots_models[loc], other, juveniles ? other : 1.0)
    end
    return factors, other
end

"""Coral-coupled outside supply factor (protocol 2026-10-21, M3): `min(1, current / reference)`;
1 when the reference is not positive."""
cots_external_coral_factor(current::Real, reference::Real)::Float64 =
    reference > 0.0 ? clamp(Float64(current) / Float64(reference), 0.0, 1.0) : 1.0

function apply_cots_initial_state!(models::CotsState, state::AbstractMatrix{<:Real})
    size(state) == (length(models), 3) || error("COTS initial-state dimensions differ")
    all(isfinite, state) && all(>=(0.0), state) ||
        error("COTS initial-state densities must be finite and non-negative")
    for (model, row) in zip(models, eachrow(state))
        model.N .= row
        model.large_adults = model.params.size_weighted_fecundity ?
            row[3] * model.params.size_initial_large_fraction : 0.0
        model.last_effective_breeders = 0.0
        model.last_small_adult_growth = 0.0
        model.burden_ha = 0.0
        model.last_burden_ha = 0.0
        model.last_hazard = 0.0
        model.last_hazard_deaths_ha = 0.0
        # Spread adults over the stable age profile (when enabled) and clear food memory.
        COTSMod.reset_adult_ages!(model)
    end
    return nothing
end

function apply_predation!(coral_cover::AbstractArray{Float64,3}, cots_state::CotsState,
    prey_map::CotsPreyMap; habitable_fraction=nothing)::Nothing
    return COTSMod.apply_predation!(coral_cover, cots_state, prey_map;
        habitable_fraction=habitable_fraction)
end

"""Experimental grazing bypass: advance COTS normally but leave coral untouched in a labelled year window."""
function apply_cots_predation_diagnostic!(coral_cover::AbstractArray{Float64,3},
    cots_state::CotsState, prey_map::CotsPreyMap, calendar_year::Int,
    grazing_off_years::Union{Nothing,AbstractUnitRange{Int}};
    habitable_fraction=nothing)::Nothing
    if !isnothing(grazing_off_years) && calendar_year in grazing_off_years
        apply_predation!(copy(coral_cover), cots_state, prey_map;
            habitable_fraction=habitable_fraction)
    else
        apply_predation!(coral_cover, cots_state, prey_map;
            habitable_fraction=habitable_fraction)
    end
    return nothing
end

function disperse_larvae!(cots_state::CotsState, conn::SparseMatrixCSC{Float64,Int64}; scalar::Float64=1.0)::Nothing
    return COTSMod.disperse_larvae!(cots_state, conn; scalar=scalar)
end

function disperse_larvae_by_reef!(cots_state::CotsState, reef_conn::AbstractMatrix, site_to_reef::AbstractVector{<:Integer}; scalar::Float64=1.0, source_survival=nothing, reef_areas=nothing)::Nothing
    return COTSMod.disperse_larvae_by_reef!(cots_state, reef_conn, site_to_reef; scalar=scalar, source_survival=source_survival, reef_areas=reef_areas)
end

function apply_external_supply!(cots_state::CotsState, pulse_locs::Set{Int}; pulse_val::Float64=2.0)::Nothing
    return COTSMod.apply_external_supply!(cots_state, pulse_locs; pulse_val=pulse_val)
end

function apply_external_supply!(cots_state::CotsState, pulse_locs::Set{Int}, pulse_vals::AbstractVector{Float64})::Nothing
    return COTSMod.apply_external_supply!(cots_state, pulse_locs, pulse_vals)
end

function apply_external_supply!(cots_state::CotsState, pulse_vals::AbstractVector{<:Real})::Nothing
    return COTSMod.apply_external_supply!(cots_state, pulse_vals)
end

function apply_external_larval_supply!(
    cots_state::CotsState,
    potential_supply::AbstractVector{<:Real},
    pelagic_supply::AbstractVector{<:Real}
)::Nothing
    return COTSMod.apply_external_larval_supply!(
        cots_state, potential_supply, pelagic_supply
    )
end

cots_flow_diagnostics(cots_state::CotsState) = COTSMod.cots_flow_diagnostics(cots_state)
cots_hazard_diagnostics(cots_state::CotsState) = COTSMod.cots_hazard_diagnostics(cots_state)
cots_mechanism_diagnostics(cots_state::CotsState) = COTSMod.cots_mechanism_diagnostics(cots_state)
