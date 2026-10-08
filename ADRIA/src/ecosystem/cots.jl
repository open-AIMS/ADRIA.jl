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
        low_cover_strength=Float64(get(params, :low_cover_strength, 0.0))
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

function apply_cots_initial_state!(models::CotsState, state::AbstractMatrix{<:Real})
    size(state) == (length(models), 3) || error("COTS initial-state dimensions differ")
    all(isfinite, state) && all(>=(0.0), state) ||
        error("COTS initial-state densities must be finite and non-negative")
    for (model, row) in zip(models, eachrow(state))
        model.N .= row
    end
    return nothing
end

function apply_predation!(coral_cover::AbstractArray{Float64,3}, cots_state::CotsState, prey_map::CotsPreyMap)::Nothing
    return COTSMod.apply_predation!(coral_cover, cots_state, prey_map)
end

"""Experimental grazing bypass: advance COTS normally but leave coral untouched in a labelled year window."""
function apply_cots_predation_diagnostic!(coral_cover::AbstractArray{Float64,3},
    cots_state::CotsState, prey_map::CotsPreyMap, calendar_year::Int,
    grazing_off_years::Union{Nothing,AbstractUnitRange{Int}})::Nothing
    if !isnothing(grazing_off_years) && calendar_year in grazing_off_years
        apply_predation!(copy(coral_cover), cots_state, prey_map)
    else
        apply_predation!(coral_cover, cots_state, prey_map)
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
