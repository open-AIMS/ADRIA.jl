# ADRIA adapter for the crown-of-thorns starfish (COTS) submodel in COTSMod.jl.
#
# ADRIA owns the parameter table, domain wiring, coral-cover tensors, result logging and
# scenario timing. COTSMod owns the COTS state transition, predation, larval dispersal and
# initialisation. Units follow COTSMod: COTS densities are COTS ha^-1.

const AbstractCotsModel = COTSMod.AbstractCOTSModel
const CotsPreyMap = COTSMod.COTSPreyMap
const CotsPreyState = COTSMod.COTSPreyState
const CotsHuman = COTSMod.COTSHuman
const CotsState = COTSMod.COTSState

"""
    cots_runtime_params(param_set::YAXArray)::COTSMod.COTSParams

Build COTSMod's runtime parameters from a scenario's `CotsParams` factors, following
`COTS_PARAM_MAP`. Mechanisms that are not exposed as factors keep COTSMod's defaults (off).
"""
function cots_runtime_params(param_set::YAXArray)::COTSMod.COTSParams
    kwargs = (
        field => _cots_param_value(factor, param_set[At(String(factor))])
        for (factor, field) in pairs(COTS_PARAM_MAP)
    )
    return COTSMod.COTSParams(; kwargs...)
end

# Convert a factor value to the type COTSMod expects, from the factor's `ptype`
function _cots_param_value(factor::Symbol, val)
    ptype = getproperty(COTS_FACTOR_SPECS, factor).ptype
    ptype == "unordered categorical" && return val == 1
    ptype == "ordered discrete" && return round(Int, val)
    return Float64(val)
end

cots_timestep!(model::CotsHuman, F::Float64, S::Float64) = COTSMod.cots_timestep!(
    model, F, S
)

function init_cots_populations(
    n_locs::Int, params::COTSMod.COTSParams; outbreak_fraction::Float64=0.25,
    seed_locs::Union{Set{Int},Nothing}=nothing, init_density::Float64=0.1,
    rng=Random.GLOBAL_RNG
)::CotsState
    return COTSMod.init_cots_populations(
        n_locs, params; outbreak_fraction=outbreak_fraction, seed_locs=seed_locs,
        init_density=init_density, rng=rng
    )
end

function initialize_cots(
    n_locs::Int, params::COTSMod.COTSParams; enabled::Bool=true,
    spatial_initial_density::Union{AbstractVector{<:Real},Nothing}=nothing,
    seed_locs::Union{Set{Int},Nothing}=nothing, init_density::Float64=0.1,
    initial_multiplier::Float64=1.5, rng=Random.GLOBAL_RNG
)::CotsState
    return COTSMod.initialize_cots(
        n_locs, params; enabled=enabled, spatial_initial_density=spatial_initial_density,
        seed_locs=seed_locs, init_density=init_density,
        initial_multiplier=initial_multiplier, rng=rng
    )
end

function cots_mortality!(
    C_cover_t::AbstractArray{Float64,3}, cots_models::CotsState, prey_map::CotsPreyMap
)::Nothing
    return COTSMod.cots_mortality!(C_cover_t, cots_models, prey_map)
end

function disperse_cots_larvae!(
    cots_models::CotsState, conn::SparseMatrixCSC{Float64,Int64};
    immigration_scalar::Float64=1.0
)::Nothing
    return COTSMod.disperse_cots_larvae!(
        cots_models, conn; immigration_scalar=immigration_scalar
    )
end

"""
    apply_predation!(coral_cover, cots_state, prey_map; habitable_fraction=nothing)

Advance each location's COTS population one year and remove the coral it eats from
`coral_cover` (groups x sizes x locations, relative cover). `habitable_fraction` converts
per-capita grazing (m^2 per adult) into habitable-cover fractions when that mechanism is on.
"""
function apply_predation!(
    coral_cover::AbstractArray{Float64,3}, cots_state::CotsState, prey_map::CotsPreyMap;
    habitable_fraction=nothing
)::Nothing
    return COTSMod.apply_predation!(
        coral_cover, cots_state, prey_map; habitable_fraction=habitable_fraction
    )
end

"""
    disperse_larvae!(cots_state, conn; scalar=1.0)

Disperse COTS larvae between locations through the source-to-sink connectivity `conn`.
"""
function disperse_larvae!(
    cots_state::CotsState, conn::SparseMatrixCSC{Float64,Int64}; scalar::Float64=1.0
)::Nothing
    return COTSMod.disperse_larvae!(cots_state, conn; scalar=scalar)
end

"""
    disperse_larvae_by_reef!(cots_state, reef_conn, site_to_reef; scalar=1.0, reef_areas=nothing)

Disperse COTS larvae between reefs through the reef x reef connectivity `reef_conn`, then
share each reef's settlers among its locations. `reef_areas` (habitable ha per reef) is
required by separated recruitment.
"""
function disperse_larvae_by_reef!(
    cots_state::CotsState, reef_conn::AbstractMatrix,
    site_to_reef::AbstractVector{<:Integer}; scalar::Float64=1.0, reef_areas=nothing
)::Nothing
    return COTSMod.disperse_larvae_by_reef!(
        cots_state, reef_conn, site_to_reef; scalar=scalar, reef_areas=reef_areas
    )
end

"""
    apply_external_larval_supply!(cots_state, supply::AbstractVector{<:Real})

Add larvae arriving from outside the domain to each location's settlers.
"""
function apply_external_larval_supply!(
    cots_state::CotsState, supply::AbstractVector{<:Real}
)::Nothing
    return COTSMod.apply_external_larval_supply!(cots_state, supply, supply)
end

"""
    set_initial_state!(cots_state, state::AbstractMatrix{<:Real})

Set each location's COTS densities (COTS ha^-1; columns juvenile, subadult, adult) and
clear the state COTSMod derives from earlier years.
"""
function set_initial_state!(cots_state::CotsState, state::AbstractMatrix{<:Real})::Nothing
    size(state) == (length(cots_state), 3) ||
        throw(DimensionMismatch("COTS initial state must be locations x 3 stages"))
    for (model, row) in zip(cots_state, eachrow(state))
        model.N .= row
        model.large_adults =
            model.params.size_weighted_fecundity ?
            row[3] * model.params.size_initial_large_fraction : 0.0
        model.last_effective_breeders = 0.0
        model.last_small_adult_growth = 0.0
        model.burden_ha = 0.0
        model.last_burden_ha = 0.0
        model.last_hazard = 0.0
        model.last_hazard_deaths_ha = 0.0
        # Spread adults over the stable age profile (when enabled) and clear food memory
        COTSMod.reset_adult_ages!(model)
    end
    return nothing
end

"""
    apply_larval_survival!(cots_state, multipliers::AbstractVector{<:Real})

Scale each location's larval production (separated recruitment) by its survival
multiplier, at the source and before dispersal.
"""
function apply_larval_survival!(
    cots_state::CotsState, multipliers::AbstractVector{<:Real}
)::Nothing
    length(multipliers) == length(cots_state) ||
        throw(DimensionMismatch("Larval-survival multipliers must match locations"))
    for (model, m) in zip(cots_state, multipliers)
        model.last_local_fecundity *= m
    end
    return nothing
end

"""
    apply_settlement_modifier!(cots_state, multipliers::AbstractVector{<:Real})

Scale each location's settlement (separated recruitment) by its habitat multiplier, after
the COTS year and before dispersal. Larvae settling later in the year (local, from other
locations and from outside) pass through the scaled settlement gate; recruits already
settled from background immigration are scaled here.
"""
function apply_settlement_modifier!(
    cots_state::CotsState, multipliers::AbstractVector{<:Real}
)::Nothing
    length(multipliers) == length(cots_state) ||
        throw(DimensionMismatch("Settlement modifiers must match locations"))
    for (model, m) in zip(cots_state, multipliers)
        model.last_settlement_gate *= m
        model.N[1] *= m
        model.last_background_immigration *= m
    end
    return nothing
end
