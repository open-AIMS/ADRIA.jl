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
        allee_threshold=Float64(get(params, :allee_threshold, 1.0))
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

function inject_upstream_pulse!(cots_models::CotsState, pulse_locs::Set{Int}; pulse_val::Float64=2.0)::Nothing
    return COTSMod.inject_upstream_pulse!(cots_models, pulse_locs; pulse_val=pulse_val)
end

function inject_upstream_pulse!(cots_models::CotsState, pulse_locs::Set{Int}, pulse_vals::AbstractVector{Float64})::Nothing
    return COTSMod.inject_upstream_pulse!(cots_models, pulse_locs, pulse_vals)
end

function initialize_cots(n_locs::Int, params; enabled::Bool=true, spatial_initial_density::Union{AbstractVector{<:Real},Nothing}=nothing, seed_locs::Union{Set{Int},Nothing}=nothing, init_density::Float64=0.1, initial_multiplier::Float64=parse(Float64, get(ENV, "COTS_INITIAL_MULTIPLIER", "1.5")), rng=Random.GLOBAL_RNG)::CotsState
    return COTSMod.initialize_cots(n_locs, _cots_runtime_params(params); enabled=enabled, spatial_initial_density=spatial_initial_density, seed_locs=seed_locs, init_density=init_density, initial_multiplier=initial_multiplier, rng=rng)
end

function apply_predation!(coral_cover::AbstractArray{Float64,3}, cots_state::CotsState, prey_map::CotsPreyMap)::Nothing
    return COTSMod.apply_predation!(coral_cover, cots_state, prey_map)
end

function disperse_larvae!(cots_state::CotsState, conn::SparseMatrixCSC{Float64,Int64}; scalar::Float64=1.0)::Nothing
    return COTSMod.disperse_larvae!(cots_state, conn; scalar=scalar)
end

function apply_external_supply!(cots_state::CotsState, pulse_locs::Set{Int}; pulse_val::Float64=2.0)::Nothing
    return COTSMod.apply_external_supply!(cots_state, pulse_locs; pulse_val=pulse_val)
end

function apply_external_supply!(cots_state::CotsState, pulse_locs::Set{Int}, pulse_vals::AbstractVector{Float64})::Nothing
    return COTSMod.apply_external_supply!(cots_state, pulse_locs, pulse_vals)
end