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

Build COTSMod's runtime parameters from a scenario's `CotsParams` factors. Mechanisms that
are not exposed as factors keep COTSMod's defaults (off).
"""
function cots_runtime_params(param_set::YAXArray)::COTSMod.COTSParams
    f(name) = Float64(param_set[At(name)])
    on(name) = param_set[At(name)] == 1
    return COTSMod.COTSParams(;
        a=f("a"),
        b=f("b"),
        IMM=f("IMM"),
        p_tilde=f("p_tilde"),
        C_max=f("C_max"),
        m1=f("m1"),
        m2=f("m2"),
        m3=f("m3"),
        a_F=f("a_F"),
        a_S=f("a_S"),
        h=f("h"),
        eta_F=f("eta_F"),
        eta_S=f("eta_S"),
        eta_starve=f("eta_starve"),
        eta_imm=f("eta_imm"),
        imm_threshold=f("imm_threshold"),
        fecundity_gate=on("fecundity_gate"),
        a_ricker=f("a_ricker"),
        b_ricker=f("b_ricker"),
        tau_condition=f("tau_condition"),
        allee_threshold=f("allee_threshold"),
        adult_senescence=on("adult_senescence"),
        senescence_age=round(Int, f("senescence_age")),
        senescence_mortality=f("senescence_mortality"),
        allee_stretched=on("allee_stretched"),
        allee_exponent=f("allee_exponent"),
        juvenile_predation=on("juvenile_predation"),
        juvenile_predation_capacity=f("juvenile_predation_capacity"),
        per_capita_consumption=on("per_capita_consumption"),
        consumption_m2_per_adult_year=f("consumption_m2_per_adult_year"),
        preferred_condition=on("preferred_condition"),
        condition_saturation=f("condition_saturation"),
        condition_memory_years=f("condition_memory_years"),
        preferred_prey_starvation=on("preferred_prey_starvation"),
        starvation_preferred_threshold=f("starvation_preferred_threshold"),
        starvation_memory_years=f("starvation_memory_years")
    )
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
