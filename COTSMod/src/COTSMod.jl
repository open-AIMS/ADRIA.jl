module COTSMod

using Random
using SparseArrays
using StaticArrays

export AbstractCOTSModel,
    COTSParams,
    COTSHuman,
    COTSState,
    COTSPreyMap,
    COTSPreyState,
    cots_timestep!,
    init_cots_populations,
    init_cots_from_spatial,
    cots_mortality!,
    disperse_cots_larvae!,
    inject_upstream_pulse!,
    initialize_cots,
    apply_predation!,
    disperse_larvae!,
    apply_external_supply!

abstract type AbstractCOTSModel end

Base.@kwdef struct COTSParams
    a::Float64 = 1.5
    b::Float64 = 0.5
    IMM::Float64 = 0.002
    p_tilde::Float64 = 0.97
    C_max::Float64 = 0.8
    m1::Float64 = 0.4
    m2::Float64 = 0.2
    m3::Float64 = 0.1
    a_F::Float64 = 0.6
    a_S::Float64 = 0.15
    h::Float64 = 0.0
    eta_F::Float64 = 1.0
    eta_S::Float64 = 1.0
    eta_starve::Float64 = 2.0
    eta_imm::Float64 = 2.0
    imm_threshold::Float64 = 0.35
    fecundity_gate::Bool = false
    a_ricker::Float64 = 6.0
    b_ricker::Float64 = 0.1
    tau_condition::Float64 = 5.0
    allee_threshold::Float64 = 1.0
end

struct COTSPreyMap
    fast_group_indices::Vector{Int}
    slow_group_indices::Vector{Int}
end

struct COTSPreyState
    fast_cover::Vector{Float64}
    slow_cover::Vector{Float64}
end

mutable struct COTSHuman <: AbstractCOTSModel
    N::MVector{3,Float64}
    body_condition::Float64
    params::COTSParams
end

const COTSState = Vector{COTSHuman}

function cots_timestep!(model::COTSHuman, F::Float64, S::Float64)
    p = model.params
    N = model.N

    F = isnan(F) || isinf(F) ? 0.0 : max(0.0, F)
    S = isnan(S) || isinf(S) ? 0.0 : max(0.0, S)

    alpha = 1.0 / p.tau_condition
    food_signal = min(1.0, (F + S) / (p.C_max * 0.5))
    model.body_condition = clamp((1.0 - alpha) * model.body_condition + alpha * food_signal, 0.0, 1.0)

    starve_threshold = p.C_max * 0.15
    f = if (F + S) > starve_threshold
        1.0
    else
        frac = (F + S) / starve_threshold
        (1.0 - p.p_tilde) + p.p_tilde * frac^3
    end

    fecundity = p.a_ricker * model.body_condition^2
    allee_effect = (N[3]^2) / (p.allee_threshold^2 + N[3]^2)
    R = fecundity * N[3] * exp(-p.b_ricker * N[3]) * allee_effect

    s1 = 1.0 - p.m1
    s2 = 1.0 - p.m2
    s3 = 1.0 - p.m3
    imm_gate = min(1.0, ((F + S) / (p.imm_threshold * p.C_max))^p.eta_imm)
    imm = p.IMM * imm_gate

    N_next = MVector{3,Float64}(0.0, 0.0, 0.0)
    N_next[3] = N[2] * (s2 * f) + N[3] * (s3 * f)
    N_next[2] = N[1] * s1
    N_next[1] = R + imm

    for i in 1:3
        val = N_next[i]
        N_next[i] = isnan(val) || isinf(val) ? 0.0 : max(0.0, val)
    end
    model.N .= N_next

    F_pow = F^p.eta_F
    S_pow = S^p.eta_S
    denom = 1.0 + p.h * (p.a_F * F_pow + p.a_S * S_pow)
    return N[3] * (p.a_F * F_pow) / denom, N[3] * (p.a_S * S_pow) / denom
end

function disperse_cots_larvae!(
    cots_models::COTSState,
    conn::SparseMatrixCSC{Float64,Int64};
    immigration_scalar::Float64=1.0
)::Nothing
    n_locs = length(cots_models)
    local_recruits = [cots_models[loc].N[1] for loc in 1:n_locs]

    for sink in 1:n_locs
        immigration = 0.0
        for k in SparseArrays.nzrange(conn, sink)
            src = SparseArrays.rowvals(conn)[k]
            src == sink && continue
            p = SparseArrays.nonzeros(conn)[k]
            val = local_recruits[src] * p * immigration_scalar
            immigration += isnan(val) || isinf(val) ? 0.0 : val
        end
        total = local_recruits[sink] + immigration
        cots_models[sink].N[1] = isnan(total) || isinf(total) ? 0.0 : max(0.0, total)
    end
    return nothing
end

function init_cots_populations(
    n_locs::Int,
    params::COTSParams;
    outbreak_fraction::Float64=0.25,
    seed_locs::Union{Set{Int},Nothing}=nothing,
    init_density::Float64=0.1,
    rng=Random.GLOBAL_RNG
)::COTSState
    seeded_locs = if !isnothing(seed_locs)
        seed_locs
    else
        n_seeded = max(1, round(Int, n_locs * outbreak_fraction))
        Set(randperm(rng, n_locs)[1:n_seeded])
    end

    return [
        COTSHuman(
            loc in seeded_locs ? MVector{3,Float64}(0.02, 0.02, init_density) : MVector{3,Float64}(0.0, 0.0, 0.0),
            0.8,
            params
        ) for loc in 1:n_locs
    ]
end

function init_cots_from_spatial(
    n_locs::Int,
    params::COTSParams,
    probabilities::AbstractVector{<:Real};
    initial_multiplier::Float64=1.5
)::COTSState
    @assert length(probabilities) == n_locs "Probability vector length must match n_locs"
    models = Vector{COTSHuman}(undef, n_locs)
    for loc in 1:n_locs
        p = Float64(probabilities[loc])
        N_init = if p > 0.1
            scale = min(p / 0.5, 1.0)
            MVector{3,Float64}(0.6 * scale * initial_multiplier, 0.3 * scale * initial_multiplier, 0.1 * scale * initial_multiplier)
        else
            MVector{3,Float64}(0.0, 0.0, 0.0)
        end
        models[loc] = COTSHuman(N_init, 0.8, params)
    end
    return models
end

function cots_mortality!(
    C_cover_t::AbstractArray{Float64,3},
    cots_models::COTSState,
    prey_map::COTSPreyMap
)::Nothing
    _, n_sizes, n_locs = size(C_cover_t)
    @inbounds for loc in 1:n_locs
        F_cover = 0.0
        for g in prey_map.fast_group_indices, s in 1:n_sizes
            F_cover += C_cover_t[g, s, loc]
        end
        S_cover = 0.0
        for g in prey_map.slow_group_indices, s in 1:n_sizes
            S_cover += C_cover_t[g, s, loc]
        end

        Cons_F, Cons_S = cots_timestep!(cots_models[loc], F_cover, S_cover)
        Cons_F = min(Cons_F, F_cover)
        Cons_S = min(Cons_S, S_cover)
        surv_F = F_cover > 0.0 ? max(0.0, 1.0 - Cons_F / F_cover) : 1.0
        surv_S = S_cover > 0.0 ? max(0.0, 1.0 - Cons_S / S_cover) : 1.0
        surv_F = isnan(surv_F) || isinf(surv_F) ? 1.0 : clamp(surv_F, 0.0, 1.0)
        surv_S = isnan(surv_S) || isinf(surv_S) ? 1.0 : clamp(surv_S, 0.0, 1.0)

        for g in prey_map.fast_group_indices, s in 1:n_sizes
            C_cover_t[g, s, loc] *= surv_F
        end
        for g in prey_map.slow_group_indices, s in 1:n_sizes
            C_cover_t[g, s, loc] *= surv_S
        end
    end
    clamp!(C_cover_t, 0.0, 1.0)
    return nothing
end

function inject_upstream_pulse!(cots_models::COTSState, pulse_locs::Set{Int}; pulse_val::Float64=2.0)::Nothing
    for loc in pulse_locs
        loc <= length(cots_models) && (cots_models[loc].N[1] += pulse_val)
    end
    return nothing
end

function inject_upstream_pulse!(cots_models::COTSState, pulse_locs::Set{Int}, pulse_vals::AbstractVector{Float64})::Nothing
    for loc in pulse_locs
        loc <= length(cots_models) && loc <= length(pulse_vals) && (cots_models[loc].N[1] += pulse_vals[loc])
    end
    return nothing
end

function initialize_cots(
    n_locs::Int,
    params::COTSParams;
    enabled::Bool=true,
    spatial_initial_density::Union{AbstractVector{<:Real},Nothing}=nothing,
    seed_locs::Union{Set{Int},Nothing}=nothing,
    init_density::Float64=0.1,
    initial_multiplier::Float64=1.5,
    rng=Random.GLOBAL_RNG
)::COTSState
    !enabled && return [COTSHuman(MVector{3,Float64}(0.0, 0.0, 0.0), 0.8, params) for _ in 1:n_locs]
    !isnothing(spatial_initial_density) && return init_cots_from_spatial(n_locs, params, spatial_initial_density; initial_multiplier=initial_multiplier)
    return init_cots_populations(n_locs, params; seed_locs=seed_locs, init_density=init_density, rng=rng)
end

apply_predation!(coral_cover::AbstractArray{Float64,3}, cots_state::COTSState, prey_map::COTSPreyMap)::Nothing =
    cots_mortality!(coral_cover, cots_state, prey_map)

disperse_larvae!(cots_state::COTSState, conn::SparseMatrixCSC{Float64,Int64}; scalar::Float64=1.0)::Nothing =
    disperse_cots_larvae!(cots_state, conn; immigration_scalar=scalar)

apply_external_supply!(cots_state::COTSState, pulse_locs::Set{Int}; pulse_val::Float64=2.0)::Nothing =
    inject_upstream_pulse!(cots_state, pulse_locs; pulse_val=pulse_val)

apply_external_supply!(cots_state::COTSState, pulse_locs::Set{Int}, pulse_vals::AbstractVector{Float64})::Nothing =
    inject_upstream_pulse!(cots_state, pulse_locs, pulse_vals)

end