# COTS (Crown-of-Thorns Starfish) predation submodel for ADRIA.
#
# Implements a stage-structured COTS population model that feeds back on coral
# cover through predation mortality. COTS preferentially consume fast-growing
# Acropora species while consuming slow-growing massive corals at a lower rate.
#
# The model tracks three age classes per location:
#   - Age 0: Recruits (larvae that have settled)
#   - Age 1: Juveniles
#   - Age 2+: Adults (the predatory stage)
#
# Adult survival is modulated by coral availability (starvation feedback),
# creating natural predator-prey oscillation dynamics.

abstract type AbstractCotsModel end

# Maps ADRIA functional groups to COTS prey categories.
# Fast-growing corals (Acropora spp.) are preferentially consumed.
# Slow-growing corals (massives) are consumed at a lower rate.
#
# In ADRIA with 5 functional groups:
#   - Fast: [1, 2, 3] = Tabular Acropora, Corymbose Acropora, Corymbose non-Acropora
#   - Slow: [4, 5]    = Small massives, Large massives
struct CotsPreyMap
    fast_group_indices::Vector{Int}
    slow_group_indices::Vector{Int}
end

# Stage-structured COTS population model with food-modulated starvation mortality.
#
# Fields:
#   N      - MVector{3, Float64}: Population density by age class [Age 0, Age 1, Age 2+]
#   params - NamedTuple with model parameters:
#     a       : Beverton-Holt recruitment parameter
#     b       : Beverton-Holt density dependence
#     IMM     : Background immigration rate (recruits per timestep)
#     p_tilde : Maximum fraction of mortality attributable to starvation (0-1)
#     C_max   : Coral cover level at which starvation mortality is zero (relative cover)
#     m1      : Age-0 mortality rate
#     m2      : Age-1 mortality rate
#     m3      : Adult (age 2+) mortality rate
#     a_F     : Fast coral consumption rate coefficient
#     a_S     : Slow coral consumption rate coefficient
mutable struct CotsHuman <: AbstractCotsModel
    N::MVector{3, Float64}  # Age 0 (recruits), Age 1 (juveniles), Age 2+ (adults)
    body_condition::Float64 # Maternal nutrition index (0-1), gates larval production
    params::NamedTuple
end

# Run one timestep of the COTS population model.
# Updates the COTS population state in-place and returns the consumption
# of fast and slow coral cover.
#
# Arguments:
#   model - CotsHuman model instance (modified in-place)
#   F     - Total relative coral cover for fast-growing prey types
#   S     - Total relative coral cover for slow-growing prey types
#
# Returns:
#   Tuple (Cons_F, Cons_S): consumption of fast and slow coral cover (relative units).
function cots_timestep!(model::CotsHuman, F::Float64, S::Float64)
    p = model.params
    N = model.N

    # Check for NaNs or negative numbers
    if isnan(F) || isnan(S) || any(isnan, N)
        # println("[DEBUG COTS] NaN detected: F=$F, S=$S, N=$N")
    end
    if F < 0.0 || S < 0.0 || any(x -> x < 0.0, N)
        # Avoid flood but print first occurrences
        # println("[DEBUG COTS] Negative detected: F=$F, S=$S, N=$N")
    end

    # Ensure coral cover is non-negative and finite to prevent complex NaNs during fractional exponentiation
    F = isnan(F) || isinf(F) ? 0.0 : max(0.0, F)
    S = isnan(S) || isinf(S) ? 0.0 : max(0.0, S)

    # --- Body Condition Update ---
    # Exponential moving average of food availability (maternal nutrition)
    tau = get(p, :tau_condition, 3.0)
    alpha = 1.0 / tau
    # Normalise food availability (0.5 * C_max represents "well-fed")
    food_signal = min(1.0, (F + S) / (p.C_max * 0.5))
    new_condition = (1.0 - alpha) * model.body_condition + alpha * food_signal
    model.body_condition = clamp(new_condition, 0.0, 1.0)

    # --- Threshold Starvation ---
    # Mortality stays low while prey is available, crashes drastically below threshold
    starve_threshold = p.C_max * 0.15
    if (F + S) > starve_threshold
        f = 1.0
    else
        frac = (F + S) / starve_threshold
        f = (1.0 - p.p_tilde) + p.p_tilde * frac^3
    end

    # --- Ricker Recruitment ---
    # Allows for overcompensation (massive larval pulses) gated by maternal nutrition
    fecundity = get(p, :a_ricker, 4.0) * model.body_condition^2
    b_ricker = get(p, :b_ricker, 0.8)
    
    # Allee effect: fertilization success crashes when population is too sparse
    allee_A = get(p, :allee_threshold, 1.0)
    allee_effect = (N[3]^2) / (allee_A^2 + N[3]^2)
    
    R = fecundity * N[3] * exp(-b_ricker * N[3]) * allee_effect

    # Base survival rates
    s1 = 1.0 - p.m1
    s2 = 1.0 - p.m2
    s3 = 1.0 - p.m3

    # Coral-gated immigration
    imm_threshold = get(p, :imm_threshold, 0.35)
    eta_imm = get(p, :eta_imm, 2.0)
    imm_gate = min(1.0, ((F + S) / (imm_threshold * p.C_max))^eta_imm)
    imm = p.IMM * imm_gate

    # State update
    N_next = MVector{3, Float64}(0.0, 0.0, 0.0)
    N_next[3] = N[2] * (s2 * f) + N[3] * (s3 * f)
    N_next[2] = N[1] * s1
    N_next[1] = R + imm

    # Ensure all elements of N_next are finite and non-negative
    for i in 1:3
        val = N_next[i]
        N_next[i] = isnan(val) || isinf(val) ? 0.0 : max(0.0, val)
    end

    model.N .= N_next

    # Generalized Holling Type II/III Functional Response
    h = get(p, :h, 0.0)
    eta_F = get(p, :eta_F, 1.0)
    eta_S = get(p, :eta_S, 1.0)

    # Use exponents to support Type III sigmoidal drop-off
    F_pow = F^eta_F
    S_pow = S^eta_S
    denom = 1.0 + h * (p.a_F * F_pow + p.a_S * S_pow)

    Cons_F = N[3] * (p.a_F * F_pow) / denom
    Cons_S = N[3] * (p.a_S * S_pow) / denom

    return Cons_F, Cons_S
end

# Disperse COTS larvae (age-0 recruits) between locations using a connectivity matrix.
# Local Ricker recruitment is retained at the source reef. Cross-reef immigration
# from other locations is added on top, scaled by `immigration_scalar` to represent
# the massive larval production of COTS (each adult female produces ~60 million eggs).
#
# The connectivity matrix was calibrated for coral larval transport. COTS produce
# orders of magnitude more larvae per unit biomass than corals, so `immigration_scalar`
# bridges that gap. Self-connectivity (conn[i,i]) is skipped because local recruitment
# is already handled by the Ricker model in `cots_timestep!`.
#
# Arguments:
#   cots_models        - Vector of CotsHuman models, one per location
#   conn               - Sparse connectivity matrix [sources x sinks]
#   immigration_scalar - Multiplier for cross-reef larval immigration (default 1.0).
#                        Higher values represent greater COTS larval production relative
#                        to the coral baseline the connectivity matrix was built for.
#
# TODO: Replace coral connectivity matrix with a COTS-specific connectivity matrix
# when available. The coral connectivity is used as a reasonable starting proxy
# for COTS larval dispersal.
function disperse_cots_larvae!(
    cots_models::Vector{CotsHuman},
    conn::SparseMatrixCSC{Float64, Int64};
    immigration_scalar::Float64=1.0
)
    n_locs = length(cots_models)

    # Snapshot current recruits (age-0) from all locations before dispersal
    local_recruits = [cots_models[loc].N[1] for loc in 1:n_locs]

    for sink in 1:n_locs
        # Immigration from OTHER reefs only (self-connectivity skipped;
        # local recruitment is already captured by the Ricker model)
        immigration = 0.0
        for k in SparseArrays.nzrange(conn, sink)
            src = SparseArrays.rowvals(conn)[k]
            if src == sink
                continue  # Local recruitment handled separately
            end
            p = SparseArrays.nonzeros(conn)[k]

            val = local_recruits[src] * p * immigration_scalar
            immigration += isnan(val) || isinf(val) ? 0.0 : val
        end

        # Keep full local recruitment + add cross-reef immigration
        total = local_recruits[sink] + immigration
        cots_models[sink].N[1] = isnan(total) || isinf(total) ? 0.0 : max(0.0, total)
    end

    return nothing
end

# Initialize COTS populations across locations. Approximately `outbreak_fraction`
# of locations are seeded with non-zero COTS densities.
#
# Arguments:
#   n_locs            - Number of reef locations
#   params            - NamedTuple of COTS parameters
#   outbreak_fraction - Fraction of locations to seed with initial COTS (default 0.25)
#   seed_locs         - Optional explicit set of location indices to seed (overrides outbreak_fraction)
#   init_density      - Initial adult COTS density at seeded locations (default 0.1)
#   rng               - Random number generator for reproducibility
#
# Returns:
#   Vector{CotsHuman} of length n_locs
#
# TODO: Replace random initialization with empirical COTS distribution data
# when provided. The user will supply a spatial distribution of observed
# COTS outbreaks for initialization.
function init_cots_populations(
    n_locs::Int, params::NamedTuple;
    outbreak_fraction::Float64=0.25,
    seed_locs::Union{Set{Int}, Nothing}=nothing,
    init_density::Float64=0.1,
    rng=Random.GLOBAL_RNG
)
    models = Vector{CotsHuman}(undef, n_locs)

    # Use explicit seed locations if provided, otherwise randomly select
    seeded_locs = if !isnothing(seed_locs)
        seed_locs
    else
        n_seeded = max(1, round(Int, n_locs * outbreak_fraction))
        Set(randperm(rng, n_locs)[1:n_seeded])
    end

    for loc in 1:n_locs
        if loc in seeded_locs
            # Seed with moderate initial COTS population
            N_init = MVector{3, Float64}(0.02, 0.02, init_density)
        else
            # Zero initial COTS - immigration will seed these over time
            N_init = MVector{3, Float64}(0.0, 0.0, 0.0)
        end
        # Initialize body condition to 0.8 (well-fed)
        models[loc] = CotsHuman(N_init, 0.8, params)
    end

    return models
end

# Initialize COTS populations from a vector of empirical probabilities/suitabilities.
# Uses locations where probability > 0.7 to seed COTS.
#
# Arguments:
#   n_locs        - Number of reef locations
#   params        - NamedTuple of COTS parameters
#   probabilities - Vector of initial COTS probability/suitability for each location
#
# Returns:
#   Vector{CotsHuman} of length n_locs
function init_cots_from_spatial(
    n_locs::Int, params::NamedTuple, probabilities::Vector{Float64}
)
    @assert length(probabilities) == n_locs "Probability vector length must match n_locs"
    
    models = Vector{CotsHuman}(undef, n_locs)
    for loc in 1:n_locs
        p = probabilities[loc]
        # Scale initial density based on habitat suitability p
        if p > 0.1
            scale = min(p / 0.5, 1.0)
            multiplier_str = get(ENV, "COTS_INITIAL_MULTIPLIER", "1.5")
            multiplier = parse(Float64, multiplier_str)
            N_init = MVector{3, Float64}(0.6 * scale * multiplier, 0.3 * scale * multiplier, 0.1 * scale * multiplier)
        else
            N_init = MVector{3, Float64}(0.0, 0.0, 0.0)
        end
        models[loc] = CotsHuman(N_init, 0.8, params)
    end
    
    return models
end

# Apply COTS predation mortality to coral cover.
# Follows the same pattern as cyclone_mortality!(): modifies C_cover_t
# in-place by reducing cover proportionally based on COTS consumption.
# The consumption is distributed uniformly across all size classes within
# each prey category (fast/slow).
#
# Arguments:
#   C_cover_t   - Coral cover array [n_groups, n_sizes, n_locs] in relative units (0.0-1.0)
#   cots_models - Vector of CotsHuman models, one per location
#   prey_map    - CotsPreyMap mapping functional groups to prey categories
function cots_mortality!(
    C_cover_t::AbstractArray{Float64,3},
    cots_models::Vector{CotsHuman},
    prey_map::CotsPreyMap
)
    n_groups, n_sizes, n_locs = size(C_cover_t)

    @inbounds for loc in 1:n_locs
        # Aggregate relative cover into Fast and Slow prey types
        F_cover = 0.0
        for g in prey_map.fast_group_indices
            for s in 1:n_sizes
                F_cover += C_cover_t[g, s, loc]
            end
        end
        S_cover = 0.0
        for g in prey_map.slow_group_indices
            for s in 1:n_sizes
                S_cover += C_cover_t[g, s, loc]
            end
        end

        # Run COTS population dynamics and get consumption
        Cons_F, Cons_S = cots_timestep!(cots_models[loc], F_cover, S_cover)

        # Cap consumption to available cover
        Cons_F = min(Cons_F, F_cover)
        Cons_S = min(Cons_S, S_cover)

        # Calculate proportional survival (same pattern as cyclone_mortality!)
        surv_F = F_cover > 0.0 ? max(0.0, 1.0 - Cons_F / F_cover) : 1.0
        surv_S = S_cover > 0.0 ? max(0.0, 1.0 - Cons_S / S_cover) : 1.0

        # Robust safeguards against NaNs or Inf
        surv_F = isnan(surv_F) || isinf(surv_F) ? 1.0 : clamp(surv_F, 0.0, 1.0)
        surv_S = isnan(surv_S) || isinf(surv_S) ? 1.0 : clamp(surv_S, 0.0, 1.0)

        # Apply proportional mortality uniformly across all size classes
        for g in prey_map.fast_group_indices
            for s in 1:n_sizes
                C_cover_t[g, s, loc] *= surv_F
            end
        end
        for g in prey_map.slow_group_indices
            for s in 1:n_sizes
                C_cover_t[g, s, loc] *= surv_S
            end
        end
    end

    clamp!(C_cover_t, 0.0, 1.0)
    return nothing
end

# Inject an upstream recruitment pulse (larvae) to specified locations.
#
# Arguments:
#   cots_models - Vector of CotsHuman models, one per location
#   pulse_locs  - Set of location indices to receive the pulse
#   pulse_val   - Larval density to add to age class 0 (default 2.0)
function inject_upstream_pulse!(
    cots_models::Vector{CotsHuman},
    pulse_locs::Set{Int};
    pulse_val::Float64=2.0
)
    for loc in pulse_locs
        if loc <= length(cots_models)
            cots_models[loc].N[1] += pulse_val
        end
    end
    return nothing
end

# Inject location-specific upstream recruitment pulses (larvae).
#
# `pulse_vals` is indexed by location and lets external supply scale with each
# location's expected larval supply rather than adding the same value everywhere.
function inject_upstream_pulse!(
    cots_models::Vector{CotsHuman},
    pulse_locs::Set{Int},
    pulse_vals::AbstractVector{Float64}
)
    for loc in pulse_locs
        if loc <= length(cots_models) && loc <= length(pulse_vals)
            cots_models[loc].N[1] += pulse_vals[loc]
        end
    end
    return nothing
end