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

    # Stock-Recruitment (Beverton-Holt)
    R = (p.a * N[3]) / (1.0 + p.b * N[3])

    # Coral-modulated mortality factor (starvation feedback)
    # When coral cover is high (F+S >= C_max), f = 1.0 (no starvation).
    # When coral cover is zero, f = 1 - p_tilde (maximum starvation mortality).
    rho = min(1.0, (F + S) / p.C_max)
    f = (1.0 - p.p_tilde) + p.p_tilde * rho

    # Base survival rates
    s1 = 1.0 - p.m1
    s2 = 1.0 - p.m2
    s3 = 1.0 - p.m3

    # State update
    N_next = MVector{3, Float64}(0.0, 0.0, 0.0)
    N_next[3] = N[2] * (s2 * f) + N[3] * (s3 * f)
    N_next[2] = N[1] * s1
    N_next[1] = R + p.IMM

    model.N .= N_next

    # Consumption (proportional to adult density and available cover)
    Cons_F = N[3] * p.a_F * F
    Cons_S = N[3] * p.a_S * S

    return Cons_F, Cons_S
end

# Disperse COTS larvae (age-0 recruits) between locations using a connectivity matrix.
# Redistributes age-class 0 across sink locations based on connectivity-weighted
# contributions from all source locations.
#
# Arguments:
#   cots_models - Vector of CotsHuman models, one per location
#   conn        - Sparse connectivity matrix [sources x sinks]
#
# TODO: Replace coral connectivity matrix with a COTS-specific connectivity matrix
# when available. The coral connectivity is used as a reasonable starting proxy
# for COTS larval dispersal.
function disperse_cots_larvae!(
    cots_models::Vector{CotsHuman},
    conn::SparseMatrixCSC{Float64, Int64}
)
    n_locs = length(cots_models)

    # Snapshot current recruits (age-0) from all locations before redistribution
    source_recruits = [cots_models[loc].N[1] for loc in 1:n_locs]

    # Normalize connectivity: columns = sink locations, rows = source locations
    col_sums = vec(sum(conn; dims=1))

    for sink in 1:n_locs
        if col_sums[sink] == 0.0
            continue
        end

        # Weighted sum of source recruits arriving at this sink
        incoming = 0.0
        for k in SparseArrays.nzrange(conn, sink)
            src = SparseArrays.rowvals(conn)[k]
            w = SparseArrays.nonzeros(conn)[k] / col_sums[sink]
            incoming += source_recruits[src] * w
        end

        cots_models[sink].N[1] = incoming
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
    rng=Random.GLOBAL_RNG
)
    models = Vector{CotsHuman}(undef, n_locs)
    n_seeded = max(1, round(Int, n_locs * outbreak_fraction))
    seeded_locs = Set(randperm(rng, n_locs)[1:n_seeded])

    for loc in 1:n_locs
        if loc in seeded_locs
            # Seed with moderate initial COTS population
            N_init = MVector{3, Float64}(0.02, 0.02, 0.1)
        else
            # Zero initial COTS - immigration will seed these over time
            N_init = MVector{3, Float64}(0.0, 0.0, 0.0)
        end
        models[loc] = CotsHuman(N_init, params)
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
