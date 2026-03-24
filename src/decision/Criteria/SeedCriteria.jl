"""
    SeedCriteriaWeights <: DecisionWeights

Criteria weights for seeding interventions.
"""
Base.@kwdef struct SeedCriteriaWeights <: DecisionWeights
    seed_heat_stress::Param = Factor(
        0.9;
        ptype="continuous",
        dist=Uniform,
        dist_params=(0.8, 1.0),
        direction=minimum,
        name="Seed Heat Stress",
        description="Importance of avoiding heat stress when seeding. Prefer locations with lower heat stress."
    )
    seed_wave_stress::Param = Factor(
        0.5;
        ptype="continuous",
        dist=Uniform,
        dist_params=(0.0, 1.0),
        direction=maximum,
        name="Seed Wave Stress",
        description="Prefer locations with higher wave activity."
    )
    seed_in_connectivity::Param = Factor(
        0.5;
        ptype="continuous",
        dist=Uniform,
        dist_params=(0.5, 1.0),
        direction=maximum,
        name="Incoming Connectivity (Seed)",
        description="Give preference to locations with high incoming connectivity (i.e., receives larvae from other sites) for coral deployments."
    )
    seed_out_connectivity::Param = Factor(
        0.80;
        ptype="continuous",
        dist=Uniform,
        dist_params=(0.5, 1.0),
        direction=maximum,
        name="Outgoing Connectivity (Seed)",
        description="Give preference to locations with high outgoing connectivity (i.e., provides larvae to other sites) for coral deployments."
    )
    seed_depth::Param = Factor(
        1.0;
        ptype="continuous",
        dist=Uniform,
        dist_params=(0.8, 1.0),
        direction=maximum,
        name="Depth (Seed)",
        description="Give preference to deeper locations for coral deployments."
    )
    seed_coral_cover::Param = Factor(
        0.7;
        ptype="continuous",
        dist=Uniform,
        dist_params=(0.0, 1.0),
        direction=minimum,
        name="Seed Coral Cover",
        description="Preference locations with lower coral cover (higher available space) for seeding deployments."
    )
    seed_cluster_diversity::Param = Factor(
        0.7;
        ptype="continuous",
        dist=Uniform,
        dist_params=(0.0, 1.0),
        direction=maximum,
        name="Cluster Diversity",
        description="Prefer locations from clusters that are under-represented."
    )
    seed_geographic_separation::Param = Factor(
        0.8;
        ptype="continuous",
        dist=Uniform,
        dist_params=(0.0, 1.0),
        direction=minimum,
        name="Geographic Separation",
        description="Prefer locations that are distant (when maximized) or closer (when minimized; the default) to their neighbors."
    )
    # Disabled as they are currently unnecessary
    # seed_priority::Param = Factor(
    #     1.0;
    #     ptype="continuous",
    #     dist=Uniform,
    #     dist_params=(0.0, 1.0),
    #     direction=maximum,
    #     name="Predecessor Priority (Seed)",
    #     description="Preference locations that provide larvae to priority reefs.",
    # )
    # seed_zone::Param = Factor(
    #     0.0;
    #     ptype="continuous",
    #     dist=Uniform,
    #     dist_params=(0.0, 1.0),
    #     direction=maximum,
    #     name="Zone Predecessor (Seed)",
    #     description="Preference locations that provide larvae to priority (target) zones.",
    # )
end

"""
    SeedPreferences <: DecisionPreference

Preference type specific for seeding interventions to allow seeding-specific routines to
be handled.
"""
struct SeedPreferences <: DecisionPreference
    names::Vector{Symbol}
    weights::Vector{Float64}
    directions::Vector{Function}
end

function SeedPreferences(dom, params::YAXArray)::SeedPreferences
    w::DataFrame = component_params(dom.model, SeedCriteriaWeights)
    cn = Symbol[Symbol(join(split(string(cn), "_")[2:end], "_")) for cn in w.fieldname]

    return SeedPreferences(cn, params[factors=At(string.(w.fieldname))], w.direction)
end
function SeedPreferences(dom, params...)::SeedPreferences
    w::DataFrame = component_params(dom.model, SeedCriteriaWeights)
    for (k, v) in params
        w[w.fieldname .== k, :val] .= v
    end

    return SeedPreferences(w.fieldname, w.val, w.direction)
end
function SeedPreferences(dom)
    w::DataFrame = component_params(dom.model, SeedCriteriaWeights)
    return SeedPreferences(w.fieldname, w.val, w.direction)
end

"""
    select_locations(
        sp::SeedPreferences,
        dm::YAXArray,
        method::Union{Function,DataType},
        considered_locs::Vector{<:Union{Int64,String,Symbol}},
        min_locs::Int64
    )::Vector{<:Union{String,Symbol,Int64}}

Select locations for seeding interventions based on multiple criteria, including spatial
distribution.

# Example
```julia
seed_pref = SeedPreferences(domain, param_set)
decision_mat = decision_matrix(domain.loc_ids, seed_pref.names, criteria_values)
mcda_method = mcda_methods()[1]  # Use first method from available MCDA methods
valid_locs = domain.loc_ids

selected_locs = select_locations(
    seed_pref,
    decision_mat,
    mcda_method,
    valid_locs,
    5  # Select at least 5 locations
)
```

# Arguments
- `sp`: SeedPreferences containing criteria names, weights, and optimization directions
- `dm`: Decision matrix with locations as rows and criteria as columns
- `method`: MCDA method to use for ranking (from the JMcDM package)
- `considered_locs`: Vector of location identifiers to consider for selection
- `min_locs`: Minimum number of locations to select

# Returns
Vector of selected location identifiers, ordered by their ranks
"""
function select_locations(
    sp::SeedPreferences,
    dm::YAXArray,
    method::Union{Function,DataType},
    considered_locs::Vector{<:Union{Int64,String,Symbol}},
    min_locs::Int64
)::Vector{<:Union{String,Symbol,Int64}}
    if length(considered_locs) == 0
        return String[]
    end

    loc_names = collect(getAxis(:location, dm))

    # Continue with existing ranking process
    local rank_ordered_idx
    try
        rank_ordered_idx = rank_by_index(sp, dm, method)
    catch err
        if err isa DomainError
            # Return empty vector to signify no ranks
            return String[]
        end
        rethrow(err)
    end

    # Take top n_locs from the ranked list
    n_locs = min(min_locs, length(rank_ordered_idx))

    return collect(loc_names[rank_ordered_idx][1:n_locs])
end
