"""
    CAqCriteriaWeights <: DecisionWeights

Criteria weights for coral aquaculture interventions.
"""
Base.@kwdef struct CAqCriteriaWeights <: DecisionWeights
    iv_CAq_heat_stress::Param = Factor(
        0.9;
        ptype="continuous",
        dist=Uniform,
        dist_params=(0.8, 1.0),
        direction=minimum,
        name="Coral Aquaculture Heat Stress",
        description="Importance of avoiding heat stress for coral aquaculture. Prefer locations with lower heat stress."
    )
    iv_CAq_wave_stress::Param = Factor(
        0.5;
        ptype="continuous",
        dist=Uniform,
        dist_params=(0.0, 1.0),
        direction=maximum,
        name="Coral Aquaculture Wave Stress",
        description="Prefer locations with higher wave activity."
    )
    iv_CAq_in_connectivity::Param = Factor(
        0.5;
        ptype="continuous",
        dist=Uniform,
        dist_params=(0.2, 1.0),
        direction=maximum,
        name="Incoming Connectivity (Coral Aquaculture)",
        description="Give preference to locations with high incoming connectivity (i.e., receives larvae from other sites) for coral deployments."
    )
    iv_CAq_out_connectivity::Param = Factor(
        0.80;
        ptype="continuous",
        dist=Uniform,
        dist_params=(0.2, 1.0),
        direction=maximum,
        name="Outgoing Connectivity (Coral Aquaculture)",
        description="Give preference to locations with high outgoing connectivity (i.e., provides larvae to other sites) for coral deployments."
    )
    iv_CAq_depth::Param = Factor(
        1.0;
        ptype="continuous",
        dist=Uniform,
        dist_params=(0.8, 1.0),
        direction=maximum,
        name="Depth (Coral Aquaculture)",
        description="Give preference to deeper locations for coral deployments."
    )
    iv_CAq_coral_cover::Param = Factor(
        0.7;
        ptype="continuous",
        dist=Uniform,
        dist_params=(0.0, 1.0),
        direction=minimum,
        name="Coral Aquaculture Coral Cover",
        description="Preference locations with lower coral cover (higher available space) for coral aquaculture deployments."
    )
    iv_CAq_cluster_diversity::Param = Factor(
        0.7;
        ptype="continuous",
        dist=Uniform,
        dist_params=(0.0, 1.0),
        direction=maximum,
        name="Cluster Diversity",
        description="Prefer locations from clusters that are under-represented."
    )
    iv_CAq_geographic_separation::Param = Factor(
        0.8;
        ptype="continuous",
        dist=Uniform,
        dist_params=(0.0, 1.0),
        direction=minimum,
        name="Geographic Separation",
        description="Prefer locations that are distant (when maximized) or closer (when minimized; the default) to their neighbors."
    )
    # Disabled as they are currently unnecessary
    # iv_CAq_priority::Param = Factor(
    #     1.0;
    #     ptype="continuous",
    #     dist=Uniform,
    #     dist_params=(0.0, 1.0),
    #     direction=maximum,
    #     name="Predecessor Priority (Coral Aquaculture)",
    #     description="Preference locations that provide larvae to priority reefs.",
    # )
    # iv_CAq_zone::Param = Factor(
    #     0.0;
    #     ptype="continuous",
    #     dist=Uniform,
    #     dist_params=(0.0, 1.0),
    #     direction=maximum,
    #     name="Zone Predecessor (Coral Aquaculture)",
    #     description="Preference locations that provide larvae to priority (target) zones.",
    # )
end

"""
    CAqPreferences <: DecisionPreference

Preference type specific for coral aquaculture interventions to allow
coral-aquaculture-specific routines to be handled.
"""
struct CAqPreferences <: DecisionPreference
    names::Vector{Symbol}
    weights::Vector{Float64}
    directions::Vector{Function}
end

function CAqPreferences(dom, params::YAXArray)::CAqPreferences
    w::DataFrame = component_params(dom.model, CAqCriteriaWeights)
    # Strip the `iv_CAq_` prefix (2 segments) to recover the bare criterion name
    # (e.g. `heat_stress`) so it matches the generic kwargs used in `setup_guided_intervention`.
    cn = Symbol[Symbol(join(split(string(cn), "_")[3:end], "_")) for cn in w.fieldname]

    return CAqPreferences(cn, params[factors = At(string.(w.fieldname))], w.direction)
end
function CAqPreferences(dom, params...)::CAqPreferences
    w::DataFrame = component_params(dom.model, CAqCriteriaWeights)
    for (k, v) in params
        w[w.fieldname .== k, :val] .= v
    end

    return CAqPreferences(w.fieldname, w.val, w.direction)
end
function CAqPreferences(dom)
    w::DataFrame = component_params(dom.model, CAqCriteriaWeights)
    return CAqPreferences(w.fieldname, w.val, w.direction)
end

"""
    select_locations(
        sp::CAqPreferences,
        dm::YAXArray,
        method::Union{Function,DataType},
        considered_locs::Vector{<:Union{Int64,String,Symbol}},
        min_locs::Int64
    )::Vector{<:Union{String,Symbol,Int64}}

Select locations for coral aquaculture interventions based on multiple criteria,
including spatial distribution.

# Example
```julia
caq_pref = CAqPreferences(domain, param_set)
decision_mat = decision_matrix(domain.loc_ids, caq_pref.names, criteria_values)
mcda_method = mcda_methods()[1]  # Use first method from available MCDA methods
valid_locs = domain.loc_ids

selected_locs = select_locations(
    caq_pref,
    decision_mat,
    mcda_method,
    valid_locs,
    5  # Select at least 5 locations
)
```

# Arguments
- `sp`: CAqPreferences containing criteria names, weights, and optimization directions
- `dm`: Decision matrix with locations as rows and criteria as columns
- `method`: MCDA method to use for ranking (from the JMcDM package)
- `considered_locs`: Vector of location identifiers to consider for selection
- `min_locs`: Minimum number of locations to select

# Returns
Vector of selected location identifiers, ordered by their ranks
"""
function select_locations(
    sp::CAqPreferences,
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
