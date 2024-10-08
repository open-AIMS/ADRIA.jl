"""
    SeedCriteriaWeights <: DecisionWeights

Criteria weights for seeding interventions.
"""
Base.@kwdef struct SeedCriteriaWeights <: DecisionWeights
    seed_heat_stress::Param = Factor(
        1.0;
        ptype="continuous",
        dist=Uniform,
        dist_params=(0.8, 1.0),
        direction=minimum,
        name="Seed Heat Stress",
        description="Importance of avoiding heat stress when seeding. Prefer locations with lower heat stress."
    )
    seed_wave_stress::Param = Factor(
        0.3;
        ptype="continuous",
        dist=Uniform,
        dist_params=(0.0, 1.0),
        direction=maximum,
        name="Seed Wave Stress",
        description="Prefer locations with higher wave activity."
    )
    seed_in_connectivity::Param = Factor(
        0.85;
        ptype="continuous",
        dist=Uniform,
        dist_params=(0.5, 1.0),
        direction=maximum,
        name="Incoming Connectivity (Seed)",
        description="Give preference to locations with high incoming connectivity (i.e., receives larvae from other sites) for coral deployments."
    )
    seed_out_connectivity::Param = Factor(
        0.90;
        ptype="continuous",
        dist=Uniform,
        dist_params=(0.5, 1.0),
        direction=maximum,
        name="Outgoing Connectivity (Seed)",
        description="Give preference to locations with high outgoing connectivity (i.e., provides larvae to other sites) for coral deployments."
    )
    seed_depth::Param = Factor(
        0.95;
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
        cluster_ids::Vector{<:Union{Int64,String,Symbol}},
        area_to_seed::Float64,
        considered_locs::Vector{<:Union{Int64,String,Symbol}},
        available_space::Vector{Float64},
        min_locs::Int64,
        max_members::Int64
    )::Vector{<:Union{String,Symbol,Int64}}

Selection locations for seeding deployments.
Attempts to spread deployment locations based on location clusters.

Note: If a subset of locations is passed in, the returned indices
will relate to the position of the subset, not the canonical dataset.

# Arguments
- `sp` : SeedPreferences
- `dm` : The decision matrix to assess pertaining to the locations to be considered
- `method` : MCDA method from JMcDM.jl
- `cluster_ids` : Cluster membership for *all* locations
- `area_to_seed` : total area to be seeded in absolute units (n_corals * mean juvenile area)
- `considered_locs` : Locations to assess, specified by name or indices
- `available_space` : available space for *all* location in absolute units (e.g., m²)
- `min_locs` : Minimum number of locations to consider
- `max_members` : Maximum number of deployment locations per cluster

# Returns
Selected location names/ids in order of their ranks.

# Examples
```julia
dom = ADRIA.load_domain("... path to domain ...")

seed_pref = ADRIA.SeedPreferences(dom)

location_names = cluster_ids
cluster_ids = dom.loc_data.cluster_id
dm = ADRIA.decision_matrix(location_names, seed_pref.names, rand(3806, 6))
available_space = dom.loc_data.k .* dom.loc_data.area
mcda_method = ADRIA.mcda_methods()[1]
min_locs = 5
max_members = 2

ADRIA.select_locations(seed_pref, dm, mcda_method, cluster_ids, 2.0, collect(1:3806), available_space, min_locs, max_members)
```
"""
function select_locations(
    sp::SeedPreferences,
    dm::YAXArray,
    method::Union{Function,DataType},
    cluster_ids::Vector{<:Union{Int64,String,Symbol}},
    area_to_seed::Float64,
    considered_locs::Vector{<:Union{Int64,String,Symbol}},
    available_space::Vector{Float64},
    min_locs::Int64,
    max_members::Int64
)::Vector{<:Union{String,Symbol,Int64}}
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

    if length(considered_locs) == 0
        return String[]
    end

    # Disperse selected locations to avoid "clumping" deployment locations
    dispersed_rank_order, _, n_locs = disperse_locations(
        rank_ordered_idx,
        cluster_ids[considered_locs][rank_ordered_idx],  # Reorder to match ranked location order
        available_space[considered_locs][rank_ordered_idx],
        area_to_seed,
        min_locs,
        max_members
    )

    loc_names = collect(getAxis(:location, dm))

    return loc_names[dispersed_rank_order][1:n_locs]
end

"""
    disperse_locations(ranked_locs::Vector, cluster_ids::Vector, available_space::Vector, area_to_seed::Float64, n_iv_locs::Int64, max_members::Int64; max_iter=5)

Disperse selected deployment locations, ensuring deployment locations are not "clumped"
together.

# Examples
```julia
ranked_locs = string.(collect(1:15))
available_space = rand(15)
cluster_ids = [1,4,4,4,4,3,6,2,5,5,6,6,7,7,1]
area_to_seed = 5.0
n_iv_locs = 5
max_members = 2
disperse_locations(ranked_locs, cluster_ids, available_space, area_to_seed, n_iv_locs, max_members)
```

# Arguments
- `ranked_locs` : Name of locations in order of their rank
- `cluster_ids` : the associated cluster for each location, by location rank order
- `available_space` : available space for each location in order of their ranks
- `area_to_seed` : area to be seeded. Must be in same unit as `available_space`.
- `n_iv_locs` : minimum number of intervention locations
- `max_members` : maximum number of allowable locations per cluster
- `max_iter` : maximum number of attempts before giving up

# Returns
Tuple{Vector,Vector,Int64}
- Locations in preferred order
- their corresponding cluster ids, and
- the number of locations selected
"""
function disperse_locations(
    ranked_locs::Vector,
    cluster_ids::Vector,
    available_space::Vector,
    area_to_seed::Float64,
    n_iv_locs::Int64,
    max_members::Int64;
    max_iter=3
)::Tuple{Vector,Vector,Int64}
    # Only expand the number of locations considered if there is not enough space for
    # deployed corals.
    enough_space = findfirst(>=(area_to_seed), cumsum(available_space[ranked_locs]))

    num_locs = n_iv_locs
    if isnothing(enough_space) || (enough_space < n_iv_locs)
        # If no combinations >= area_to_seed then consider all potential locations
        num_locs = min(length(ranked_locs), n_iv_locs)
    end

    # Need to ensure max_members rule is not breached
    # Count the number of times each selected cluster appears
    # then identify clusters that breach the max membership rule
    selected_clusters, cluster_frequency,
    rule_violators_idx, exceeded_clusters,
    potential_alternatives = _update_state(cluster_ids, num_locs, max_members)

    # If no cluster breaches the rule, then nothing to do!
    if length(rule_violators_idx) == 0
        return ranked_locs, cluster_ids, num_locs
    end

    # Otherwise, swap out locations that violate the rule
    # attempt it `max_iter` times, then abort
    c = 0
    while length(exceeded_clusters) > 0
        for rule_violator in exceeded_clusters
            if !(rule_violator in keys(cluster_frequency))
                continue
            end

            # For each cluster that violates the rule find out how much the rule was
            # violated by. The difference is how many alternate locations we need to find.
            freq = cluster_frequency[rule_violator]
            reduce_by = freq - max_members

            for _ in 1:reduce_by
                # Identify viable clusters that do not breach the rule
                alternates = vcat(
                    fill(false, num_locs),
                    (cluster_ids .∈ Ref(potential_alternatives))[(num_locs + 1):end]
                )

                if count(alternates) == 0
                    # No other options available
                    @debug "No other alternatives. Some grouping remains."
                    return ranked_locs, cluster_ids, num_locs
                end

                # Swap the worst/lowest location for the cluster that breaches the rule
                # for the next best location
                worst_loc_idx = findlast(selected_clusters .== rule_violator)
                next_best_loc_idx = findfirst(alternates .> 0)

                # Swap selection for corresponding location, cluster ids and available space
                worst_loc = ranked_locs[worst_loc_idx]
                ranked_locs[worst_loc_idx] = ranked_locs[next_best_loc_idx]
                ranked_locs[next_best_loc_idx] = worst_loc

                worst_member_id = cluster_ids[worst_loc_idx]
                cluster_ids[worst_loc_idx] = cluster_ids[next_best_loc_idx]
                cluster_ids[next_best_loc_idx] = worst_member_id

                space_at_worst_loc = available_space[worst_loc_idx]
                available_space[worst_loc_idx] = available_space[next_best_loc_idx]
                available_space[next_best_loc_idx] = space_at_worst_loc

                # Swapping out a location may include a location with not much space
                # so we reconsider how many locations we need
                num_locs = max(
                    findfirst(>=(area_to_seed), cumsum(available_space[ranked_locs])),
                    n_iv_locs
                )

                # Update state
                selected_clusters, cluster_frequency,
                rule_violators_idx, exceeded_clusters,
                potential_alternatives = _update_state(cluster_ids, num_locs, max_members)
            end

            if length(exceeded_clusters) == 0
                # Break out of for loop early if solution found
                break
            end
        end

        c += 1
        if c > max_iter
            @debug "Could not reduce clusters down to max_members. Exceeded `max_iter`."
            break
        end
    end

    return ranked_locs, cluster_ids, num_locs
end

"""
    _update_state(cluster_ids::Vector, num_locs::Int64, max_members::Int64)

Update list of:
- selected clusters
- number of locations selected for each cluster
- the indices of clusters which exceed the membership rule
- IDs of the above
- and potential suitable alternative locations

# Arguments
- `cluster_ids` : ID of clusters
- `num_locs` : Number of locations to select
- `max_members` : Maximum number of members per cluster
"""
function _update_state(cluster_ids::Vector, num_locs::Int64, max_members::Int64)
    selected_clusters = cluster_ids[1:min(num_locs, length(cluster_ids))]
    cluster_frequency = countmap(selected_clusters)
    rule_violators_idx = findall(values(cluster_frequency) .> max_members)
    exceeded_clusters = collect(keys(cluster_frequency))[rule_violators_idx]
    potential_alternatives = cluster_ids[cluster_ids .∉ Ref(exceeded_clusters)]

    return selected_clusters,
    cluster_frequency,
    rule_violators_idx,
    exceeded_clusters,
    potential_alternatives
end
