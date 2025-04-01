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
    seed_cluster_diversity::Param = Factor(
        0.5;
        ptype="continuous",
        dist=Uniform,
        dist_params=(0.0, 1.0),
        direction=maximum,
        name="Cluster Diversity",
        description="Prefer locations from clusters that are under-represented."
    )
    seed_geographic_separation::Param = Factor(
        0.5;
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

    if length(considered_locs) == 0
        return String[]
    end

    # Take top n_locs from the ranked list
    n_locs = min(min_locs, length(rank_ordered_idx))

    return collect(loc_names[rank_ordered_idx][1:n_locs])
end

"""
    cluster_diversity(cluster_ids::Vector{<:Union{Int64,String,Symbol}})::Vector{Float64}

Calculate how selecting each location would improve cluster diversity.
Higher scores indicate locations from underrepresented clusters.
"""
function cluster_diversity(
    cluster_ids::Vector{<:Union{Int64,String,Symbol}}
)::Vector{Float64}
    # Count unique clusters and their frequencies
    clusters = unique(cluster_ids)
    n_clusters = length(clusters)
    n_locations = length(cluster_ids)

    # Calculate the ideal count per cluster (perfect distribution)
    ideal_count = n_locations / n_clusters

    # Create a map of cluster to its count
    cluster_counts = countmap(cluster_ids)

    # Calculate diversity score for each unique cluster
    cluster_scores = Dict{eltype(clusters),Float64}()

    for cluster in clusters
        # Calculate how underrepresented this cluster is
        count_ratio = ideal_count / cluster_counts[cluster]

        # Scale so the maximum possible value is 1.0
        # A fully represented cluster gets 0.0
        # An underrepresented cluster gets a score based on its deficit
        if count_ratio > 1.0  # Underrepresented
            cluster_scores[cluster] = min(1.0, (count_ratio - 1.0) / (n_clusters - 1.0))
        else  # Overrepresented
            cluster_scores[cluster] = 0.0
        end
    end

    # Map cluster scores back to each location
    diversity_scores = zeros(n_locations)
    for i in 1:n_locations
        diversity_scores[i] = cluster_scores[cluster_ids[i]]
    end

    return diversity_scores
end

"""
    geographic_separation(mean_distances::Matrix{Float64})::Vector{Float64}

Calculate how spatially separated each location is from others.
Higher scores indicate locations that are more distant from their neighbors.
"""
function geographic_separation(mean_distances::Vector{Float64})::Vector{Float64}
    # Normalize to [0,1]
    max_dist = maximum(mean_distances)
    if max_dist > 0
        return mean_distances ./ max_dist
    end

    return zeros(length(mean_distances))
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
