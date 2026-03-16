Base.@kwdef struct Intervention <: EcoModel
    # Intervention Factors
    # Bounds are defined as floats to maintain type stability
    guided::Param = Factor(
        0;
        ptype="unordered categorical",
        dist=CategoricalDistribution,
        dist_params=(Tuple(-1:Float64(length(decision.mcda_methods())))),
        name="Guided",
        description="Choice of MCDA approach."
    )
    N_seed_TA::Param = Factor(
        0;
        ptype="ordered discrete",
        dist=DiscreteOrderedUniformDist,
        dist_params=(0.0, 1000000.0, 50000.0),  # increase in steps of 50K
        name="Seeded Tabular Acropora",
        description="Number of Tabular Acropora to seed per deployment event."
    )
    N_seed_CA::Param = Factor(
        0;
        ptype="ordered discrete",
        dist=DiscreteOrderedUniformDist,
        dist_params=(0.0, 1000000.0, 50000.0),  # increase in steps of 50K
        name="Seeded Corymbose Acropora",
        description="Number of Corymbose Acropora to seed per deployment event."
    )
    N_seed_CNA::Param = Factor(
        0;
        ptype="ordered discrete",
        dist=DiscreteOrderedUniformDist,
        dist_params=(0.0, 1000000.0, 50000.0),  # increase in steps of 50K
        name="Seeded Corymbose non-Acropora",
        description="Number of Corymbose non-Acropora to seed per deployment event."
    )
    N_seed_SM::Param = Factor(
        0;
        ptype="ordered discrete",
        dist=DiscreteOrderedUniformDist,
        dist_params=(0.0, 1000000.0, 50000.0),  # increase in steps of 50K
        name="Seeded Small Massives",
        description="Number of small massives/encrusting to seed per deployment event."
    )
    N_seed_LM::Param = Factor(
        0;
        ptype="ordered discrete",
        dist=DiscreteOrderedUniformDist,
        dist_params=(0.0, 1000000.0, 50000.0),  # increase in steps of 50K
        name="Seeded Large Massives",
        description="Number of large massives/encrusting to seed per deployment event."
    )
    N_mc_settlers::Param = Factor(
        0;
        ptype="ordered discrete",
        dist=DiscreteOrderedUniformDist,
        dist_params=(0.0, 1000000.0, 50000.0),  # increase in steps of 50K
        name="Moving corals settlers",
        description="Number of moving coral settlers added per deployment event."
    )
    min_iv_locations::Param = Factor(
        5;
        ptype="ordered discrete",
        dist=DiscreteUniform,
        dist_params=(5.0, 20.0),
        name="Minimum moving corals locations",
        description="Minimum number of locations to perform moving corals intervention"
    )
    fogging::Param = Factor(
        0.0;
        ptype="continuous",
        dist=TriangularDist,
        dist_params=(0.0, 0.3, 0.0),
        name="Fogging",
        description="Assumed reduction in bleaching mortality."
    )
    SRM::Param = Factor(
        0.0;
        ptype="continuous",
        dist=TriangularDist,
        dist_params=(0.0, 7.0, 0.0),
        name="SRM",
        description="Reduction in DHWs due to shading."
    )
    a_adapt::Param = Factor(
        0.0;
        ptype="ordered discrete",
        dist=DiscreteOrderedUniformDist,
        dist_params=(0.0, 15.0, 0.5),  # increase in steps of 0.5 DHW enhancement
        name="Assisted Adaptation",
        description="Assisted adaptation in terms of DHW resistance."
    )
    seed_years::Param = Factor(
        10;
        ptype="ordered categorical",
        dist=DiscreteUniform,
        dist_params=(5.0, 75.0),
        name="Years to Seed",
        description="Number of years to seed for."
    )
    shade_years::Param = Factor(
        10;
        ptype="ordered categorical",
        dist=DiscreteUniform,
        dist_params=(5.0, 75.0),
        name="Years to Shade",
        description="Number of years to shade for."
    )
    fog_years::Param = Factor(
        10;
        ptype="ordered categorical",
        dist=DiscreteUniform,
        dist_params=(5.0, 75.0),
        name="Years to fog",
        description="Number of years to fog for."
    )
    plan_horizon::Param = Factor(
        5;
        ptype="ordered categorical",
        dist=DiscreteUniform,
        dist_params=(0.0, 20.0),
        name="Planning Horizon",
        description="How many years of projected data to take into account when selecting intervention locations (0 only accounts for current deployment year)."
    )
    seed_deployment_freq::Param = Factor(
        5;
        ptype="ordered categorical",
        dist=DiscreteUniform,
        dist_params=(0.0, 15.0),
        name="Selection Frequency (Seed)",
        description="Frequency of seeding deployments (0 deploys once)."
    )
    fog_deployment_freq::Param = Factor(
        5;
        ptype="ordered categorical",
        dist=DiscreteUniform,
        dist_params=(0.0, 15.0),
        name="Selection Frequency (Fog)",
        description="Frequency of fogging deployments (0 deploys once)."
    )
    shade_deployment_freq::Param = Factor(
        1;
        ptype="ordered categorical",
        dist=DiscreteUniform,
        dist_params=(1.0, 15.0),
        name="Deployment Frequency (Shading)",
        description="Frequency of shading deployments."
    )
    seed_year_start::Param = Factor(
        2;
        ptype="ordered categorical",
        dist=DiscreteUniform,
        dist_params=(0.0, 25.0),
        name="Seeding Start Year",
        description="Start seeding deployments after this number of years has elapsed."
    )
    shade_year_start::Param = Factor(
        2;
        ptype="ordered categorical",
        dist=DiscreteUniform,
        dist_params=(2.0, 25.0),
        name="Shading Start Year",
        description="Start of shading deployments after this number of years has elapsed."
    )
    fog_year_start::Param = Factor(
        2;
        ptype="ordered categorical",
        dist=DiscreteUniform,
        dist_params=(2.0, 25.0),
        name="Fogging Start Year",
        description="Start of fogging deployments after this number of years has elapsed."
    )
    mc_year_start::Param = Factor(
        2;
        ptype="ordered categorical",
        dist=DiscreteUniform,
        dist_params=(0.0, 25.0),
        name="Moving corals Start Year",
        description="Start moving corals deployments after this number of years has elapsed."
    )
    mc_years::Param = Factor(
        10;
        ptype="ordered categorical",
        dist=DiscreteUniform,
        dist_params=(5.0, 75.0),
        name="Years to deploy moving corals",
        description="Number of years to deploy moving corals."
    )
    mc_deployment_freq::Param = Factor(
        1;
        ptype="ordered categorical",
        dist=DiscreteUniform,
        dist_params=(1.0, 15.0),
        name="Deployment Frequency (Moving corals)",
        description="Frequency of moving corals deployments."
    )
    # Intervention strategy parameters
    seed_strategy::Param = Factor(
        2;
        ptype="ordered categorical",
        dist=CategoricalDistribution,
        dist_params=(1.0, 2.0),
        name="Seed Strategy Type",
        description="Deployment strategy: 1=Periodic (time-based), 2=Reactive (condition-based); 0 is off"
    )
    fog_strategy::Param = Factor(
        2;
        ptype="ordered categorical",
        dist=CategoricalDistribution,
        dist_params=(1.0, 2.0),
        name="Fog Strategy Type",
        description="Deployment strategy: 1=Periodic (time-based), 2=Reactive (condition-based); 0 is off"
    )
    mc_strategy::Param = Factor(
        2;
        ptype="ordered categorical",
        dist=CategoricalDistribution,
        dist_params=(1.0, 2.0),
        name="Moving Corals Strategy Type",
        description="Deployment strategy: 1=Periodic (time-based), 2=Reactive (condition-based); 0 is off"
    )
    reactive_absolute_threshold::Param = Factor(
        0.95;
        ptype="ordered discrete",
        dist=DiscreteOrderedUniformDist,
        dist_params=(0.7, 0.95, 0.05),
        name="Cover Absolute Threshold",
        description="Deploy when coral cover falls below this proportion (for reactive strategy)"
    )
    reactive_loss_threshold::Param = Factor(
        0.30;
        ptype="ordered discrete",
        dist=DiscreteOrderedUniformDist,
        dist_params=(0.1, 0.50, 0.05),
        name="Cover Loss Threshold",
        description="Deploy when proportional cover loss exceeds this value (for reactive strategy)"
    )
    reactive_min_cover_remaining::Param = Factor(
        0.05;
        ptype="ordered discrete",
        dist=DiscreteOrderedUniformDist,
        dist_params=(0.0, 0.15, 0.025),
        name="Minimum Viable Cover",
        description="Do not deploy to locations with less than this cover proportion (for reactive strategy)"
    )
    reactive_response_delay::Param = Factor(
        1.0;
        ptype="ordered categorical",
        dist=DiscreteUniform,
        dist_params=(0.0, 5.0),
        name="Response Delay",
        description="Timesteps to wait after trigger before deployment (for reactive strategy)"
    )
    reactive_cooldown_period::Param = Factor(
        4.0;
        ptype="ordered categorical",
        dist=DiscreteUniform,
        dist_params=(0.0, 10.0),
        name="Reactive Cooldown Period",
        description="Timesteps before location becomes eligible again; 0=no cooldown (for reactive strategy)"
    )
end

function interventions()
    return [:seed, :fog]
end

function year_start_factors(dom::Domain)::DataFrame
    ms = model_spec(dom)
    return ms[occursin.(Ref("year_start"), string.(ms.fieldname)), :]
end

function setup_guided_intervention(
    domain::Domain,
    param_set::YAXArray,
    depth_criteria::BitVector,
    preference,
    target_locs::Vector{String},
    is_intervention::Bool,
    build_strategy::Function
)
    # Remove locations that cannot support corals or are out of depth bounds
    # from consideration
    valid_locs_mask =
        (location_k(domain) .> 0.0) .& depth_criteria .& (domain.loc_ids .∈ [target_locs])

    # Calculate cluster diversity and geographic separation scores
    diversity_scores = cluster_diversity(domain.loc_data.cluster_id)
    separation_scores = geographic_separation(domain.loc_data.mean_to_neighbor)

    pref = preference(domain, param_set)
    decision_mat = decision_matrix(
        domain.loc_ids[valid_locs_mask],
        pref.names;
        depth=domain.loc_data.depth_med[valid_locs_mask],
        cluster_diversity=diversity_scores[valid_locs_mask],
        geographic_separation=separation_scores[valid_locs_mask]
    )
    strategy =
        is_intervention ?
        build_strategy(
            param_set, domain, domain.loc_ids[valid_locs_mask]
        ) :
        nothing
    return pref, decision_mat, strategy
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

# """
#     _update_state(cluster_ids::Vector, num_locs::Int64, max_members::Int64)

# Update list of:
# - selected clusters
# - number of locations selected for each cluster
# - the indices of clusters which exceed the membership rule
# - IDs of the above
# - and potential suitable alternative locations

# # Arguments
# - `cluster_ids` : ID of clusters
# - `num_locs` : Number of locations to select
# - `max_members` : Maximum number of members per cluster
# """
# function _update_state(cluster_ids::Vector, num_locs::Int64, max_members::Int64)
#     selected_clusters = cluster_ids[1:min(num_locs, length(cluster_ids))]
#     cluster_frequency = countmap(selected_clusters)
#     rule_violators_idx = findall(values(cluster_frequency) .> max_members)
#     exceeded_clusters = collect(keys(cluster_frequency))[rule_violators_idx]
#     potential_alternatives = cluster_ids[cluster_ids .∉ Ref(exceeded_clusters)]

#     return selected_clusters,
#     cluster_frequency,
#     rule_violators_idx,
#     exceeded_clusters,
#     potential_alternatives
# end
