Base.@kwdef struct Intervention <: EcoModel
    # Intervention Factors
    # Bounds are defined as floats to maintain type stability
    guided::Param = Factor(
        0;
        ptype="ordered categorical",
        dist=CategoricalDistribution,
        dist_params=(-1.0, 0.0, 1.0),
        name="Guided",
        description="Intervention effort level: -1 no intervention (counterfactual), 0 unguided, 1 guided (MCDA-driven)."
    )
    mcda_method::Param = Factor(
        1;
        ptype="unordered categorical",
        dist=CategoricalDistribution,
        dist_params=(Tuple(1:Float64(length(decision.mcda_methods())))),
        name="MCDA Method",
        description="Which multi-criteria decision analysis method to use for location selection (only active when `guided` > 0)."
    )
    N_CAq_TA::Param = Factor(
        0;
        ptype="ordered discrete",
        dist=DiscreteOrderedUniformDist,
        dist_params=(0.0, 1_000_000.0, 50_000.0),  # increase in steps of 50K
        name="Seeded Tabular Acropora",
        description="Number of Tabular Acropora to seed per deployment event."
    )
    N_CAq_CA::Param = Factor(
        0;
        ptype="ordered discrete",
        dist=DiscreteOrderedUniformDist,
        dist_params=(0.0, 1_000_000.0, 50_000.0),  # increase in steps of 50K
        name="Seeded Corymbose Acropora",
        description="Number of Corymbose Acropora to seed per deployment event."
    )
    N_CAq_CNA::Param = Factor(
        0;
        ptype="ordered discrete",
        dist=DiscreteOrderedUniformDist,
        dist_params=(0.0, 1_000_000.0, 50_000.0),  # increase in steps of 50K
        name="Seeded Corymbose non-Acropora",
        description="Number of Corymbose non-Acropora to seed per deployment event."
    )
    N_CAq_SM::Param = Factor(
        0;
        ptype="ordered discrete",
        dist=DiscreteOrderedUniformDist,
        dist_params=(0.0, 1_000_000.0, 50_000.0),  # increase in steps of 50K
        name="Seeded Small Massives",
        description="Number of small massives/encrusting to seed per deployment event."
    )
    N_CAq_LM::Param = Factor(
        0;
        ptype="ordered discrete",
        dist=DiscreteOrderedUniformDist,
        dist_params=(0.0, 1_000_000.0, 50_000.0),  # increase in steps of 50K
        name="Seeded Large Massives",
        description="Number of large massives/encrusting to seed per deployment event."
    )
    N_LvM_settlers::Param = Factor(
        0;
        ptype="ordered discrete",
        dist=DiscreteOrderedUniformDist,
        dist_params=(0.0, 25_000_000.0, 1_000_000.0),  # increase in steps of 50K
        name="Moving corals settlers",
        description="Number of moving coral settlers added per deployment event."
    )
    CAq_devices_per_m2::Param = Factor(
        0;
        ptype="unordered categorical",
        dist=CategoricalDistribution,
        dist_params=(3, 5, 6, 9),
        name="seeding_devices_per_m2",
        description="Number of seeding devices per m²."
    )
    min_iv_locations::Param = Factor(
        5;
        ptype="ordered discrete",
        dist=DiscreteUniform,
        dist_params=(5.0, 100.0),
        name="Minimum intervention locations",
        description="Minimum number of locations to perform intervention"
    )
    LvM_min_iv_locations::Param = Factor(
        5;
        ptype="ordered discrete",
        dist=DiscreteUniform,
        dist_params=(1.0, 400.0),
        name="Moving corals minimum intervention locations",
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
    Shd::Param = Factor(
        0.0;
        ptype="continuous",
        dist=TriangularDist,
        dist_params=(0.0, 10.0, 0.0),
        name="SRM",
        description="Reduction in DHWs due to shading."
    )
    CAq_a_adapt::Param = Factor(
        0.0;
        ptype="ordered discrete",
        dist=DiscreteOrderedUniformDist,
        dist_params=(0.0, 15.0, 0.5),  # increase in steps of 0.5 DHW enhancement
        name="Assisted Adaptation",
        description="Assisted adaptation in terms of DHW resistance."
    )
    CAq_a_adapt_ref::Param = Factor(
        5.0;        # If CAq_a_adapt_ref == 0 uses first year as c_mean reference for entire run
        ptype="ordered discrete",
        dist=DiscreteOrderedUniformDist,
        dist_params=(0.0, 15.0, 1.0),
        name="Assisted adaptation reference",
        description="Distance from current year used as reference for assisted adaptation."
    )
    CAq_years::Param = Factor(
        10;
        ptype="ordered categorical",
        dist=DiscreteUniform,
        dist_params=(5.0, 75.0),
        name="Years to Seed",
        description="Number of years to seed for."
    )
    Shd_years::Param = Factor(
        10;
        ptype="ordered categorical",
        dist=DiscreteUniform,
        dist_params=(5.0, 75.0),
        name="Years to Shade",
        description="Number of years to shade for."
    )
    Fog_years::Param = Factor(
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
    projection_confidence::Param = Factor(
        0.0;
        ptype="continuous",
        dist=Uniform,
        dist_params=(0.0, 1.0),
        name="Projection Confidence",
        description="Confidence in far-future environmental projections when weighting them for site selection (only active when `guided` > 0)."
    )
    CAq_deployment_freq::Param = Factor(
        5;
        ptype="ordered categorical",
        dist=DiscreteUniform,
        dist_params=(0.0, 15.0),
        name="Selection Frequency (Seed)",
        description="Frequency of seeding deployments (0 deploys once)."
    )
    CAq_revisit_cadence::Param = Factor(
        0;
        ptype="ordered categorical",
        dist=DiscreteUniform,
        dist_params=(0.0, 15.0),
        name="Revisit Cadence (Seed)",
        description="Minimum number of years before a location can be re-selected for seed deployment (0 = no restriction)."
    )
    Fog_deployment_freq::Param = Factor(
        5;
        ptype="ordered categorical",
        dist=DiscreteUniform,
        dist_params=(0.0, 15.0),
        name="Selection Frequency (Fog)",
        description="Frequency of fogging deployments (0 deploys once)."
    )
    Fog_revisit_cadence::Param = Factor(
        0;
        ptype="ordered categorical",
        dist=DiscreteUniform,
        dist_params=(0.0, 15.0),
        name="Revisit Cadence (Fog)",
        description="Minimum number of years before a location can be re-selected for fog deployment (0 = no restriction)."
    )
    Shd_deployment_freq::Param = Factor(
        1;
        ptype="ordered categorical",
        dist=DiscreteUniform,
        dist_params=(0.0, 15.0),
        name="Deployment Frequency (Shading)",
        description="Frequency of shading deployments."
    )
    LvM_deployment_freq::Param = Factor(
        1;
        ptype="ordered categorical",
        dist=DiscreteUniform,
        dist_params=(0.0, 15.0),
        name="Deployment Frequency (Moving corals)",
        description="Frequency of moving corals deployments."
    )
    LvM_revisit_cadence::Param = Factor(
        0;
        ptype="ordered categorical",
        dist=DiscreteUniform,
        dist_params=(0.0, 15.0),
        name="Revisit Cadence (Moving corals)",
        description="Minimum number of years before a location can be re-selected for mc deployment (0 = no restriction)."
    )
    CAq_year_start::Param = Factor(
        2;
        ptype="ordered categorical",
        dist=DiscreteUniform,
        dist_params=(0.0, 25.0),
        name="Seeding Start Year",
        description="Start seeding deployments after this number of years has elapsed."
    )
    Shd_year_start::Param = Factor(
        2;
        ptype="ordered categorical",
        dist=DiscreteUniform,
        dist_params=(2.0, 25.0),
        name="Shading Start Year",
        description="Start of shading deployments after this number of years has elapsed."
    )
    Fog_year_start::Param = Factor(
        2;
        ptype="ordered categorical",
        dist=DiscreteUniform,
        dist_params=(2.0, 25.0),
        name="Fogging Start Year",
        description="Start of fogging deployments after this number of years has elapsed."
    )
    LvM_year_start::Param = Factor(
        2;
        ptype="ordered categorical",
        dist=DiscreteUniform,
        dist_params=(0.0, 25.0),
        name="Moving corals Start Year",
        description="Start moving corals deployments after this number of years has elapsed."
    )
    LvM_years::Param = Factor(
        10;
        ptype="ordered categorical",
        dist=DiscreteUniform,
        dist_params=(5.0, 75.0),
        name="Years to deploy moving corals",
        description="Number of years to deploy moving corals."
    )
    MCB_albedo::Param = Factor(
        0.0;
        ptype="ordered categorical",
        dist=CategoricalDistribution,
        dist_params=(0.0, 0.0),
        name="MCB Albedo",
        description="Albedo level to use from 5D DHW dataset."
    )
    MCB_duration::Param = Factor(
        0.0;
        ptype="ordered categorical",
        dist=CategoricalDistribution,
        dist_params=(0.0, 0.0),
        name="MCB Duration",
        description="Duration level (yearly days) to use from 5D DHW dataset."
    )
    MCB_deployment_freq::Param = Factor(
        1;
        ptype="ordered discrete",
        dist=DiscreteUniform,
        dist_params=(1.0, 50.0),
        name="MCB Deployment Frequency",
        description="Frequency of MCB deployment in years (e.g. 1 is every year, 2 is every 2nd year)."
    )

    # Intervention strategy parameters
    CAq_strategy::Param = Factor(
        2;
        ptype="ordered categorical",
        dist=CategoricalDistribution,
        dist_params=(1.0, 2.0),
        name="Seed Strategy Type",
        description="Deployment strategy: 1=Periodic (time-based), 2=Reactive (condition-based); 0 is off"
    )
    Fog_strategy::Param = Factor(
        2;
        ptype="ordered categorical",
        dist=CategoricalDistribution,
        dist_params=(1.0, 2.0),
        name="Fog Strategy Type",
        description="Deployment strategy: 1=Periodic (time-based), 2=Reactive (condition-based); 0 is off"
    )
    LvM_strategy::Param = Factor(
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
        dist_params=(0.2, 0.95, 0.05),
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
end

function interventions()
    return [:caq, :fog, :lvm]
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

    pref = preference(domain, param_set)

    # No locations meet the criteria (e.g. depth bounds exclude all candidates for this
    # scenario/year) - skip deployment rather than build an empty decision matrix
    if !any(valid_locs_mask)
        return pref, nothing, nothing
    end

    # Calculate cluster diversity and geographic separation scores
    diversity_scores = decision.cluster_diversity(domain.loc_data.cluster_id)
    separation_scores = decision.geographic_separation(domain.loc_data.mean_to_neighbor)

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
