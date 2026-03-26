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
        dist_params=(0.0, 1_000_000.0, 50_000.0),  # increase in steps of 50K
        name="Seeded Tabular Acropora",
        description="Number of Tabular Acropora to seed per deployment event."
    )
    N_seed_CA::Param = Factor(
        0;
        ptype="ordered discrete",
        dist=DiscreteOrderedUniformDist,
        dist_params=(0.0, 1_000_000.0, 50_000.0),  # increase in steps of 50K
        name="Seeded Corymbose Acropora",
        description="Number of Corymbose Acropora to seed per deployment event."
    )
    N_seed_CNA::Param = Factor(
        0;
        ptype="ordered discrete",
        dist=DiscreteOrderedUniformDist,
        dist_params=(0.0, 1_000_000.0, 50_000.0),  # increase in steps of 50K
        name="Seeded Corymbose non-Acropora",
        description="Number of Corymbose non-Acropora to seed per deployment event."
    )
    N_seed_SM::Param = Factor(
        0;
        ptype="ordered discrete",
        dist=DiscreteOrderedUniformDist,
        dist_params=(0.0, 1_000_000.0, 50_000.0),  # increase in steps of 50K
        name="Seeded Small Massives",
        description="Number of small massives/encrusting to seed per deployment event."
    )
    N_seed_LM::Param = Factor(
        0;
        ptype="ordered discrete",
        dist=DiscreteOrderedUniformDist,
        dist_params=(0.0, 1_000_000.0, 50_000.0),  # increase in steps of 50K
        name="Seeded Large Massives",
        description="Number of large massives/encrusting to seed per deployment event."
    )
    N_mc_settlers::Param = Factor(
        0;
        ptype="ordered discrete",
        dist=DiscreteOrderedUniformDist,
        dist_params=(0.0, 25_000_000.0, 1_000_000.0),  # increase in steps of 50K
        name="Moving corals settlers",
        description="Number of moving coral settlers added per deployment event."
    )
    seeding_devices_per_m2::Param = Factor(
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
        dist_params=(5.0, 20.0),
        name="Minimum intervention locations",
        description="Minimum number of locations to perform intervention"
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
    a_adapt_ref::Param = Factor(
        0.0;        # If a_adapt_ref == 0 uses first year as c_mean reference for entire run
        ptype="ordered discrete",
        dist=DiscreteOrderedUniformDist,
        dist_params=(0.0, 15.0, 1.0),
        name="Assisted adaptation reference",
        description="Distance from current year used as referece for assisted adaptation."
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
        dist_params=(0.0, 15.0),
        name="Deployment Frequency (Shading)",
        description="Frequency of shading deployments."
    )
    mc_deployment_freq::Param = Factor(
        1;
        ptype="ordered categorical",
        dist=DiscreteUniform,
        dist_params=(0.0, 15.0),
        name="Deployment Frequency (Moving corals)",
        description="Frequency of moving corals deployments."
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
    return [:seed, :fog, :mc]
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
    diversity_scores = decision.cluster_diversity(domain.loc_data.cluster_id)
    separation_scores = decision.geographic_separation(domain.loc_data.mean_to_neighbor)

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
