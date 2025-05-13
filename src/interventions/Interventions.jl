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
    seeding_strategy::Param = Factor(
        0;
        ptype="unordered categorical",
        dist=CategoricalDistribution,
        dist_params=(Tuple(0:length(decision.seeding_strategies()))),
        name="Seeding Strategy",
        description="Seeding strategy to deploy corals."
    )
    seeding_density::Param = Factor(
        5.0;
        ptype="continuous",
        dist=TriangularDist,
        dist_params=(2.0, 12.0, 5.0),
        name="Seeding Density",
        description="Density to seed corals."
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
    N_seed_SM::Param = Factor(
        0;
        ptype="ordered discrete",
        dist=DiscreteOrderedUniformDist,
        dist_params=(0.0, 1000000.0, 50000.0),  # increase in steps of 50K
        name="Seeded Small Massives",
        description="Number of small massives/encrusting to seed per deployment event."
    )
    min_iv_locations::Param = Factor(
        5;
        ptype="ordered discrete",
        dist=DiscreteUniform,
        dist_params=(5.0, 20.0),
        name="Intervention Locations",
        description="Minimum number of deployment locations (for both seeding and fogging)."
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
        dist=DiscreteTriangularDist,
        dist_params=(5.0, 75.0, 5.0),
        name="Years to Seed",
        description="Number of years to seed for."
    )
    shade_years::Param = Factor(
        10;
        ptype="ordered categorical",
        dist=DiscreteTriangularDist,
        dist_params=(5.0, 75.0, 5.0),
        name="Years to Shade",
        description="Number of years to shade for."
    )
    fog_years::Param = Factor(
        10;
        ptype="ordered categorical",
        dist=DiscreteTriangularDist,
        dist_params=(5.0, 75.0, 5.0),
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
end

function interventions()
    return [:seed, :fog]
end
