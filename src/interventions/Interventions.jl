Base.@kwdef struct Intervention{DU,T,DT} <: EcoModel
    # Intervention Factors
    # Bounds are defined as floats to maintain type stability
    guided::DU = Factor(
        0;
        ptype="unordered categorical",
        dist=DiscreteUniform,
        dist_params=(-1.0, Float64(length(decision.mcda_methods()))),
        name="Guided",
        description="Choice of MCDA approach.",
    )
    N_seed_TA::DU = Factor(
        0;
        ptype="ordered categorical",
        dist=DiscreteUniform,
        dist_params=(0.0, 1000000.0),
        name="Seeded Tabular Acropora",
        description="Number of Tabular Acropora to seed per deployment year.",
    )
    N_seed_CA::DU = Factor(
        0;
        ptype="ordered categorical",
        dist=DiscreteUniform,
        dist_params=(0.0, 1000000.0),
        name="Seeded Corymbose Acropora",
        description="Number of Corymbose Acropora to seed per deployment year.",
    )
    N_seed_SM::DU = Factor(
        0;
        ptype="ordered categorical",
        dist=DiscreteUniform,
        dist_params=(0.0, 1000000.0),
        name="Seeded Small Massives",
        description="Number of small massives/encrusting to seed per deployment year.",
    )
    fogging::T = Factor(
        0.16;
        ptype="continuous",
        dist=TriangularDist,
        dist_params=(0.0, 0.3, 0.16),
        name="Fogging",
        description="Assumed reduction in bleaching mortality.",
    )
    SRM::T = Factor(
        0.0;
        ptype="continuous",
        dist=TriangularDist,
        dist_params=(0.0, 7.0, 0.0),
        name="SRM",
        description="Reduction in DHWs due to shading.",
    )
    a_adapt::T = Factor(
        0.0;
        ptype="continuous",
        dist=TriangularDist,
        dist_params=(0.0, 8.0, 0.0),
        name="Assisted Adaptation",
        description="Assisted adaptation in terms of DHW resistance.",
    )
    seed_years::DT = Factor(
        10;
        ptype="ordered categorical",
        dist=DiscreteTriangularDist,
        dist_params=(5.0, 74.0, 5.0),
        name="Years to Seed",
        description="Number of years to seed for.",
    )
    shade_years::DT = Factor(
        10;
        ptype="ordered categorical",
        dist=DiscreteTriangularDist,
        dist_params=(5.0, 74.0, 5.0),
        name="Years to Shade",
        description="Number of years to shade for.",
    )
    fog_years::DT = Factor(
        10;
        ptype="ordered categorical",
        dist=DiscreteTriangularDist,
        dist_params=(5.0, 74.0, 5.0),
        name="Years to fog",
        description="Number of years to fog for.",
    )
    plan_horizon::DU = Factor(
        5;
        ptype="ordered categorical",
        dist=DiscreteUniform,
        dist_params=(0.0, 40.0),
        name="Planning Horizon",
        description="How many years of projected data to take into account when selecting intervention locations (0 only accounts for current deployment year).",
    )
    seed_freq::DU = Factor(
        5;
        ptype="ordered categorical",
        dist=DiscreteUniform,
        dist_params=(0.0, 15.0),
        name="Seeding Frequency",
        description="Frequency of seeding site selection (0 is set and forget).",
    )
    shade_freq::DU = Factor(
        1;
        ptype="ordered categorical",
        dist=DiscreteUniform,
        dist_params=(0.0, 15.0),
        name="Shading Frequency",
        description="Frequency of shading site selection (0 is set and forget).",
    )
    fog_freq::DU = Factor(
        1;
        ptype="ordered categorical",
        dist=DiscreteUniform,
        dist_params=(0.0, 15.0),
        name="Fogging Frequency",
        description="Frequency of fogging site selection (0 is set and forget).",
    )
    seed_year_start::DU = Factor(
        2;
        ptype="ordered categorical",
        dist=DiscreteUniform,
        dist_params=(0.0, 25.0),
        name="Seeding Start Year",
        description="Start seeding deployments after this number of years has elapsed.",
    )
    shade_year_start::DU = Factor(
        2;
        ptype="ordered categorical",
        dist=DiscreteUniform,
        dist_params=(2.0, 25.0),
        name="Shading Start Year",
        description="Start of shading deployments after this number of years has elapsed.",
    )
    fog_year_start::DU = Factor(
        2;
        ptype="ordered categorical",
        dist=DiscreteUniform,
        dist_params=(2.0, 25.0),
        name="Fogging Start Year",
        description="Start of fogging deployments after this number of years has elapsed.",
    )
end
