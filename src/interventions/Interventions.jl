Base.@kwdef struct Intervention{N,P,P2} <: EcoModel
    # Intervention Parameters
    # Integer values have a +1 offset to allow for discrete value mapping
    # (see `set()` and `map_to_discrete()` methods)
    # Bounds are defined as floats to maintain type stability
    guided::N = Param(
        0;
        ptype="categorical",
        bounds=(-1.0, length(decision.mcda_methods()) + 1.0),
        default_bounds=(-1.0, length(decision.mcda_methods()) + 1.0),
        dists="unif",
        criteria_keywords=(""),
        name="Guided",
        description="Choice of MCDA approach.",
    )
    N_seed_TA::N = Param(
        0;
        ptype="integer",
        bounds=(0.0, 1000000.0 + 1.0),
        default_bounds=(0.0, 1000000.0 + 1.0),
        dists="unif",
        criteria_keywords=(""),
        name="Seeded Tabular Acropora",
        description="Number of Tabular Acropora to seed per deployment year.",
    )
    N_seed_CA::N = Param(
        0;
        ptype="integer",
        bounds=(0.0, 1000000.0 + 1.0),
        default_bounds=(0.0, 1000000.0 + 1.0),
        dists="unif",
        criteria_keywords=(""),
        name="Seeded Corymbose Acropora",
        description="Number of Corymbose Acropora to seed per deployment year.",
    )
    N_seed_SM::N = Param(
        0;
        ptype="integer",
        bounds=(0.0, 1000000.0 + 1.0),
        default_bounds=(0.0, 1000000.0 + 1.0),
        dists="unif",
        criteria_keywords=(""),
        name="Seeded Small Massives",
        description="Number of small massives/encrusting to seed per deployment year.",
    )
    fogging::P = Param(
        0.16;
        ptype="real",
        bounds=(0.0, 0.3, 0.16 / 0.3),
        default_bounds=(0.0, 0.3, 0.16 / 0.3),
        dists="triang",
        criteria_keywords=(""),
        name="Fogging",
        description="Assumed reduction in bleaching mortality.",
    )
    SRM::P = Param(
        0.0;
        ptype="real",
        bounds=(0.0, 7.0, 0.0),
        default_bounds=(0.0, 7.0, 0.0),
        dists="triang",
        criteria_keywords=(""),
        name="SRM",
        description="Reduction in DHWs due to shading.",
    )
    a_adapt::P = Param(
        0.0;
        ptype="real",
        bounds=(0.0, 8.0, 0.0),
        default_bounds=(0.0, 8.0, 0.0),
        dists="triang",
        criteria_keywords=(""),
        name="Assisted Adaptation",
        description="Assisted adaptation in terms of DHW resistance.",
    )
    seed_years::P2 = Param(
        10;
        ptype="integer",
        bounds=(5.0, 74.0 + 1.0, 5 / 70),
        default_bounds=(5.0, 74.0 + 1.0, 5 / 70),
        dists="triang",
        criteria_keywords=(""),
        name="Years to Seed",
        description="Number of years to seed for.",
    )
    shade_years::P2 = Param(
        10;
        ptype="integer",
        bounds=(5.0, 74.0 + 1.0, 5 / 70),
        default_bounds=(5.0, 74.0 + 1.0, 5 / 70),
        dists="triang",
        criteria_keywords=(""),
        name="Years to Shade",
        description="Number of years to shade for.",
    )
    fog_years::P2 = Param(
        10;
        ptype="integer",
        bounds=(5.0, 74.0 + 1.0, 5 / 70),
        default_bounds=(5.0, 74.0 + 1.0, 5 / 70),
        dists="triang",
        criteria_keywords=(""),
        name="Years to fog",
        description="Number of years to fog for.",
    )
    plan_horizon::N = Param(
        5;
        ptype="integer",
        bounds=(0.0, 40.0 + 1.0),
        default_bounds=(0.0, 40.0 + 1.0),
        dists="unif",
        criteria_keywords=(""),
        name="Planning Horizon",
        description="How many years of projected data to take into account when selecting intervention locations (0 only accounts for current year).",
    )
    seed_freq::N = Param(
        5;
        ptype="integer",
        bounds=(0.0, 15.0 + 1.0),
        default_bounds=(0.0, 15.0 + 1.0),
        dists="unif",
        criteria_keywords=(""),
        name="Seeding Frequency",
        description="Frequency of seeding site selection (0 is set and forget).",
    )
    shade_freq::N = Param(
        1;
        ptype="integer",
        bounds=(0.0, 15.0 + 1.0),
        default_bounds=(0.0, 15.0 + 1.0),
        dists="unif",
        criteria_keywords=(""),
        name="Shading Frequency",
        description="Frequency of shading site selection (0 is set and forget).",
    )
    fog_freq::N = Param(
        1;
        ptype="integer",
        bounds=(0.0, 15.0 + 1.0),
        default_bounds=(0.0, 15.0 + 1.0),
        dists="unif",
        criteria_keywords=(""),
        name="Fogging Frequency",
        description="Frequency of fogging site selection (0 is set and forget).",
    )
    seed_year_start::N = Param(
        2;
        ptype="integer",
        bounds=(2.0, 25.0 + 1.0),
        default_bounds=(2.0, 25.0 + 1.0),
        dists="unif",
        criteria_keywords=(""),
        name="Seeding Start Year",
        description="Start seeding deployments after this number of years has elapsed.",
    )
    shade_year_start::N = Param(
        2;
        ptype="integer",
        bounds=(2.0, 25.0 + 1.0),
        default_bounds=(2.0, 25.0 + 1.0),
        dists="unif",
        criteria_keywords=(""),
        name="Shading Start Year",
        description="Start of shading deployments after this number of years has elapsed.",
    )
    fog_year_start::N = Param(
        2;
        ptype="integer",
        bounds=(2.0, 25.0 + 1.0),
        default_bounds=(2.0, 25.0 + 1.0),
        dists="unif",
        criteria_keywords=(""),
        name="Fogginging Start Year",
        description="Start of fogging deployments after this number of years has elapsed.",
    )
end
