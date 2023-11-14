Base.@kwdef struct Criteria{P,P2,P3,N} <: EcoModel
    wave_stress::P = Param(
        1.0;
        ptype="real",
        bounds=(0.0, 1.0),
        default_bounds=(0.0, 1.0),
        dists="unif",
        criteria_keywords=("seed", "shade", "weight"),
        name="Wave Stress",
        description="Importance of avoiding wave stress. Higher values places more weight on areas with low wave stress.",
    )
    heat_stress::P = Param(
        1.0;
        ptype="real",
        bounds=(0.0, 1.0),
        default_bounds=(0.0, 1.0),
        dists="unif",
        criteria_keywords=("seed", "shade", "weight"),
        name="Heat Stress",
        description="Importance of avoiding heat stress. Higher values places more weight on areas with low heat stress.",
    )
    shade_connectivity::P2 = Param(
        0.0;
        ptype="real",
        bounds=(0.0, 1.0),
        default_bounds=(0.0, 1.0),
        dists="unif",
        criteria_keywords=("shade", "weight"),
        name="Shade Connectivity",
        description="Higher values give preference to locations with high connectivity for shading deployments.",
    )
    in_seed_connectivity::P2 = Param(
        1.0;
        ptype="real",
        bounds=(0.0, 1.0),
        default_bounds=(0.0, 1.0),
        dists="unif",
        criteria_keywords=("seed", "weight"),
        name="Incoming Connectivity (Seed)",
        description="Higher values give preference to locations with high incoming connectivity (i.e., receives larvae from other sites) for enhanced coral deployments.",
    )
    out_seed_connectivity::P2 = Param(
        1.0;
        ptype="real",
        bounds=(0.0, 1.0),
        default_bounds=(0.0, 1.0),
        dists="unif",
        criteria_keywords=("seed", "weight"),
        name="Outgoing Connectivity (Seed)",
        description="Higher values give preference to locations with high outgoing connectivity (i.e., provides larvae to other sites) for enhanced coral deployments.",
    )
    coral_cover_low::P2 = Param(
        0.0;
        ptype="real",
        bounds=(0.0, 1.0),
        default_bounds=(0.0, 1.0),
        dists="unif",
        criteria_keywords=("seed", "weight"),
        name="Low Coral Cover",
        description="Higher values give greater preference to sites with low coral cover for seeding deployments.",
    )
    coral_cover_high::P2 = Param(
        0.0;
        ptype="real",
        bounds=(0.0, 1.0),
        default_bounds=(0.0, 1.0),
        dists="unif",
        criteria_keywords=("shade", "weight"),
        name="High Coral Cover",
        description="Higher values give preference to sites with high coral cover for shading deployments.",
    )
    seed_priority::P2 = Param(
        1.0;
        ptype="real",
        bounds=(0.0, 1.0),
        default_bounds=(0.0, 1.0),
        dists="unif",
        criteria_keywords=("seed", "weight"),
        name="Predecessor Priority (Seed)",
        description="Importance of seeding sites that provide larvae to priority reefs.",
    )
    shade_priority::P2 = Param(
        0.0;
        ptype="real",
        bounds=(0.0, 1.0),
        default_bounds=(0.0, 1.0),
        dists="unif",
        criteria_keywords=("shade", "weight"),
        name="Predecessor Priority (Shade)",
        description="Importance of shading sites that provide larvae to priority reefs.",
    )
    zone_seed::P2 = Param(
        0.0;
        ptype="real",
        bounds=(0.0, 1.0),
        default_bounds=(0.0, 1.0),
        dists="unif",
        criteria_keywords=("seed", "weight"),
        name="Zone Predecessor (Seed)",
        description="Importance of seeding sites that provide larvae to priority (target) zones.",
    )
    zone_shade::P2 = Param(
        0.0;
        ptype="real",
        bounds=(0.0, 1.0),
        default_bounds=(0.0, 1.0),
        dists="unif",
        criteria_keywords=("shade", "weight"),
        name="Zone Predecessor (Shade)",
        description="Importance of shading sites that provide larvae to priority (target) zones.",
    )
    coral_cover_tol::P2 = Param(
        0.2;
        ptype="real",
        bounds=(0.0, 1.0),
        default_bounds=(0.0, 1.0),
        dists="unif",
        criteria_keywords=("threshold", "seed"),
        name="Low Area Tolerance",
        description="Tolerance for low proportional space for seeding deployments.",
    )
    deployed_coral_risk_tol::P = Param(
        1.0;
        ptype="real",
        bounds=(0.75, 1.0),
        default_bounds=(0.0, 1.0),
        dists="unif",
        criteria_keywords=("threshold", "seed", "shade"),
        name="Risk Tolerance",
        description="Filters out sites with heat/wave stress above threshold.",
    )
    use_dist::N = Param(
        1;
        ptype="categorical",
        bounds=(0.0, 1.0 + 1.0),
        default_bounds=(0.0, 1.0 + 1.0),
        dists="unif",
        criteria_keywords=("threshold",),
        name="Use Distance Threshold",
        description="Turns distance sorting on or off.",
    )
    dist_thresh::P3 = Param(
        0.1;
        ptype="real",
        bounds=(0.0, 1.0),
        default_bounds=(0.0, 1.0),
        dists="unif",
        criteria_keywords=("threshold",),
        name="Distance Threshold",
        description="Sites selected by MCDA must be further apart than median(dist)-dist_thresh*median(dist).",
    )
    depth_min::P3 = Param(
        5.0;
        ptype="real",
        bounds=(3.0, 5.0),
        default_bounds=(3.0, 5.0),
        dists="unif",
        criteria_keywords=("threshold",),
        name="Minimum Depth",
        description="Minimum depth for a site to be included for consideration.\nNote: This value will be replaced with the shallowest depth value found if all sites are found to be deeper than `depth_min + depth_offset`.",
    )
    depth_offset::P3 = Param(
        10.0;
        ptype="real",
        bounds=(10.0, 25.0),
        default_bounds=(10.0, 25.0),
        dists="unif",
        criteria_keywords=("threshold",),
        name="Depth Offset",
        description="Offset from minimum depth, used to indicate maximum depth.",
    )
end
