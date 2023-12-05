Base.@kwdef struct CriteriaWeights{P,P2,P3,N} <: EcoModel
    seed_wave_stress::P = Param(
        1.0;
        ptype="real",
        bounds=(0.0, 1.0),
        default_bounds=(0.0, 1.0),
        dists="unif",
        criteria_keywords=(:seed, :weight),
        name="Seed Wave Stress",
        description="Importance of avoiding wave stress when seeding. Higher values places more weight on areas with low wave stress.",
    )
    seed_heat_stress::P = Param(
        1.0;
        ptype="real",
        bounds=(0.0, 1.0),
        default_bounds=(0.0, 1.0),
        dists="unif",
        criteria_keywords=(:seed, :weight),
        name="Seed Heat Stress",
        description="Importance of avoiding heat stress when seeding. Higher values places more weight on areas with low heat stress.",
    )
    fog_wave_stress::P = Param(
        1.0;
        ptype="real",
        bounds=(0.0, 1.0),
        default_bounds=(0.0, 1.0),
        dists="unif",
        criteria_keywords=(:fog, :weight),
        name="Shade Wave Stress",
        description="Importance of avoiding wave stress when fogging. Higher values places more weight on areas with low wave stress.",
    )
    fog_heat_stress::P = Param(
        1.0;
        ptype="real",
        bounds=(0.0, 1.0),
        default_bounds=(0.0, 1.0),
        dists="unif",
        criteria_keywords=(:fog, :weight),
        name="Shade Heat Stress",
        description="Importance of avoiding heat stress when fogging. Higher values places more weight on areas with low heat stress.",
    )
    fog_connectivity::P = Param(
        0.0;
        ptype="real",
        bounds=(0.0, 1.0),
        default_bounds=(0.0, 1.0),
        dists="unif",
        criteria_keywords=(:fog, :weight),
        name="Shade Connectivity",
        description="Higher values give preference to locations with high connectivity for shading deployments.",
    )
    seed_in_connectivity::P = Param(
        1.0;
        ptype="real",
        bounds=(0.0, 1.0),
        default_bounds=(0.0, 1.0),
        dists="unif",
        criteria_keywords=(:seed, :weight),
        name="Incoming Connectivity (Seed)",
        description="Higher values give preference to locations with high incoming connectivity (i.e., receives larvae from other sites) for enhanced coral deployments.",
    )
    seed_out_connectivity::P = Param(
        1.0;
        ptype="real",
        bounds=(0.0, 1.0),
        default_bounds=(0.0, 1.0),
        dists="unif",
        criteria_keywords=(:seed, :weight),
        name="Outgoing Connectivity (Seed)",
        description="Higher values give preference to locations with high outgoing connectivity (i.e., provides larvae to other sites) for enhanced coral deployments.",
    )
    seed_depth::P = Param(
        1.0;
        ptype="real",
        bounds=(0.0, 1.0),
        default_bounds=(0.0, 1.0),
        dists="unif",
        criteria_keywords=(:seed, :weight),
        name="Depth (Seed)",
        description="Higher values give preference to depper locations for enhanced coral deployments.",
    )
    coral_cover_low::P = Param(
        0.0;
        ptype="real",
        bounds=(0.0, 1.0),
        default_bounds=(0.0, 1.0),
        dists="unif",
        criteria_keywords=(:seed, :weight),
        name="Low Coral Cover",
        description="Higher values give greater preference to sites with low coral cover for seeding deployments.",
    )
    coral_cover_high::P = Param(
        0.0;
        ptype="real",
        bounds=(0.0, 1.0),
        default_bounds=(0.0, 1.0),
        dists="unif",
        criteria_keywords=(:fog, :weight),
        name="High Coral Cover",
        description="Higher values give preference to sites with high coral cover for shading deployments.",
    )
    seed_priority::P = Param(
        1.0;
        ptype="real",
        bounds=(0.0, 1.0),
        default_bounds=(0.0, 1.0),
        dists="unif",
        criteria_keywords=(:seed, :weight),
        name="Predecessor Priority (Seed)",
        description="Importance of seeding sites that provide larvae to priority reefs.",
    )
    fog_priority::P = Param(
        0.0;
        ptype="real",
        bounds=(0.0, 1.0),
        default_bounds=(0.0, 1.0),
        dists="unif",
        criteria_keywords=(:fog, :weight),
        name="Predecessor Priority (Shade)",
        description="Importance of shading sites that provide larvae to priority reefs.",
    )
    seed_zone::P = Param(
        0.0;
        ptype="real",
        bounds=(0.0, 1.0),
        default_bounds=(0.0, 1.0),
        dists="unif",
        criteria_keywords=(:seed, :weight),
        name="Zone Predecessor (Seed)",
        description="Importance of seeding sites that provide larvae to priority (target) zones.",
    )
    fog_zone::P = Param(
        0.0;
        ptype="real",
        bounds=(0.0, 1.0),
        default_bounds=(0.0, 1.0),
        dists="unif",
        criteria_keywords=(:fog, :weight),
        name="Zone Predecessor (Shade)",
        description="Importance of shading sites that provide larvae to priority (target) zones.",
    )
    deployed_coral_risk_tol::P2 = Param(
        1.0;
        ptype="real",
        bounds=(0.75, 1.0),
        default_bounds=(0.0, 1.0),
        dists="unif",
        criteria_keywords=(:threshold, :seed, :fog),
        name="Risk Tolerance",
        description="Filters out sites with heat/wave stress above threshold.",
    )
    depth_min::P3 = Param(
        5.0;
        ptype="real",
        bounds=(3.0, 5.0),
        default_bounds=(3.0, 5.0),
        dists="unif",
        criteria_keywords=(:threshold,),
        name="Minimum Depth",
        description="Minimum depth for a site to be included for consideration.\nNote: This value will be replaced with the shallowest depth value found if all sites are found to be deeper than `depth_min + depth_offset`.",
    )
    depth_offset::P3 = Param(
        10.0;
        ptype="real",
        bounds=(10.0, 25.0),
        default_bounds=(10.0, 25.0),
        dists="unif",
        criteria_keywords=(:threshold,),
        name="Depth Offset",
        description="Offset from minimum depth, used to indicate maximum depth.",
    )
end

function criteria_params(
    crit::DataFrame,
    criteria_keywords::Union{Tuple{Symbol,Symbol},Tuple{Symbol,Symbol,Symbol}},
)
    sel_crit = vec(
        all(
            hcat([in.(c_k, crit.criteria_keywords) for c_k in criteria_keywords]...); dims=2
        );
    )
    return crit[sel_crit, :]
end
