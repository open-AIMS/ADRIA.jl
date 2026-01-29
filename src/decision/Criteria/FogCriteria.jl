"""
    FogCriteriaWeights <: DecisionWeights

Criteria weights for fogging interventions.
"""
Base.@kwdef struct FogCriteriaWeights <: DecisionWeights
    fog_heat_stress::Param = Factor(
        1.0;
        ptype="continuous",
        dist=Uniform,
        dist_params=(0.0, 1.0),
        direction=minimum,
        name="Fog Heat Stress",
        description="Preference locations with lower heat stress for fogging."
    )
    fog_wave_stress::Param = Factor(
        1.0;
        ptype="continuous",
        dist=Uniform,
        dist_params=(0.0, 1.0),
        direction=minimum,
        name="Fog Wave Stress",
        description="Preference locations with lower wave activity for fogging."
    )
    fog_in_connectivity::Param = Factor(
        1.0;
        ptype="continuous",
        dist=Uniform,
        dist_params=(0.0, 1.0),
        direction=maximum,
        name="Incoming Connectivity (Fog)",
        description="Give preference to locations with high incoming connectivity (i.e., receives larvae from other sites) for fogging deployments."
    )
    fog_out_connectivity::Param = Factor(
        1.0;
        ptype="continuous",
        dist=Uniform,
        dist_params=(0.0, 1.0),
        direction=maximum,
        name="Outgoing Connectivity (Fog)",
        description="Give preference to locations with high outgoing connectivity (i.e., provides larvae to other sites) for fogging deployments."
    )
    fog_depth::Param = Factor(
        1.0;
        ptype="continuous",
        dist=Uniform,
        dist_params=(0.0, 1.0),
        direction=minimum,
        name="Depth (Fog)",
        description="Give preference to shallower locations for fogging deployments."
    )
    fog_coral_cover::Param = Factor(
        0.0;
        ptype="continuous",
        dist=Uniform,
        dist_params=(0.0, 1.0),
        direction=maximum,
        name="Fog Coral Cover",
        description="Higher values give preference to sites with high coral cover for fogging deployments."
    )
    fog_cluster_diversity::Param = Factor(
        0.5;
        ptype="continuous",
        dist=Uniform,
        dist_params=(0.0, 1.0),
        direction=maximum,
        name="Cluster Diversity",
        description="Prefer to fog locations in clusters that are under-represented."
    )
    fog_geographic_separation::Param = Factor(
        0.5;
        ptype="continuous",
        dist=Uniform,
        dist_params=(0.0, 1.0),
        direction=minimum,
        name="Geographic Separation",
        description="Fog locations that are distant (when maximized) or closer (when minimized) to their neighbors."
    )
    # Disabled as they are currently unnecessary
    # fog_priority::Param= Factor(
    #     0.0;
    #     ptype="continuous",
    #     dist=Uniform,
    #     dist_params=(0.0, 1.0),
    #     direction=maximum,
    #     name="Predecessor Priority (Fog)",
    #     description="Importance of fogging sites that provide larvae to priority reefs.",
    # )
    # fog_zone::Param= Factor(
    #     0.0;
    #     ptype="continuous",
    #     dist=Uniform,
    #     dist_params=(0.0, 1.0),
    #     direction=maximum,
    #     name="Zone Predecessor (Fog)",
    #     description="Importance of fogging sites that provide larvae to priority (target) zones.",
    # )
end

# Alias default constructor
FogPreferences(names, criteria, directions) = DecisionPreferences(
    names, criteria, directions
)

function FogPreferences(dom, params::YAXArray)::DecisionPreferences
    w::DataFrame = component_params(dom.model, FogCriteriaWeights)
    cn = Symbol[Symbol(join(split(string(cn), "_")[2:end], "_")) for cn in w.fieldname]

    return DecisionPreferences(cn, params[factors = At(string.(w.fieldname))], w.direction)
end
