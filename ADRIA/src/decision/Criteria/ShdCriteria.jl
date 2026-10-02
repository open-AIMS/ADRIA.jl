"""
    ShdCriteriaWeights <: DecisionWeights

Weights for shading interventions.
"""
Base.@kwdef struct ShdCriteriaWeights <: DecisionWeights
    iv_Shd_heat_stress::Param = Factor(
        1.0;
        ptype="continuous",
        dist=Uniform,
        dist_params=(0.0, 1.0),
        direction=minimum,
        name="Shading Heat Stress",
        description="Preference locations with lower heat stress for shading."
    )
    iv_Shd_wave_stress::Param = Factor(
        1.0;
        ptype="continuous",
        dist=Uniform,
        dist_params=(0.0, 1.0),
        direction=minimum,
        name="Shading Wave Stress",
        description="Prefer locations with lower wave stress for shading."
    )
    iv_Shd_connectivity::Param = Factor(
        0.0;
        ptype="continuous",
        dist=Uniform,
        dist_params=(0.0, 1.0),
        direction=maximum,
        name="Shading Connectivity",
        description="Preference locations with higher outgoing connectivity for shading."
    )
    iv_Shd_coral_cover::Param = Factor(
        0.0;
        ptype="continuous",
        dist=Uniform,
        dist_params=(0.0, 1.0),
        direction=maximum,
        name="Coral Cover (Shading)",
        description="Give greater weight to locations with higher coral cover for shading."
    )
    iv_Shd_priority::Param = Factor(
        0.0;
        ptype="continuous",
        dist=Uniform,
        dist_params=(0.0, 1.0),
        direction=maximum,
        name="Predecessor Priority (Shading)",
        description="Relative importance of locations with higher outgoing connectivity to priority locations."
    )
    iv_Shd_zone::Param = Factor(
        0.0;
        ptype="continuous",
        dist=Uniform,
        dist_params=(0.0, 1.0),
        direction=maximum,
        name="Zone Predecessor (Shading)",
        description="Relative importance of locations with higher outgoing connectivitiy to priority (target) zones."
    )
end

function ShdPreferences(
    dom, params::YAXArray
)::DecisionPreferences
    w::DataFrame = component_params(dom.model, ShdCriteriaWeights)

    return DecisionPreferences(
        string.(w.fieldname), params[At(string.(w.fieldname))], w.direction
    )
end
