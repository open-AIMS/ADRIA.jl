"""
    SRMCriteriaWeights <: DecisionWeights

Weights for shading (Solar Radiation Management) interventions.
"""
Base.@kwdef struct SRMCriteriaWeights <: DecisionWeights
    srm_heat_stress::Param= Factor(
        1.0;
        ptype="continuous",
        dist=Uniform,
        dist_params=(0.0, 1.0),
        direction=minimum,
        name="Shade Heat Stress",
        description="Preference locations with lower heat stress for SRM.",
    )
    srm_wave_stress::Param= Factor(
        1.0;
        ptype="continuous",
        dist=Uniform,
        dist_params=(0.0, 1.0),
        direction=minimum,
        name="Shade Wave Stress",
        description="Prefer locations with lower wave stress for SRM.",
    )
    srm_connectivity::Param= Factor(
        0.0;
        ptype="continuous",
        dist=Uniform,
        dist_params=(0.0, 1.0),
        direction=maximum,
        name="SRM Connectivity",
        description="Preference locations with higher outgoing connectivity for SRM.",
    )
    srm_coral_cover::Param= Factor(
        0.0;
        ptype="continuous",
        dist=Uniform,
        dist_params=(0.0, 1.0),
        direction=maximum,
        name="Coral Cover (SRM)",
        description="Give greater weight to locations with higher coral cover for SRM.",
    )
    srm_priority::Param= Factor(
        0.0;
        ptype="continuous",
        dist=Uniform,
        dist_params=(0.0, 1.0),
        direction=maximum,
        name="Predecessor Priority (SRM)",
        description="Relative importance of locations with higher outgoing connectivity to priority locations.",
    )
    srm_zone::Param= Factor(
        0.0;
        ptype="continuous",
        dist=Uniform,
        dist_params=(0.0, 1.0),
        direction=maximum,
        name="Zone Predecessor (SRM)",
        description="Relative importance of locations with higher outgoing connectivitiy to priority (target) zones.",
    )
end

function SRMPreferences(
    dom, params::YAXArray
)::DecisionPreferences
    w::DataFrame = component_params(dom.model, SRMCriteriaWeights)

    return DecisionPreferences(string.(w.fieldname), params[At(string.(w.fieldname))], w.direction)
end
