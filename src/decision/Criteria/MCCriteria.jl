"""
    MCCriteriaWeights <: DecisionWeights

Criteria weights for moving corals (MC) intervention.
"""
Base.@kwdef struct MCCriteriaWeights <: DecisionWeights
    mc_heat_stress::Param = Factor(
        0.9;
        ptype="continuous",
        dist=Uniform,
        dist_params=(0.8, 1.0),
        direction=minimum,
        name="MC Heat Stress",
        description="Importance of avoiding heat stress. Prefer locations with lower heat stress."
    )
    mc_wave_stress::Param = Factor(
        0.5;
        ptype="continuous",
        dist=Uniform,
        dist_params=(0.0, 1.0),
        direction=maximum,
        name="MC Wave Stress",
        description="Prefer locations with higher wave activity."
    )
    mc_in_connectivity::Param = Factor(
        0.5;
        ptype="continuous",
        dist=Uniform,
        dist_params=(0.5, 1.0),
        direction=maximum,
        name="MC Incoming Connectivity",
        description="Give preference to locations with high incoming connectivity (i.e., receives larvae from other sites) for coral deployments."
    )
    mc_out_connectivity::Param = Factor(
        0.80;
        ptype="continuous",
        dist=Uniform,
        dist_params=(0.5, 1.0),
        direction=maximum,
        name="MC Outgoing Connectivity",
        description="Give preference to locations with high outgoing connectivity (i.e., provides larvae to other sites) for coral deployments."
    )
    mc_depth::Param = Factor(
        1.0;
        ptype="continuous",
        dist=Uniform,
        dist_params=(0.8, 1.0),
        direction=maximum,
        name="MC Depth",
        description="Give preference to deeper locations for coral deployments."
    )
    mc_coral_cover::Param = Factor(
        0.7;
        ptype="continuous",
        dist=Uniform,
        dist_params=(0.0, 1.0),
        direction=minimum,
        name="MC Coral Cover",
        description="Preference locations with lower coral cover (higher available space)."
    )
    mc_cluster_diversity::Param = Factor(
        0.7;
        ptype="continuous",
        dist=Uniform,
        dist_params=(0.0, 1.0),
        direction=maximum,
        name="MC Cluster Diversity",
        description="Prefer locations from clusters that are under-represented."
    )
    mc_geographic_separation::Param = Factor(
        0.8;
        ptype="continuous",
        dist=Uniform,
        dist_params=(0.0, 1.0),
        direction=minimum,
        name="MC Geographic Separation",
        description="Prefer locations that are distant (when maximized) or closer (when minimized; the default) to their neighbors."
    )
end

"""
    MCPreferences <: DecisionPreference

Preference type specific for moving corals interventions to allow specific routines.
"""
struct MCPreferences <: DecisionPreference
    names::Vector{Symbol}
    weights::Vector{Float64}
    directions::Vector{Function}
end

function MCPreferences(dom, params::YAXArray)::MCPreferences
    w::DataFrame = component_params(dom.model, MCCriteriaWeights)
    cn = Symbol[Symbol(join(split(string(cn), "_")[2:end], "_")) for cn in w.fieldname]

    return MCPreferences(cn, params[factors=At(string.(w.fieldname))], w.direction)
end
function MCPreferences(dom, params...)::MCPreferences
    w::DataFrame = component_params(dom.model, MCCriteriaWeights)
    for (k, v) in params
        w[w.fieldname .== k, :val] .= v
    end

    return MCPreferences(w.fieldname, w.val, w.direction)
end
function MCPreferences(dom)
    w::DataFrame = component_params(dom.model, MCCriteriaWeights)
    return MCPreferences(w.fieldname, w.val, w.direction)
end
