const Shd_CRITERIA = (
    Criterion(:heat_stress; val=1.0),
    Criterion(:wave_stress; val=1.0),
    Criterion(:connectivity; val=0.0),
    Criterion(:coral_cover; val=0.0, direction=maximum),
    Criterion(:priority; val=0.0),
    Criterion(:zone; val=0.0)
)

_define_criteria_weights!(@__MODULE__, :Shd, "Shading", Shd_CRITERIA)

function ShdPreferences(
    dom, params::YAXArray
)::DecisionPreferences
    w::DataFrame = component_params(dom, ShdCriteriaWeights)

    return DecisionPreferences(
        string.(w.fieldname), params[At(string.(w.fieldname))], w.direction
    )
end
