const Fog_CRITERIA = (
    Criterion(:heat_stress; val=1.0),
    Criterion(:wave_stress; val=1.0),
    Criterion(:in_connectivity; val=1.0),
    Criterion(:out_connectivity; val=1.0),
    Criterion(:depth; val=1.0, direction=minimum),
    Criterion(:coral_cover; val=0.0, direction=maximum),
    Criterion(:cluster_diversity; val=0.5),
    Criterion(:geographic_separation; val=0.5)
)

_define_criteria_weights!(@__MODULE__, :Fog, "Fogging", Fog_CRITERIA)
