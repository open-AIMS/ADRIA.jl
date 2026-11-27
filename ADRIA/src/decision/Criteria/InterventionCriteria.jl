# Coral aquaculture
const CAq_CRITERIA = (
    Criterion(:heat_stress; val=0.9, dist_params=(0.8, 1.0)),
    Criterion(:wave_stress; val=0.5, direction=maximum),
    Criterion(:in_connectivity; val=0.5, dist_params=(0.2, 1.0)),
    Criterion(:out_connectivity; val=0.80, dist_params=(0.2, 1.0)),
    Criterion(:out_spread; val=0.80, dist_params=(0.2, 1.0)),
    Criterion(:depth; val=1.0, dist_params=(0.8, 1.0)),
    Criterion(:coral_cover; val=0.7),
    Criterion(:cluster_diversity; val=0.7),
    Criterion(:geographic_separation; val=0.8)
    # Disabled as they are currently unnecessary:
    # Criterion(:priority; val=1.0),
    # Criterion(:zone; val=0.0)
)

_define_criteria_weights!(@__MODULE__, :CAq, "Coral Aquaculture", CAq_CRITERIA)

# Fogging
const Fog_CRITERIA = (
    Criterion(:heat_stress; val=1.0),
    Criterion(:wave_stress; val=1.0),
    Criterion(:in_connectivity; val=1.0),
    Criterion(:out_connectivity; val=1.0),
    Criterion(:out_spread; val=1.0),
    Criterion(:depth; val=1.0, direction=minimum),
    Criterion(:coral_cover; val=0.0, direction=maximum),
    Criterion(:cluster_diversity; val=0.5),
    Criterion(:geographic_separation; val=0.5)
    # Disabled as they are currently unnecessary:
    # Criterion(:priority; val=0.0),
    # Criterion(:zone; val=0.0)
)

_define_criteria_weights!(@__MODULE__, :Fog, "Fogging", Fog_CRITERIA)

# Larval methods
const LvM_CRITERIA = (
    Criterion(:heat_stress; val=0.9, dist_params=(0.8, 1.0)),
    Criterion(:wave_stress; val=0.5, direction=maximum),
    Criterion(:in_connectivity; val=0.5, dist_params=(0.5, 1.0)),
    Criterion(:out_connectivity; val=0.80, dist_params=(0.5, 1.0)),
    Criterion(:out_spread; val=0.80, dist_params=(0.5, 1.0)),
    Criterion(:depth; val=1.0, dist_params=(0.8, 1.0)),
    Criterion(:coral_cover; val=0.7),
    Criterion(:cluster_diversity; val=0.7),
    Criterion(:geographic_separation; val=0.8)
)

_define_criteria_weights!(@__MODULE__, :LvM, "Larval Methods", LvM_CRITERIA)
