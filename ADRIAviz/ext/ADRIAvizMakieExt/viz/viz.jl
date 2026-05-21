using
    GeoMakie,
    GeoMakie.Colors

const COLORMAP_TYPE = Union{Symbol,RGB,RGBA,Vector{Symbol},Vector{RGBA},Vector{RGB}}

include("scenarios.jl")
include("sensitivity.jl")
include("clustering.jl")
include("rule_extraction.jl")
include("location_selection.jl")
include("spatial.jl")
include("taxa_dynamics.jl")
include("environment/dhw.jl")
include("environment/cyclones.jl")
include("data_envelopment.jl")
