"""
    Criterion

Specification of one decision criterion of an intervention: its bare `name` (e.g.
`:heat_stress`) and the factor spec of its weight. Build one with
[`Criterion`](@ref)`(name; kwargs...)`.
"""
struct Criterion
    name::Symbol
    spec::NamedTuple
end

"""
Definition of each known criterion: the label and description it has in every intervention,
and its default weight (`val`), range (`dist_params`) and `direction`. An intervention only
overrides `val`, `dist_params` and `direction`.

`direction` is the function that picks the preferred value: `maximum` prefers higher values
of the criterion, `minimum` prefers lower ones.
"""
const CRITERIA_CATALOG = Dict{Symbol,NamedTuple}(
    :heat_stress => (
        val=1.0,
        ptype="continuous",
        dist=Uniform,
        dist_params=(0.0, 1.0),
        direction=minimum,
        name="Heat Stress",
        description="Importance of heat stress when selecting locations."
    ),
    :wave_stress => (
        val=1.0,
        ptype="continuous",
        dist=Uniform,
        dist_params=(0.0, 1.0),
        direction=minimum,
        name="Wave Stress",
        description="Importance of wave stress when selecting locations."
    ),
    :in_connectivity => (
        val=1.0,
        ptype="continuous",
        dist=Uniform,
        dist_params=(0.0, 1.0),
        direction=maximum,
        name="Incoming Connectivity",
        description="Importance of incoming connectivity (i.e., receiving larvae from other locations)."
    ),
    :out_connectivity => (
        val=1.0,
        ptype="continuous",
        dist=Uniform,
        dist_params=(0.0, 1.0),
        direction=maximum,
        name="Outgoing Connectivity",
        description="Importance of outgoing connectivity (i.e., providing larvae to other locations)."
    ),
    :depth => (
        val=1.0,
        ptype="continuous",
        dist=Uniform,
        dist_params=(0.0, 1.0),
        direction=maximum,
        name="Depth",
        description="Importance of depth when selecting locations."
    ),
    :coral_cover => (
        val=1.0,
        ptype="continuous",
        dist=Uniform,
        dist_params=(0.0, 1.0),
        direction=minimum,
        name="Coral Cover",
        description="Importance of coral cover (lower cover means more available space)."
    ),
    :cluster_diversity => (
        val=1.0,
        ptype="continuous",
        dist=Uniform,
        dist_params=(0.0, 1.0),
        direction=maximum,
        name="Cluster Diversity",
        description="Importance of selecting locations from clusters that are under-represented."
    ),
    :geographic_separation => (
        val=1.0,
        ptype="continuous",
        dist=Uniform,
        dist_params=(0.0, 1.0),
        direction=minimum,
        name="Geographic Separation",
        description="Importance of the mean distance of a location to its neighbors."
    ),
    :connectivity => (
        val=1.0,
        ptype="continuous",
        dist=Uniform,
        dist_params=(0.0, 1.0),
        direction=maximum,
        name="Connectivity",
        description="Importance of outgoing connectivity (i.e., providing larvae to other locations)."
    ),
    :priority => (
        val=1.0,
        ptype="continuous",
        dist=Uniform,
        dist_params=(0.0, 1.0),
        direction=maximum,
        name="Predecessor Priority",
        description="Importance of locations that provide larvae to priority locations."
    ),
    :zone => (
        val=1.0,
        ptype="continuous",
        dist=Uniform,
        dist_params=(0.0, 1.0),
        direction=maximum,
        name="Zone Predecessor",
        description="Importance of locations that provide larvae to priority (target) zones."
    )
)

"""
    Criterion(name::Symbol; kwargs...)::Criterion

Criterion `name` as defined in `CRITERIA_CATALOG`, with `val`, `dist_params` or `direction`
(and any other [`Factor`](@ref) keyword except `name` and `description`) overridden by
`kwargs`.

# Arguments
- `name` : Criterion to build. Must be a key of `CRITERIA_CATALOG`.
- `kwargs` : Overrides of the catalog spec, e.g. `val`, `dist_params`, `direction`.
    `name` and `description` are not accepted, as they are set in the catalog.

# Returns
The `Criterion`, with its factor spec validated.

# Throws
`ArgumentError` if `name` is not in `CRITERIA_CATALOG`, or if `kwargs` contains `name` or
`description`.

# Example
```julia
Criterion(:heat_stress)  # catalog defaults
Criterion(:heat_stress; val=0.9, dist_params=(0.8, 1.0))
```
"""
function Criterion(name::Symbol; kwargs...)::Criterion
    haskey(CRITERIA_CATALOG, name) || throw(
        ArgumentError("Unknown criterion `$name`. Define it in `CRITERIA_CATALOG` first.")
    )
    overrides = NamedTuple(kwargs)
    for k in (:name, :description)
        haskey(overrides, k) && throw(
            ArgumentError("`$k` of `$name` is set in `CRITERIA_CATALOG`, not per intervention")
        )
    end

    spec = merge(CRITERIA_CATALOG[name], overrides)
    _factor_spec(spec.val; (k => v for (k, v) in pairs(spec) if k != :val)...)  # validate
    return Criterion(name, spec)
end
