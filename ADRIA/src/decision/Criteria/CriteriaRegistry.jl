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

const CRITERIA_SETS = Dict{DataType,@NamedTuple{prefix::String, criteria::Vector{Criterion}}}()

"""
    _factor_expr(spec::NamedTuple, label::AbstractString)::Expr

Expression `Factor(val; ...)` for the weight of a criterion. The factor name carries the
intervention (`"<name> (<label>)"`) and the description states the direction in use.

# Arguments
- `spec` : Spec of a [`Criterion`](@ref), including its catalog `name` and `description`.
- `label` : Human-readable name of the intervention, e.g. "Coral Aquaculture".

# Returns
`Expr` that evaluates to the `Param` of the weight.
"""
function _factor_expr(spec::NamedTuple, label::AbstractString)::Expr
    spec = merge(
        spec,
        (
            name="$(spec.name) ($label)",
            description="$(spec.description) Prefers $(spec.direction === maximum ? "higher" : "lower") values."
        )
    )
    kw = [Expr(:kw, k, v) for (k, v) in pairs(spec) if k != :val]
    return Expr(:call, :Factor, Expr(:parameters, kw...), spec.val)
end

"""
    _define_criteria_weights!(mod::Module, intervention::Symbol, label::AbstractString, criteria)

Define the `DecisionWeights` struct `<intervention>CriteriaWeights` in `mod` (e.g. `:CAq` gives
`CAqCriteriaWeights`), with one `Param` field `iv_<intervention>_<criterion>` per element of
`criteria` (in that order), and register the criteria of the struct for
[`criterion_names`](@ref) and [`criterion_fields`](@ref). `label` names the intervention in the
label of each factor, e.g. "Heat Stress (Coral Aquaculture)".

# Arguments
- `mod` : Module to define the struct in.
- `intervention` : Short intervention name, e.g. `:CAq`. Sets the struct name
    (`CAqCriteriaWeights`) and the factor prefix (`iv_CAq`).
- `label` : Human-readable intervention name, e.g. "Coral Aquaculture".
- `criteria` : [`Criterion`](@ref)s of the intervention, in field order.

# Returns
The defined struct type.

# Throws
`ArgumentError` if `criteria` contains the same criterion twice.
"""
function _define_criteria_weights!(
    mod::Module, intervention::Symbol, label::AbstractString, criteria
)::DataType
    type_name = Symbol(intervention, :CriteriaWeights)
    prefix = "iv_$(intervention)"
    doc = "    $type_name <: DecisionWeights

Criteria weights for $(lowercase(label)) interventions."
    crit = collect(Criterion, criteria)
    allunique(c.name for c in crit) || throw(ArgumentError("Duplicate criteria in `$type_name`"))

    fields = [
        Expr(
            :(=),
            Expr(:(::), Symbol(prefix, "_", c.name), :Param),
            _factor_expr(c.spec, label)
        )
        for c in crit
    ]
    Core.eval(mod, quote
        Base.@kwdef struct $(type_name) <: DecisionWeights
            $(fields...)
        end
    end)

    T = Base.invokelatest(getfield, mod, type_name)
    CRITERIA_SETS[T] = (prefix=String(prefix), criteria=crit)
    Core.eval(mod, :(Core.@doc $doc $type_name))

    return T
end

"""
    criterion_names(T::Type{<:DecisionWeights})::Vector{Symbol}

Bare criterion names of `T` (e.g. `:heat_stress`), in field order.

# Arguments
- `T` : Weights type defined with [`_define_criteria_weights!`](@ref), e.g. `CAqCriteriaWeights`.

# Returns
Criterion names, matching the names used in the decision matrix.
"""
criterion_names(T::Type{<:DecisionWeights})::Vector{Symbol} =
    Symbol[c.name for c in CRITERIA_SETS[T].criteria]

"""
    criterion_fields(T::Type{<:DecisionWeights})::Vector{Symbol}

Factor names of the weights in `T` (e.g. `:iv_CAq_heat_stress`), in field order.

# Arguments
- `T` : Weights type defined with [`_define_criteria_weights!`](@ref), e.g. `CAqCriteriaWeights`.

# Returns
Factor names, matching the `fieldname` column of the model spec.
"""
criterion_fields(T::Type{<:DecisionWeights})::Vector{Symbol} =
    Symbol[Symbol(CRITERIA_SETS[T].prefix, "_", c.name) for c in CRITERIA_SETS[T].criteria]

"""
    DecisionPreferences(W::Type{<:DecisionWeights}, dom)
    DecisionPreferences(W::Type{<:DecisionWeights}, dom, params::YAXArray)

Decision preferences of the intervention with criteria weights `W`, using the default weights
of the model spec of `dom`, or those of one scenario.

# Arguments
- `W` : Weights type defined with [`_define_criteria_weights!`](@ref), e.g. `CAqCriteriaWeights`.
- `dom` : Domain, which holds the model spec.
- `params` : Parameters of one scenario, with a `factors` dimension.

# Returns
`DecisionPreferences` with bare criterion names (e.g. `:heat_stress`), as used in the decision
matrix.
"""
function DecisionPreferences(W::Type{<:DecisionWeights}, dom, params::YAXArray)
    w::DataFrame = component_params(dom, W)
    return DecisionPreferences(
        criterion_names(W), params[factors = At(string.(w.fieldname))], w.direction
    )
end
function DecisionPreferences(W::Type{<:DecisionWeights}, dom)
    w::DataFrame = component_params(dom, W)
    return DecisionPreferences(criterion_names(W), w.val, w.direction)
end
