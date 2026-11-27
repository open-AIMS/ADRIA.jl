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
