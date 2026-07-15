module ADRIAanalysisRulesExt

using ADRIAanalysis
using ADRIAanalysis: Rule
using ADRIA: ADRIA
import ADRIA: is_test_env, model_spec
using SIRUS: SIRUS
using MLJ: MLJ
using StableRNGs: StableRNG
using DataFrames: DataFrame

# ---------------------------------------------------------------------------
# Private helpers (extension-internal only)
# ---------------------------------------------------------------------------

function _strip_value(condition_subclause::Vector)
    return join(condition_subclause[1:2], "__")
end

function _remove_duplicates(
    rules::T
)::T where {T<:Vector{Rule{Vector{Vector},Vector{Float64}}}}
    subclauses = join.([_strip_value.(r.condition) for r in rules], "_&_")
    unique_subclauses = unique(subclauses)

    n_unique_rules = length(unique_subclauses)
    if n_unique_rules == length(rules)
        return rules
    end

    n_rules = length(rules)
    n_duplicates = n_rules - n_unique_rules

    if !is_test_env()
        @warn "$n_duplicates of $n_rules duplicated rules were found and will be removed."
    end

    unique_rules::Vector{Rule} = Vector{Rule}(undef, n_unique_rules)
    for (unique_idx, unique_subclause) in enumerate(unique_subclauses)
        duplicate_rules_filter = unique_subclause .== subclauses
        if sum(duplicate_rules_filter) == 1
            unique_rules[unique_idx] = rules[duplicate_rules_filter][1]
            continue
        end
        duplicate_rules = rules[duplicate_rules_filter]
        max_probability_idx = findmax([r.consequent[1] for r in duplicate_rules])[2]
        unique_rules[unique_idx] = duplicate_rules[max_probability_idx]
    end

    return unique_rules
end

function _condition(sirus_rules::SIRUS.StableRules{Int64}, index::Int64)::Vector{Vector}
    condition::Vector{Vector} = []
    for subclause in sirus_rules.rules[index].clause.subclauses
        feature_name::String = subclause.feature_name
        direction::Symbol = subclause.direction
        value::Float32 = subclause.splitval
        push!(condition, [feature_name, direction, value])
    end
    return condition
end

function _consequent(sirus_rules::SIRUS.StableRules{Int64}, index::Int64)::Vector{Float64}
    weight = sirus_rules.weights[index]
    rule = sirus_rules.rules[index]
    # SIRUS._simplify_binary_probabilities is a private API (leading underscore).
    # Compat range covers "1.3, 2, 3" -- guard against removal across major versions.
    if !isdefined(SIRUS, :_simplify_binary_probabilities)
        @warn "SIRUS._simplify_binary_probabilities not found; falling back to raw probabilities. " *
            "Update the SIRUS compat bound and check the new public API."
        return Float64[rule.then, rule.otherwise]
    end
    then_prob = SIRUS._simplify_binary_probabilities(weight, rule.then)
    otherwise_prob = SIRUS._simplify_binary_probabilities(weight, rule.otherwise)
    return [then_prob, otherwise_prob]
end

# ---------------------------------------------------------------------------
# Public API -- adds concrete methods to the stubs declared in core
# ---------------------------------------------------------------------------

function ADRIAanalysis.rules(
    sirus_rules::SIRUS.StableRules{Int64}
)::Vector{Rule{Vector{Vector},Vector{Float64}}}
    return [
        Rule(_condition(sirus_rules, i), _consequent(sirus_rules, i))
        for i in eachindex(sirus_rules.rules)
    ]
end

function ADRIAanalysis.cluster_rules(
    result_set::ADRIA.ResultSet,
    clusters::Vector{T},
    scenarios::DataFrame,
    factors::Vector{Symbol},
    max_rules::T;
    seed::Int64=123,
    remove_duplicates::Bool=true,
    kwargs...
)::Vector{Rule{Vector{Vector},Vector{Float64}}} where {T<:Int64}
    ms = model_spec(result_set)

    # Build a set of non-constant fieldnames from the model spec.
    # Using a Set + filter (rather than logical indexing back into `factors`) avoids a
    # BoundsError when one or more elements of `factors` are not present in ms.fieldname
    # (e.g. DHW-derived columns added by feature_set, or `guided`).
    non_constant_in_spec = Set(ms[.!ms.is_constant, :fieldname])
    variable_factors::Vector{Symbol} = filter(f -> f in non_constant_in_spec, factors)

    if isempty(variable_factors)
        throw(ArgumentError("Factors of interest cannot be constant"))
    end

    X = scenarios[:, variable_factors]
    rng = StableRNG(seed)
    model = SIRUS.StableRulesClassifier(; max_rules=max_rules, rng=rng, kwargs...)
    mach = MLJ.machine(model, X, clusters)

    try
        MLJ.fit!(mach)
    catch err
        if !(err isa MethodError)
            rethrow(err)
        end
        error("Failed fitting SIRUS. Try increasing the number of scenarios/samples.")
    end

    if remove_duplicates
        return _remove_duplicates(ADRIAanalysis.rules(mach.fitresult))
    else
        return ADRIAanalysis.rules(mach.fitresult)
    end
end

function ADRIAanalysis.cluster_rules(
    result_set::ADRIA.ResultSet,
    clusters::Union{BitVector,Vector{Bool}},
    scenarios::DataFrame,
    factors::Vector{Symbol},
    max_rules::T;
    seed::Int64=123,
    remove_duplicates::Bool=true,
    kwargs...
)::Vector{Rule{Vector{Vector},Vector{Float64}}} where {T<:Int64}
    return ADRIAanalysis.cluster_rules(
        result_set, convert.(Int64, clusters), scenarios, factors, max_rules;
        seed=seed, remove_duplicates=remove_duplicates, kwargs...
    )
end

end  # module ADRIAanalysisRulesExt
