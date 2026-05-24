"""
    Rule{V<:Vector{Vector},W<:Vector{Float64}}

A Rule contains a condition (vector of conditional clauses) and a consequent (vector
of probability values for then/else outcomes). Conditional clauses represent
threshold comparisons on model factors.
"""
struct Rule{V<:Vector{Vector},W<:Vector{Float64}}
    condition::V
    consequent::W

    function Rule(condition, consequent)
        return new{typeof(condition),typeof(consequent)}(
            sort(condition; by=x -> x[1]), consequent
        )
    end
end

"""
    print_rules(rules::Vector{Rule})

Print all rules in a human-readable form. Does not require SIRUS.
"""
function print_rules(rules::Vector{Rule{Vector{Vector},Vector{Float64}}})::Nothing
    condition(rule) = [
        c[2] == :L ? "$(c[1]) < $(c[3])" : "$(c[1]) >= $(c[3])" for c in rule.condition
    ]
    consequent(rule) = " then $(rule.consequent[1]) else $(rule.consequent[2])\n"
    rule_string(rule) = "if " * join(condition(rule), " & ") * consequent(rule)
    print(join([rule_string(rule) for rule in rules]))

    return nothing
end

# Stubs -- concrete methods are added by ADRIAanalysisRulesExt when SIRUS+MLJ are loaded.
# Error stubs give users an actionable message rather than a confusing MethodError.
"""
    rules(sirus_rules)

Convert SIRUS.StableRules to a Vector{Rule}. Requires `using SIRUS, MLJ`.
"""
function rules(args...; kwargs...)
    return error(
        "rules requires SIRUS and MLJ. Load them first:\n\n" *
        "    using SIRUS, MLJ\n\n" *
        "then call rules again."
    )
end

"""
    cluster_rules(result_set, clusters, scenarios, factors, max_rules; kwargs...)

Fit a SIRUS classifier and extract interpretable rules. Requires `using SIRUS, MLJ`.
"""
function cluster_rules(args...; kwargs...)
    return error(
        "cluster_rules requires SIRUS and MLJ. Load them first:\n\n" *
        "    using SIRUS, MLJ\n\n" *
        "then call cluster_rules again."
    )
end

"""
    maximum_probability(rules::Vector{Rule{Vector{Vector},Vector{Float64}}})

Sum of the maximum consequent probability for each rule. Does not require SIRUS.
"""
function maximum_probability(rules::Vector{Rule{Vector{Vector},Vector{Float64}}})
    return sum([maximum(rule.consequent) for rule in rules])
end
