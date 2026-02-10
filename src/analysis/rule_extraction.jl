using SIRUS, MLJ
using StableRNGs: StableRNG
import ..performance: temporal_variability

"""
    Rule{V<:Vector{Vector},W<:Vector{Float64}}

A Rule contains a condition (vector of conditional clauses) and consequent (vector
of values for then and else results of the condition). Conditional clauses represent
'a < b' statements. Consequents contain the probabilities of a observation begin
past of the 'positive class' given that the condition is true or false.
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
    rules(rules::SIRUS.StableRules{Int64})

Collates vector of Rule objects. These are **not** the same as SIRUS.Rule object.
See also [`Rule`](@ref).

# Arguments
- `rules` : SIRUS.StableRules object containing all rules

# Returns
Vector{ADRIA.analysis.Rule{Vector{Float64}, Vector{Vector}}}
"""
function rules(
    rules::SIRUS.StableRules{Int64}
)::Vector{Rule{Vector{Vector},Vector{Float64}}}
    return [
        Rule(_condition(rules, i), _consequent(rules, i)) for i in eachindex(rules.rules)
    ]
end

"""
    _condition(rules::SIRUS.StableRules{Int64}, index::Int64)

Vector containing condition clauses. Each condition clause is a vector with three
components: a feature_name::String, a direction::String (< or ≤) and a value:<Float64

# Arguments
- `rules` : SIRUS.StableRules object containing all rules
- `index` : Index of the rule

# Returns
Vector of Rule condition clauses (each one being a vector itself).
"""
function _condition(rules::SIRUS.StableRules{Int64}, index::Int64)::Vector{Vector}
    condition::Vector{Vector} = []

    for subclause in rules.rules[index].clause.subclauses
        feature_name::String = subclause.feature_name
        direction::Symbol = subclause.direction
        value::Float32 = subclause.splitval
        push!(condition, [feature_name, direction, value])
    end
    return condition
end

"""
    _consequent(rules::SIRUS.StableRules{Int64}, index::Int64)

Vector of Rule consequent with two components: the probability of the 'then' clause and
probability of the 'else' clause (here called otherwise)

# Arguments
- rules : SIRUS.StableRules object containing all rules
- `index` : Index of the rule

# Returns
Probabilities vector, one for Rule condition == true, one for Rule condition == false.
"""
function _consequent(rules::SIRUS.StableRules{Int64}, index::Int64)::Vector{Float64}
    weight = rules.weights[index]
    rule = rules.rules[index]
    then_probability = SIRUS._simplify_binary_probabilities(weight, rule.then)
    otherwise_probability = SIRUS._simplify_binary_probabilities(weight, rule.otherwise)
    return [then_probability, otherwise_probability]
end

"""
    print_rules(rules::Vector{Rule})

Print all rules in a human readable way. SIRUS also has this functionality but 1) it only
prints once, when the model is fitted and 2) this function prints the rules in a more
compact form.

# Arguments
- `rules` : Vector of Rule objects
"""
function print_rules(rules::Vector{Rule{Vector{Vector},Vector{Float64}}})::Nothing
    condition(rule) = [
        c[2] == :L ? "$(c[1]) < $(c[3])" : "$(c[1]) ≥ $(c[3])" for c in rule.condition
    ]
    consequent(rule) = " then $(rule.consequent[1]) else $(rule.consequent[2])\n"
    rule_string(rule) = "if " * join(condition(rule), " & ") * consequent(rule)
    print(join([rule_string(rule) for rule in rules]))

    return nothing
end

"""
    cluster_rules(result_set::ResultSet, clusters::Vector{T}, scenarios::DataFrame, factors::Vector{Symbol}, max_rules::T; seed::Int64=123, remove_duplicates::Bool=true, kwargs...)::Vector{Rule{Vector{Vector},Vector{Float64}}} where {T<:Int64}
    cluster_rules(result_set::ResultSet, clusters::Union{BitVector,Vector{Bool}}, scenarios::DataFrame, factors::Vector{Symbol}, max_rules::T; seed::Int64=123, remove_duplicates::Bool=true, kwargs...)::Vector{Rule{Vector{Vector},Vector{Float64}}} where {T<:Int64}

Use SIRUS package to extract rules from time series clusters based on some summary metric
(default is median). More information about the keyword arguments accepeted can be found in
MLJ's doc (https://juliaai.github.io/MLJ.jl/dev/models/StableRulesClassifier_SIRUS/).

# Arguments
- `result_set` : ResultSet.
- `clusters` : Vector of cluster indexes for each scenario outcome.
- `scenarios` : Scenarios DataFrame.
- `factors` : Vector of factors of interest.
- `max_rules` : Maximum number of rules, to be used as input by SIRUS.
- `seed` : Seed to be used by RGN. Defaults to 123.
- `remove_duplicates` : If true, duplicate rules will be removed from resulting ruleset. In
that case, the rule with the highest probability score is kept. Defaults to true.
- `kwargs` : Keyword arguments to be passed to `StableRulesClassifier`.

# Returns
A StableRules object (implemented by SIRUS).

# References
1. Huijzer, R.
   rikhuijzer/SIRUS.jl [Computer software].
   https://doi.org/10.5281/zenodo.7875577

2. Bénard, C., Biau, G., Da Veiga, S., & Scornet, E. 2021.
   Sirus: Stable and interpretable rule set for classification.
   Electron. J. Statist. 15 (1) 427 - 505.
   https://doi.org//10.1214/20-EJS1792
"""
function cluster_rules(
    result_set::ResultSet,
    clusters::Vector{T},
    scenarios::DataFrame,
    factors::Vector{Symbol},
    max_rules::T;
    seed::Int64=123,
    remove_duplicates::Bool=true,
    kwargs...
)::Vector{Rule{Vector{Vector},Vector{Float64}}} where {T<:Int64}
    ms = ADRIA.model_spec(result_set)

    variable_factors_filter::BitVector = .!ms[ms.fieldname .∈ [factors], :is_constant]
    variable_factors::Vector{Symbol} = factors[variable_factors_filter]

    if isempty(variable_factors)
        throw(ArgumentError("Factors of interest cannot be constant"))
    end

    X = scenarios[:, variable_factors]

    # Set seed and Random Number Generator
    rng = StableRNG(seed)

    # Use SIRUS Stable Rules Classifier model to extract the rules
    model = StableRulesClassifier(; max_rules=max_rules, rng=rng, kwargs...)
    mach = machine(model, X, clusters)

    try
        MLJ.fit!(mach)
    catch err
        if !(err isa MethodError)
            rethrow(err)
        end

        error("Failed fitting SIRUS. Try increasing the number of scenarios/samples.")
    end

    if remove_duplicates
        return _remove_duplicates(rules(mach.fitresult))
    else
        return rules(mach.fitresult)
    end
end
function cluster_rules(
    result_set::ResultSet,
    clusters::Union{BitVector,Vector{Bool}},
    scenarios::DataFrame,
    factors::Vector{Symbol},
    max_rules::T;
    seed::Int64=123,
    remove_duplicates::Bool=true,
    kwargs...
)::Vector{Rule{Vector{Vector},Vector{Float64}}} where {T<:Int64}
    return cluster_rules(
        result_set, convert.(Int64, clusters), scenarios, factors, max_rules;
        seed=seed, remove_duplicates=remove_duplicates, kwargs...
    )
end

"""
    _remove_duplicates(rules)::Vector{Rule{Vector{Vector},Vector{Float64}}}

Identifies and removes duplicate rulesets (if any are found).

The criteria to choose which rule to keep is based on the rule consequence probability (the one with the highest
probability is kept). If there are more than one rule with the same highest probability,
then the first one is chosen.

# Returns
A ruleset with duplicate rules removed
"""
function _remove_duplicates(
    rules::T
)::T where {T<:Vector{Rule{Vector{Vector},Vector{Float64}}}}
    # Extract subclauses from each rule without value
    subclauses = join.([_strip_value.(r.condition) for r in rules], "_&_")
    unique_subclauses = unique(subclauses)

    # Check if there are duplicate rules before moving on
    n_unique_rules = length(unique_subclauses)
    if n_unique_rules == length(rules)
        return rules
    end

    n_rules = length(rules)
    n_duplicates = n_rules - n_unique_rules

    if !is_test_env()
        @warn "$n_duplicates of $n_rules duplicated rules were found and are going to be removed."
    end

    unique_rules::Vector{Rule} = Vector{Rule}(undef, n_unique_rules)
    for (unique_idx, unique_subclause) in enumerate(unique_subclauses)
        duplicate_rules_filter = unique_subclause .== subclauses

        # If current rule has no duplicates go to next iteration
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

"""
    _strip_value(condition_subclause::Vector)

Helper function that extracts factor name and direction from a rule condition subclause.
Besides having just one line, this was extracted to a separate function to allow/facilitate
broadcasting this operation.
"""
function _strip_value(condition_subclause::Vector)
    return join(condition_subclause[1:2], "__")
end

"""
    maximum_probability(rules::SIRUS.StableRules{Int64})

Sum of biggest probabilities for each rule consequent

# Arguments
- `rules` : Vector of Rule objects
"""
function maximum_probability(rules::Vector{Rule{Vector{Vector},Vector{Float64}}})
    return sum([maximum(rule.consequent) for rule in rules])
end
