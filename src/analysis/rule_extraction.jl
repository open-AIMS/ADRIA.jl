using SIRUS, MLJ
using StableRNGs: StableRNG
import ..performance: temporal_variability

"""
    Rule{V<:Vector{Float64},W<:Vector{Vector{Any}}}

A Rule contains a condition (vector of conditional clauses) and consequent (vector 
of values for then and else results of the condition)
"""
struct Rule{V<:Vector{Float64},W<:Vector{Vector{Any}}}
    condition::W
    consequent::V
end

"""
    rules(rules::SIRUS.StableRules)

Vector of Rule objects. These are **not** the same as SIRUS.Rule object.

# Arguments
- `rules` : SIRUS.StableRules object

# Return
Vector{ADRIA.analysis.Rule{Vector{Float64}, Vector{Vector{Any}}}}
"""
function rules(rules::SIRUS.StableRules{Float64})
    [Rule(_condition(rules, i), _consequent(rules, i)) for i in eachindex(rules.rules)]
end

"""
    _condition(rules::SIRUS.StableRules{Float64}, index::Int64)

Vector containing condition clauses. Each condition clause is a vector with three 
components: a feature_name::String, a direction::String (< or ≤) and a value:<Float64 

# Arguments
- rules : SIRUS.StableRules object containing all rules
- `index` : Index of the rule

# Return
Vector{Vector{Any}}
Vector of vectors of rule conditions.
"""
function _condition(rules::SIRUS.StableRules{Float64}, index::Int64)
    feature_name(split) = split.splitpoint.feature_name
    direction(split) = split.direction
    value(split) = split.splitpoint.value

    [[feature_name(s), direction(s), value(s)] for s in rules.rules[index].path.splits]
end

"""
    _consequent(rules::SIRUS.StableRules{Float64}, index::Int64)

Vector of Rule consequent with two components: the probability of the 'then' clause and 
probability of the 'else' clause (here called otherwise)

# Arguments
- rules : SIRUS.StableRules object containing all rules
- `index` : Index of the rule

# Return
Vector{Float64}
"""
function _consequent(rules::SIRUS.StableRules{Float64}, index::Int64)
    weight = rules.weights[index]
    rule = rules.rules[index]
    then_probability = SIRUS._simplify_binary_probabilities(weight, rule.then)
    otherwise_probability = SIRUS._simplify_binary_probabilities(weight, rule.otherwise)
    [then_probability, otherwise_probability]
end

"""
    print_rules(rules::Vector{Rule})

Print all rules in a human readable way. SIRUS also has this functionality but 1) it only 
prints once, when the model is fitted and 2) this function prints the rules in a more 
compact form.

# Arguments
- `rules` : Vector of Rule objects
"""
function print_rules(rules::Vector{Rule{Vector{Float64},Vector{Vector{Any}}}})
    condition(rule) = [c == :L ? "$(c[1]) < $(c[3])" : "$(c[1]) ≥ $(c[3])" for c in rule.condition]
    consequent(rule) = " then $(rule.consequent[1]) else $(rule.consequent[2])\n"
    rule_string(rule) = string("if ", join(condition(rule), " "), consequent(rule))
    print(join([rule_string(rule) for rule in rules]))
end
