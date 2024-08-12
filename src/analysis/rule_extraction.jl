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
    condition(rule) =
        [c[2] == :L ? "$(c[1]) < $(c[3])" : "$(c[1]) ≥ $(c[3])" for c in rule.condition]
    consequent(rule) = " then $(rule.consequent[1]) else $(rule.consequent[2])\n"
    rule_string(rule) = "if " * join(condition(rule), " & ") * consequent(rule)
    print(join([rule_string(rule) for rule in rules]))

    return nothing
end

"""
    cluster_rules(clusters::Vector{T}, X::DataFrame, max_rules::T; seed::Int64=123, kwargs...) where {T<:Integer,F<:Real}
    cluster_rules(clusters::Union{BitVector,Vector{Bool}}, X::DataFrame, max_rules::T; kwargs...) where {T<:Int64}

Use SIRUS package to extract rules from time series clusters based on some summary metric
(default is median). More information about the keyword arguments accepeted can be found in
MLJ's doc (https://juliaai.github.io/MLJ.jl/dev/models/StableRulesClassifier_SIRUS/).

# Arguments
- `clusters` : Vector of cluster indexes for each scenario outcome
- `X` : Features to be used as input by SIRUS
- `max_rules` : Maximum number of rules, to be used as input by SIRUS
- `seed` : Seed to be used by RGN

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
function cluster_rules(clusters::Vector{T}, X::DataFrame, max_rules::T;
    seed::Int64=123, kwargs...) where {T<:Int64}
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

    return rules(mach.fitresult)
end
function cluster_rules(clusters::Union{BitVector,Vector{Bool}}, X::DataFrame, max_rules::T;
    kwargs...) where {T<:Int64}
    return cluster_rules(convert.(Int64, clusters), X, max_rules; kwargs...)
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
