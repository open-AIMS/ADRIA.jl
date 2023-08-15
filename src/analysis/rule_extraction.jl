using SIRUS, MLJ
using StableRNGs: StableRNG
import ..performance: temporal_variability

"""
    Rule{V<:Vector{Vector},W<:Vector{Float64}}

A Rule contains a condition (vector of conditional clauses) and consequent (vector
of values for then and else results of the condition)
"""
struct Rule{V<:Vector{Vector},W<:Vector{Float64}}
    condition::V
    consequent::W
end

"""
    rules(rules::SIRUS.StableRules{Float64})

Collates vector of Rule objects. These are **not** the same as SIRUS.Rule object.

# Arguments
- `rules` : SIRUS.StableRules object

# Returns
Vector{ADRIA.analysis.Rule{Vector{Float64}, Vector{Vector}}}
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

# Returns
Vector{Vector}
Vector of vectors of rule conditions.
"""
function _condition(rules::SIRUS.StableRules{Float64}, index::Int64)
    condition::Vector{Vector} = []
    for split in rules.rules[index].path.splits
        feature_name = split.splitpoint.feature_name
        direction = split.direction
        value = split.splitpoint.value
        push!(condition, [feature_name, direction, value])
    end
    return condition
end

"""
    _consequent(rules::SIRUS.StableRules{Float64}, index::Int64)

Vector of Rule consequent with two components: the probability of the 'then' clause and
probability of the 'else' clause (here called otherwise)

# Arguments
- rules : SIRUS.StableRules object containing all rules
- `index` : Index of the rule

# Returns
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
function print_rules(rules::Vector{Rule{Vector{Vector},Vector{Float64}}})
    condition(rule) = [c == :L ? "$(c[1]) < $(c[3])" : "$(c[1]) ≥ $(c[3])" for c in rule.condition]
    consequent(rule) = " then $(rule.consequent[1]) else $(rule.consequent[2])\n"
    rule_string(rule) = "if " * join(condition(rule), " ") * consequent(rule)
    print(join([rule_string(rule) for rule in rules]))
end


"""
    cluster_rules(clusters::Vector{T}, X::DataFrame, outcomes::AbstractMatrix{F}, max_rules::Int64; n_trees::Int64=1000) where {T<:Integer,F<:Real}

Use SIRUS package to extract rules from time series clusters based on some summary metric (default is median)

# Arguments
- `clusters` - Vector of cluster indexes for each scenario outcome
- `X::DataFrame` - DataFrame with factors to be used as input by SIRUS
- `outcomes` - Matrix of some metric over time for all scenarios
- `max_rules` - Maximum number of rules, to be used as input by SIRUS
- `n_trees` - Number of trees to be created by SIRUS algorithm

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
function cluster_rules(clusters::Vector{T}, X::DataFrame, outcomes::AbstractMatrix{F},
    max_rules::T; n_trees::T=1000) where {T<:Int64,F<:Real}
    # Find robust cluster index
    rc = robust_cluster(clusters, outcomes)
    y = clusters .== rc

    # Set seed and Random Number Generator
    seed = 123
    _rng(seed::Int=1) = StableRNG(seed)

    # Use SIRUS Stable Rules Classifier model to extract the rules
    model = StableRulesClassifier(; max_rules=max_rules, n_trees=n_trees)
    mach = machine(model, X, y)
    fit!(mach)
    return rules(mach.fitresult)
end

"""
    robust_cluster(clusters::Vector{T}, result_set::ResultSet) where {T<: Integer}

Find most robust cluster.

# Arguments
- `clusters` - Vector with outcomes cluster indexes
- `outcomes` - AbstractMatrix of outcomes factors

# Returns
Index of the most robust cluster
"""
function robust_cluster(clusters::Vector{T}, outcomes::AbstractMatrix{F}) where {T<:Int64,F<:Real}
    clusters_statistics::Vector{Float64} = []

    max_outcome = maximum(outcomes)
    for cluster in unique(clusters)
        normalized_outcomed = outcomes[:, clusters.==cluster] ./ max_outcome
        statistic = median(temporal_variability(normalized_outcomed; w=[0.9, 0.1]))
        push!(clusters_statistics, statistic)
    end

    # Find index of the cluster with maximum value of summary statistic
    return argmax(clusters_statistics)
end

"""
    maximum_probability(rules::SIRUS.StableRules{Float64})

Sum of biggest probabilities for each rule consequent

# Arguments
- `rules` : Vector of Rule objects
"""
function maximum_probability(rules::Vector{Rule{Vector{Vector},Vector{Float64}}})
    sum([maximum(rule.consequent) for rule in rules])
end
