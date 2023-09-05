"""
Deprecated collection of functions for use with Random Forests
(from the DecisionTree.jl package).

No longer used as RFs take too long to build in the context of the
app. Can't afford to wait 1-2mins every time the data slice is
updated.
"""

"""
    outcome_probability(data::AbstractVector)::NamedTuple

Determine probability occurrence.
"""
function outcome_probability(data::AbstractVector)::NamedTuple
    p_outcomes = cdf.(fit(Distributions.Normal, data), data)

    n = length(data)
    return (
        values=[
            count(p_outcomes .> 0.80) / n,
            count((p_outcomes .> 0.70) .& (p_outcomes .<= 0.80)) / n,
            count((p_outcomes .>= 0.50) .& (p_outcomes .<= 0.70)) / n,
            count((p_outcomes .> 0.20) .& (p_outcomes .< 0.50)) / n,
            count(p_outcomes .< 0.20) / n],
        labels=["Very High\n> 80%", "High\n70 - 80%", "Medium\n50 - 70%", "Low\n20 - 50%", "Very Low\n< 20%"]
    )
end

function probability_table(
    model, X, y;
    bstrap=BalancedSampling(100),
    c_intv=BCaConfInt(0.95),
    outcome_classes=["Very Low (< 25%)", "Low (25 - 55%)", "Medium (55 - 70%)", "High (70 - 85%)", "Very High (> 85%)"]
)
    outcome_order = unique(y)
    outcome_order_idx = indexin(outcome_classes, outcome_order)
    filter!(x -> !isnothing(x), outcome_order_idx)
    outcomes_display_order = outcome_order[outcome_order_idx]

    ŷ = DataFrame(apply_forest_proba(model, X, outcome_order), outcome_order)
    df = ŷ[:, outcomes_display_order]  # Reorder to low-high

    # Bootstrap confidence intervals
    conf_prob(x) = bootstrap(mean, x, bstrap) |> xb -> confint(xb, c_intv)

    proba_intv = DataFrame([conf_prob(col)[1] for col in eachcol(df)])
    rename!(proba_intv, [:mean, :lower, :upper])
    insertcols!(proba_intv, 1, "Outcome" => outcomes_display_order)

    return proba_intv
end

"""
    _accuracy(model, y, X)

Calculate predictive accuracy of random forest.
"""
function _accuracy(model, y, X)
    return DecisionTree.accuracy(y, apply_forest(model, X))
end


function _permutation_importance(
    trees::U,
    labels::AbstractVector{T},
    features::AbstractVecOrMat{S},
    score::Function,
    n_iter::Int=3;
    rng=Random.GLOBAL_RNG
) where {S,T,U<:Union{<:DecisionTree.Ensemble{S,T},<:DecisionTree.Root{S,T},<:DecisionTree.LeafOrNode{S,T},Tuple{<:DecisionTree.Ensemble{S,T},AbstractVector{Float64}}}}

    base = score(trees, labels, features)
    scores = Matrix{Float64}(undef, size(features, 2), n_iter)
    rng = DecisionTree.mk_rng(rng)::Random.AbstractRNG

    @floop for (i, col) in enumerate(eachcol(features))
        origin = copy(col)
        scores[i, :] = map(1:n_iter) do _
            shuffle!(rng, col)
            base - score(trees, labels, features)
        end

        features[:, i] = origin
    end

    (mean=reshape(mapslices(scores, dims=2) do im
            mean(im)
        end, :),
        std=reshape(mapslices(scores, dims=2) do im
            std(im)
        end, :),
        scores=scores)
end


"""
Extract feature importance from random forest.
"""
function ft_importance(model::Ensemble{Float64,Any}, X::DataFrame, p::Vector; rng::Int64=101)::DataFrame
    X_base = copy(X)
    insertcols!(X_base, 1, :dummy => zeros(nrow(X_base)))

    p1 = _permutation_importance(model, p, Matrix(X_base), _accuracy, 25; rng=rng)

    norm = replace(p1.mean ./ (maximum(p1.mean) - minimum(p1.mean)), Inf => 0.0, NaN => 0.0)

    feat_importance = DataFrame((param=names(X_base), mean=p1.mean, std=p1.std, norm=norm))
    sort!(feat_importance, :norm, rev=true)

    dummy_idx = findall(feat_importance.param .== "dummy")[1]
    if dummy_idx == 1
        imp = feat_importance[2:11, [:param, :norm]]
    else
        imp = feat_importance[1:min(10, dummy_idx - 1), [:param, :norm]]
    end

    return imp
end
