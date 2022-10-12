function outcome_probability(data)
    p_outcomes = cdf.(fit(Distributions.Normal, data), data)

    p = Vector{Any}(p_outcomes)
    p[p_outcomes.>0.70] .= "High (70 - 85%)"
    p[p_outcomes.>0.85] .= "Very High (> 85%)"
    p[p_outcomes.<0.55] .= "Low (15 - 55%)"
    p[p_outcomes.<0.15] .= "Very Low (< 15%)"
    p[(p_outcomes.>=0.55).&(p_outcomes.<=0.70)] .= "Medium (55 - 70%)"

    return p
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
    trees   :: U,
    labels  :: AbstractVector{T},
    features:: AbstractVecOrMat{S},
    score   :: Function,
    n_iter  :: Int = 3;
    rng     =  Random.GLOBAL_RNG
    ) where {S, T, U <: Union{<: DecisionTree.Ensemble{S, T}, <: DecisionTree.Root{S, T}, <: DecisionTree.LeafOrNode{S, T}, Tuple{<: DecisionTree.Ensemble{S, T}, AbstractVector{Float64}}}}

    base = score(trees, labels, features)
    scores = Matrix{Float64}(undef, size(features, 2), n_iter)
    rng = DecisionTree.mk_rng(rng)::Random.AbstractRNG

    ThreadsX.foreach(enumerate(eachcol(features))) do (i, col)
        origin = copy(col)
        scores[i, :] = map(1:n_iter) do _
            shuffle!(rng, col)
            base - score(trees, labels, features)
        end

        @inbounds features[:, i] = origin
    end

    # origin = similar(features[:, 1], Any)
    # non_constants = map(d -> !all(d .== d[1]), eachcol(features))
    # for (i, col) in enumerate(eachcol(features))
        # if non_constants[i] == 0
        #     scores[i, :] .= 0.0
        #     continue
        # end

    #     origin .= copy(col)
    #     scores[i, :] .= map(1:n_iter) do _
    #         shuffle!(rng, col)
    #         base - score(trees, labels, features)
    #     end

    #     features[:, i] .= origin
    # end

    # Main.@infiltrate

    (mean = reshape(mapslices(scores, dims = 2) do im
        mean(im)
    end, :),
    std = reshape(mapslices(scores, dims = 2) do im
        std(im)
    end, :),
    scores = scores)
end


"""
Extract feature importance from random forest.
"""
function ft_importance(model::Ensemble{Float64, Any}, X::DataFrame, p::Vector; rng::Int64=101)::DataFrame
    X_base = copy(X)
    insertcols!(X_base, 1, :dummy=>zeros(nrow(X_base)))

    # DecisionTree.accuracy(y, apply_forest(model, X))
    # apply_forest_proba(model, X, y)
    @time p1 = _permutation_importance(model, p, Matrix(X_base), _accuracy, 25; rng=rng)

    norm = replace(p1.mean ./ (maximum(p1.mean) - minimum(p1.mean)), Inf=>0.0, NaN=>0.0)

    feat_importance = DataFrame((param=names(X_base), mean=p1.mean, std=p1.std, norm=norm))
    sort!(feat_importance, :norm, rev=true)

    dummy_idx = findall(feat_importance.param .== "dummy")[1]
    if dummy_idx == 1
        imp = feat_importance[2:11, [:param, :norm]]
    else
        imp = feat_importance[1:min(10, dummy_idx - 1), [:param, :norm]]
    end
    # sort!(imp, :norm, rev=true)

    return imp
end
