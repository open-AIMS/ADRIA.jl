using CategoricalArrays
using Distributions

struct CategoricalDistribution
    child::CategoricalVector
    distribution::Categorical
end

function CategoricalDistribution(
    categories::AbstractVector
)::CategoricalDistribution
    n_categories::Int64 = size(categories, 1)
    cat_arr::CategoricalVector = CategoricalVector(
        categories; levels=categories, ordered=true
    )
    return CategoricalDistribution(
        cat_arr,
        Categorical(fill(1/n_categories, n_categories))
    )
end

function CategoricalDistribution(
    categories...
)::CategoricalDistribution
    _categories = collect(categories)
    return CategoricalDistribution(_categories)
end

function Distributions.quantile(categorical::CategoricalDistribution, q::Real)::Int64
    underlying_q::Int64 = Distributions.quantile(categorical.distribution, q)
    return levelcode(categorical.child[underlying_q])
end
