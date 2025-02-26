using CategoricalArrays
using Distributions

struct CategoricalDistribution
    child::CategoricalVector
    distribution::Categorical
end

function OrderedCategoricalVector(categories::AbstractVector)::CategoricalVector
    return CategoricalVector(categories; levels=categories, ordered=true)
end

function CategoricalDistribution(
    categories::CategoricalVector
)::CategoricalDistribution
    n_categories::Int64 = size(categories, 1)
    return CategoricalDistribution(
        categories,
        Categorical(fill(1/n_categories, n_categories))
    )
end

function Distributions.quantile(categorical::CategoricalDistribution, q::Real)::Int64
    underlying_q::Int64 = Distributions.quantile(categorical.distribution, q)
    return levelcode(categorical.child[underlying_q])
end
