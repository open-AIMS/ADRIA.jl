"""
    Factor(val; kwargs...)::Param

Convenience constructor for Param type with ADRIA-specific metadata.

# Arguments
- `val` : Default value for the factor
- `kwargs` : additional metadata for the factor

# Required keyword arguments
- `ptype` : String, distribution parameter type.
    - "continuous" [default]
    - "ordered categorical"
    - "unordered categorical"
    - "ordered discrete"
    - "discrete"
- `name` : String, Human-readable factor name
- `description` : String, Description of what the factor is
- `dist` : Distribution, as defined by `Distributions.jl` or one of the custom
    constructors defined by ADRIA
- `dist_params` : Tuple or Vector, parameters for `dist`. Stored internally as
    `Vector{Float64}`.

# Optional keyword arguments
- `default_dist_params` : Tuple or Vector, parameters for the default `Distribution`
    type. Defaults to `dist_params` if not provided. Stored internally as
    `Vector{Float64}`.

# Returns
Param
"""
function Factor(val; kwargs...)::Param
    nt = NamedTuple(kwargs)
    _check_has_required_info(nt)
    nt = _set_factor_defaults(nt)

    # Normalize dist_params to Vector{Float64} so all Param instances
    # share one concrete type regardless of tuple arity or element type.
    nt = (;
        nt...,
        dist_params=collect(Float64, nt.dist_params),
        default_dist_params=collect(Float64, nt.default_dist_params)
    )

    return Param((; val=val, nt...))
end

struct CategoricalDistribution{T}
    categories::Vector{T}
    distribution::Categorical
end

"""
    CategoricalDistribution(categories::Vector{T}, weights::Vector{Float64})::CategoricalDistribution{T} where {T}
    CategoricalDistribution(categories::Vector{T})::CategoricalDistribution{T} where {T}
    CategoricalDistribution(categories...)::CategoricalDistribution

Construct a categorical variable. Default to a uniform categorical variable if the
probability weightings are not provided.
"""
function CategoricalDistribution(
    categories::Vector{T}, weights::Vector{Float64}
)::CategoricalDistribution{T} where {T}
    if length(unique(categories)) != length(categories)
        throw(ArgumentError("Categories in a categorical variable must be unique."))
    end

    if length(categories) != length(weights)
        msg = "Length of categories and weightings do not match."
        msg *= " Got $(length(categories)) and $(length(weights))."
        throw(ArgumentError(msg))
    end

    if !(sum(weights) ≈ 1) && all(weights .>= 0.0)
        throw(ArgumentError("Weights must sum to one."))
    end

    return CategoricalDistribution(
        categories,
        Categorical(weights))
end
function CategoricalDistribution(
    categories::Vector{T}
)::CategoricalDistribution{T} where {T}
    n_categories::Int64 = length(categories)
    return CategoricalDistribution(categories, fill(1 / n_categories, n_categories))
end
function CategoricalDistribution(categories...)::CategoricalDistribution
    cats = collect(categories)
    return CategoricalDistribution(cats)
end

function Distributions.quantile(dist::CategoricalDistribution{T}, q::Real)::T where {T}
    underlying_idx::Int64 = Distributions.quantile(dist.distribution, q)
    return dist.categories[underlying_idx]
end

"""
    lower_bound(dist::CategoricalDistribution{T})::T where {T<:Real}

Return the lower bound of the underlying categories if they are numerical.
"""
function lower_bound(dist::CategoricalDistribution{T})::T where {T<:Real}
    return minimum(dist.categories)
end
"""
    upper_bound(dist::CategoricalDistribution{T})::T where {T<:Real}

Return the upper bound of the unerling categories if they are numerical.
"""
function upper_bound(dist::CategoricalDistribution{T})::T where {T<:Real}
    return maximum(dist.categories)
end

# Some distributions have in built methods for bounds, however that would involve
# constructing the distribution

function distribution_lower_bound(::Type{CategoricalDistribution}, dist_params)::Float64
    return minimum(dist_params)
end
function distribution_upper_bound(::Type{CategoricalDistribution}, dist_params)::Float64
    return maximum(dist_params)
end
function distribution_lower_bound(
    ::Type{T}, dist_params
)::Float64 where {T<:Union{DiscreteUniform,Uniform,TriangularDist}}
    return first(dist_params)
end
function distribution_upper_bound(
    ::Type{T}, dist_params
)::Float64 where {T<:Union{DiscreteUniform,Uniform,TriangularDist}}
    return getindex(dist_params, 2)
end
function distribution_lower_bound(
    ::T, dist_params
)::Float64 where {T}
    return first(dist_params)
end
function distribution_upper_bound(
    ::T, dist_params
)::Float64 where {T}
    return getindex(dist_params, 2)
end

function factor_lower_bounds(factor::DataFrameRow)::Float64
    return distribution_lower_bound(factor.dist, factor.dist_params)
end
function factor_upper_bounds(factor::DataFrameRow)::Float64
    return distribution_upper_bound(factor.dist, factor.dist_params)
end

function _set_factor_defaults(kwargs::NT) where {NT<:NamedTuple}
    missing_defaults = (; default_dist_params=kwargs.dist_params)

    for k in keys(missing_defaults)
        if !haskey(kwargs, k)
            kwargs = (; kwargs..., k => missing_defaults[k])
        end
    end

    return kwargs
end

function _check_has_required_info(kwargs::NT) where {NT<:NamedTuple}
    @assert haskey(kwargs, :ptype) "Missing factor field `ptype`"
    @assert haskey(kwargs, :dist) "Missing factor field `dist`"
    @assert haskey(kwargs, :dist_params) "Missing factor field `dist_params`"
    @assert haskey(kwargs, :name) "Missing factor field `name`"
    @assert haskey(kwargs, :description) "Missing factor field `description`"

    param_dist_types = [
        "continuous",
        "ordered categorical",
        "unordered categorical",
        "ordered discrete",
        "discrete"
    ]
    @assert any(occursin.(kwargs[:ptype], param_dist_types)) "`ptype` field is not one of $(param_dist_types)"
end

"""
    DiscreteTriangularDist(lb::Int64, ub::Int64, peak::Int64)::DiscreteNonParametric

Create a discrete triangular distribution.

# Extended Help

The discrete probabilities are estimated using the Cumulative Density Function (CDF),
where the difference between \$F(x)\$, where \$x\$ are the discrete values is used to
infer their probabilities.

We account for the lowest bound being assigned a zero probability by extending the lower
bound by 1, such that the distribution domain is `lb-1 ≤ x ≤ ub`.

Illustrating this principle with code:

```julia
# Triangular Distribution between 2 and 6, with a peak/mode at 3.
d = TriangularDist(2, 6, 3)

# Note that the lowest discrete value will be 0 probability.
cdf.(d, 2:6)
# 5-element Vector{Float64}:
#  0.0
#  0.25
#  0.6666666666666667
#  0.9166666666666666
#  1.0

# Extend the lower bound by 1 instead.
d = TriangularDist(1, 6, 3)

# Now we can take the difference between each F(x) to get their probabilities
cdf.(d, 1:6)
# 6-element Vector{Float64}:
#  0.0
#  0.1
#  0.4
#  0.7333333333333334
#  0.9333333333333333
#  1.0

diff(cdf.(d, 1:6))
# 5-element Vector{Float64}:
#  0.1
#  0.30000000000000004
#  0.33333333333333337
#  0.19999999999999996
#  0.06666666666666665

sum(diff(cdf.(d, 1:6)))
# 1.0
```
"""
function DiscreteTriangularDist(
    lb::T,
    ub::T,
    peak::T
)::DiscreteNonParametric where {T<:Union{Int64,Float64}}
    # The lower bound will always resolve to 0 probability
    # so we extend the lower bound to capture the edge.
    _lb, ub, peak = trunc.(Int64, [lb - 1, ub, peak])
    dist = TriangularDist(_lb, ub, peak)

    # Approximate discrete probabilities using extended lower bound.
    # We get the probabilities for each discrete option.
    probas = diff(cdf.(dist, _lb:ub))

    # Use the original bounds to create the discrete probability distribution.
    return DiscreteNonParametric(lb:ub, probas)
end

"""
    DiscreteOrderedUniformDist

Callable struct representing a uniform ordered categorical distribution factory.
Stored as `dist=DiscreteOrderedUniformDist` in Factor definitions.

Using a struct (DataType) rather than a plain function ensures `typeof(DiscreteOrderedUniformDist)`
is `DataType` — the same metatype as `DiscreteUniform` — so both share one concrete `Param`
type in the NamedTuple type parameter, reducing Flatten.jl specialization count.

# Example
```julia
# Create an ordered categorical distribution
d = DiscreteOrderedUniformDist(0.0, 10.0, 0.5)

# Sampling from the above distribution returns a vector with values ∈ [0, 10]
# including discrete steps of 0.5, e.g.
# [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, ... 9.0, 9.5, 10.0]
```
"""
struct DiscreteOrderedUniformDist end

"""
    DiscreteOrderedUniformDist(lb::T, ub::T, step::Float64)

Creates a uniform ordered categorical distribution, allowing for real-valued options.
Principally used for depth-related factors.

# Arguments
- `lb` : lower bound
- `ub` : upper bound
- `step` : step interval between `lb` and `ub`
"""
function DiscreteOrderedUniformDist(
    lb::T,
    ub::T,
    step::Float64
)::DiscreteNonParametric where T
    options = lb:step:ub
    n_opts = length(options)

    return DiscreteNonParametric(options, fill(1.0 / n_opts, n_opts))
end

"""
    _build_dist_types(vary_vars::DataFrame)

Construct a concretely-typed vector of distribution instances from the `dist`
and `dist_params` columns of a model spec DataFrame.

One instance per unique constructor is built to determine the Union type, then
all rows are constructed into a `Vector{T}` where `T` is that Union. This allows
Julia to compile specialised methods (e.g. for `quantile`) via union-splitting
rather than falling back to `Vector{Any}`.
"""
function _build_dist_types(vary_vars::DataFrame)
    unique_rows = unique(vary_vars, :dist)
    T = Union{(typeof(row.dist(row.dist_params...)) for row in eachrow(unique_rows))...}
    return T[row.dist(row.dist_params...) for row in eachrow(vary_vars)]
end

"""
    _quantile_scale(vary_dists::Vector{T}, u_samples::Matrix) where {T}

Function barrier to scale uniform samples to their target distributions using
the inverse CDF method.

Julia compiles a specialization of this method for each concrete `T` it
encounters. Without this barrier, `quantile` would be compiled generically
against `Vector{Any}` and recompiled on each call.

# Arguments
- `vary_dists` : Vector of distribution instances, one per uncertain parameter.
- `u_samples` : Matrix of uniform `[0, 1]` samples (n_params × n_samples).

# Returns
Matrix of samples scaled to each parameter's target distribution (n_samples × n_params).
"""
function _quantile_scale(vary_dists::Vector{T}, u_samples::Matrix) where {T}
    return Matrix(Distributions.quantile.(vary_dists, u_samples)')
end
