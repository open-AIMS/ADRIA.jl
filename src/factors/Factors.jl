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
    - "discrete"
- `name` : String, Human-readable factor name
- `description` : String, Description of what the factor is
- `dist` : Distribution, as defined by `Distributions.jl` or one of the custom
    constructors defined by ADRIA
- `dist_params` : Tuple, of parameters for `dist`

# Optional keyword arguments
- `default_dist_params` : Tuple, parameters for `Distribution` type provided by `Distributions.jl`
- `criteria_keywords` : "tags" specific to CriteriaWeights to indicate what intervention/criteria type the factor is related to.

# Returns
Parameter
"""
function Factor(val; kwargs...)::Param
    nt = NamedTuple(kwargs)
    _check_has_required_info(nt)
    nt = _set_factor_defaults(nt)

    return Param((; val=val, nt...))
end

function _set_factor_defaults(kwargs::NT) where {NT <: NamedTuple}
    missing_defaults = (;
        default_dist_params=kwargs.dist_params,
        criteria_keywords=("",)
    )

    for k in keys(missing_defaults)
        if !haskey(kwargs, k)
            kwargs = (; kwargs..., k=>missing_defaults[k])
        end
    end

    return kwargs
end

function _check_has_required_info(kwargs::NT) where {NT <: NamedTuple}
    @assert haskey(kwargs, :ptype) "Missing factor field `ptype`"
    @assert haskey(kwargs, :dist) "Missing factor field `dist`"
    @assert haskey(kwargs, :dist_params) "Missing factor field `dist_params`"
    @assert haskey(kwargs, :name) "Missing factor field `name`"
    @assert haskey(kwargs, :description) "Missing factor field `description`"

    param_dist_types = [
        "continuous", "ordered categorical", "unordered categorical", "discrete"
    ]
    @assert any(occursin.(kwargs[:ptype], param_dist_types)) "`ptype` field is not one of $(param_dist_types)"
end

"""
    DiscreteTriangularDist(lb::Int64, ub::Int64, peak::Int64)::DiscreteNonParametric

Custom constructor to create a discrete triangular distribution.

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
)::DiscreteNonParametric where {T<:Union{Int64, Float64}}
    # The lower bound will always resolve to 0 probability
    # so we extend the lower bound to capture the edge.
    _lb, ub, peak = trunc.(Int64, [lb-1, ub, peak])
    dist = TriangularDist(_lb, ub, peak)

    # Approximate discrete probabilities using extended lower bound.
    # We get the probabilities for each discrete option.
    probas = diff(cdf.(dist, _lb:ub))

    # Use the original bounds to create the discrete probability distribution.
    return DiscreteNonParametric(lb:ub, probas)
end

"""
    DiscreteOrderedUniformDist(lb::T, ub::T, step::Float64)

Creates an uniform ordered categorical distribution, allowing for real-valued options.
Principally used for depth-related factors.

# Arguments
- `lb` : lower bound
- `ub` : upper bound
- `step` : step interval between `lb` and `ub`

# Example
```julia
# Create a ordered categorical distribution
d = DiscreteOrderedUniformDist(0.0, 10.0, 0.5)

# Sampling from the above distribution returns a vector with values ∈ [0, 10]
# including discrete steps of 0.5, e.g.
# [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, ... 9.0, 9.5, 10.0]
```
"""
function DiscreteOrderedUniformDist(
    lb::T,
    ub::T,
    step::Float64
)::DiscreteNonParametric where T
    options = lb:step:ub
    n_opts = length(options)

    return DiscreteNonParametric(options, fill(1.0/n_opts, n_opts))
end
