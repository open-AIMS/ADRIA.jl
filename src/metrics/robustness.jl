module robustness

using Statistics, Distributions


"""
    norm(vals::AbstractArray{<:Real})

Normalizes values between 0 and 1, where 1 is the maximum value found.
"""
function norm(vals::AbstractArray{<:Real})
    if (maximum(vals) - minimum(vals)) == 0.0
        return 0.0
    end

    return (vals .- minimum(vals)) ./ (maximum(vals) - minimum(vals))
end


"""
    gini(vals::AbstractVector{<:Real})::Float64 
    gini(vals::AbstractArray{<:Real, 2})

Gini coefficient.

Lower values is greater consistency/equality.
Higher values indicates greater variability/inequality.

# References
1. Hurley, N. P., & Rickard, S. T. (2009). 
   Comparing Measures of Sparsity (arXiv:0811.4706). 
   arXiv. http://arxiv.org/abs/0811.4706

2. https://en.wikipedia.org/wiki/Gini_coefficient#Generalized_inequality_indices
"""
function gini(vals::AbstractVector{<:Real})::Float64
    sv = sort(vals)
    n = length(vals)
    g = (2.0 * sum([x * i for (i, x) in enumerate(sv)]) / sum(sv) - (n + 1)) / n
    if isnan(g)
        return 0.0
    end

    return g
end
function gini(vals::AbstractArray{<:Real,2})
    return gini.(eachcol(vals))
end


"""
    temporal_variability(x::AbstractVector{<:Real})
    temporal_variability(x::AbstractArray{<:Real, 2})
    temporal_variability(x::AbstractArray{<:Real}, func_or_data...)

The V meta-metric.

As a meta-metric, it can be applied to any combination of
metrics (including itself), assuming \$x\$ is bound between 0 and 1.
If this is not the case, consider normalizing values first.

# Examples
```julia-repl
# Apply V to a time series
julia> temporal_variability(rand(50))

# Apply V to an ensemble of of time series
julia> x = rand(50, 200)
julia> temporal_variability(x)

# Create and apply a modified V metric to an ensemble of of time series.
# Where the argument is not a function, the data is used directly
# and so it is assumed all matrices are of the same size and shape.
julia> temporal_variability(x, temporal_variabilty, temporal_variability(P(x)))
julia> temporal_variability(x, temporal_variabilty, P(x), D(x), E(x))
```
"""
function temporal_variability(x::AbstractVector{<:Real})
    return (mean(x) + (1.0 .- gini(x))) ./ 2.0
end
function temporal_variability(x::AbstractArray{<:Real,2})
    return temporal_variability.(eachcol(x))
end
function temporal_variability(x::AbstractArray{<:Real}, func_or_data...)
    return mean(map(f -> f isa Function ? f(x) : f, func_or_data))
end


"""
    intervention_effort(ms, inputs_i)

Obtain an indication of intervention effort for each scenario.
This is referred to as \$F\$.

# Arguments
- ms : model spec
- inputs_i : inputs used for scenarios of interest
"""
function intervention_effort(X, ub, lb)
    return (X .- lb) ./ (ub .- lb)
end
function intervention_effort(ms::DataFrame, X::DataFrame)
    interv_cols = ["seed_TA", "seed_CA", "fogging", "SRM"]
    interv_s = ms[findall(in(interv_cols), ms.fieldname), ["fieldname", "lower_bound", "upper_bound"]]
    ub = interv_s[:, "upper_bound"]
    lb = interv_s[:, "lower_bound"]

    return hcat([intervention_effort(values(X[:, interv_cols[i]]), ub[i], lb[i])
                 for i in eachindex(interv_cols)]...)
end


"""
    intervention_diversity(ms, inputs_i)

Obtain an indication of intervention diversity for a scenario.
Higher values indicate a greater of mix of interventions options were applied.

This is referred to as \$D\$.

# Arguments
- ms : model spec
- inputs_i : inputs used for scenarios of interest
"""
function intervention_diversity(ms, inputs_i)
    return mean(gini(intervention_effort(ms, inputs_i)))
end


"""
    environmental_diversity(ms, inputs_i)

Obtain an indication of environmental factor diversity for a scenario set.
Higher values indicate a greater of mix of environmental conditions were
experienced between scenarios.

This is referred to as \$E\$.

# Arguments
- ms : model spec
- inputs_i : inputs used for scenarios of interest
"""
function environmental_diversity(ms, inputs_i)
    env_cols = ADRIA.component_params(ms, ADRIA.EnvironmentalLayer).fieldname

    env_s = ms[findall(in(env_cols), ms.fieldname), ["fieldname", "lower_bound", "upper_bound"]]
    ub = env_s[:, "upper_bound"]
    lb = env_s[:, "lower_bound"]

    Et = gini(Matrix((inputs_i[:, env_cols] .- lb') ./ (ub .- lb)'))
    if all(isnan.(Et))
        return 0.0
    elseif any(isnan.(Et))
        replace!(Et, NaN => 0.0)
    end

    return mean(Et)
end


end  # module
