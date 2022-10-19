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
    temporal_variability(x::AbstractVector{<:Real})
    temporal_variability(x::AbstractArray{<:Real})

The V meta-metric.

Assumes \$x\$ is bound between 0 and 1.

If this is not the case, consider normalizing values first.
"""
function temporal_variability(x::AbstractVector{<:Real})
    return (mean(x) + (1.0 .- gini(x))) ./ 2
end
function temporal_variability(x::AbstractArray{<:Real,2})
    return temporal_variability.(eachcol(x))
end


"""
    intervention_effort(ms, inputs_i)

Obtain an indication of intervention effort for each scenario.
This is referred to as \$F\$.

# Arguments
- ms : model spec
- inputs_i : inputs used for scenarios of interest
"""
function intervention_effort(ms, inputs_i)
    interv_cols = ["seed_TA", "seed_CA", "fogging", "SRM"]
    interv_s = ms[findall(in(interv_cols), ms.fieldname), ["fieldname", "lower_bound", "upper_bound"]]
    ub = interv_s[:, "upper_bound"]
    lb = interv_s[:, "lower_bound"]

    return mean.(eachrow((inputs_i[:, interv_cols] .- lb') ./ (ub .- lb)'))
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
    return mean(1.0 .- gini(intervention_effort(ms, inputs_i)))
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

    Et = 1 .- gini(Matrix((inputs_i[:, env_cols] .- lb') ./ (ub .- lb)'))
    if all(isnan.(Et))
        return 0.0
    elseif any(isnan.(Et))
        replace!(Et, NaN => 0.0)
    end

    return mean(Et)
end


"""
    gini(vals::AbstractVector{<:Real})::Float64 
    gini(vals::AbstractArray{<:Real})

Gini coefficient.

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


end  # module