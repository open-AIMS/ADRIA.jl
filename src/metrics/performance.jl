module performance

using Statistics, Distributions, DataFrames, StatsBase
using ADRIA

"""
    normalize(vals::AbstractArray{<:Real})

Normalize values using feature scaling such that values are bound
between 0 and 1, where 1 is equivalent to the maximum value found.
"""
function normalize(vals::AbstractArray{<:Real})
    if (maximum(vals) - minimum(vals)) == 0.0
        return 0.0
    end

    return (vals .- minimum(vals)) ./ (maximum(vals) - minimum(vals))
end

"""Root Mean Square Error"""
RMSE(obs, sim) = (sum((sim .- obs) .^ 2) / length(sim))^0.5

"""
    probability(vals::AbstractArray{<:Real})

Calculate probability of individual trajectories, given a scenario ensemble \$S\$.
"""
function probability(S::AbstractArray{<:Real})
    return cdf.(fit(Normal, S), S)
end

"""
    gmd(vals::AbstractVector{<:Real})::Float64
    gmd(vals::AbstractMatrix{<:Real})

Gini's Mean Difference.

The absolute mean of all pairwise distances between elements in a given set.

# References
1. La Haye, R., & Zizler, P. (2019).
   The Gini mean difference and variance.
   METRON, 77(1), 43-52.
   https://doi.org/10.1007/s40300-019-00149-2

2. Yitzhaki, S. (2003).
   Gini's Mean difference: A superior measure of variability for non-normal
     distributions.
   Metron - International Journal of Statistics, LXI(2), 285-316.
   https://ideas.repec.org/a/mtn/ancoec/030208.html

3. Kashif, M., Aslam, M., Al-Marshadi, A. H., & Jun, C.-H. (2016).
   Capability Indices for Non-Normal Distribution Using Gini's Mean Difference as Measure of Variability.
   IEEE Access, 4, 7322-7330.
   https://doi.org/10.1109/ACCESS.2016.2620241
"""
function gmd(vals::AbstractVector{<:Real})::Float64
    n = length(vals)
    sv = sort(vals)
    return (2 / (n * (n - 1))) .* sum(([((2 * i) - n - 1) * sv[i] for i in 1:n]))
end
function gmd(vals::AbstractMatrix{<:Real})
    return gmd.(eachcol(vals))
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

# Apply V to an ensemble of time series
julia> x = rand(50, 200)
julia> temporal_variability(x)

# Create and apply a modified V metric to an ensemble of time series.
# Where the argument is an array and not a function, the data is used directly
# and so it is assumed all matrices are of the same size and shape.
julia> temporal_variability(x, temporal_variabilty, temporal_variability(P(x)))
julia> temporal_variability(x, temporal_variabilty, P(x), D(x), E(x))
```
"""
function temporal_variability(x::AbstractVector{<:Real}; w=[0.9, 0.1])
    return mean([median(x), 1.0 - gmd(x)], weights(w))
end
function temporal_variability(x::AbstractArray{<:Real,2}; w=[0.9, 0.1])
    return temporal_variability.(eachcol(x); w=w)
end
function temporal_variability(x::AbstractArray{<:Real}, func_or_data...)
    return mean([map(f -> f isa Function ? f(x) : f, func_or_data)...])
end

"""
    intervention_effort(ms, inputs_i)

Obtain an indication of intervention effort for each scenario and intervention type.
This is referred to as \$F\$.

# Arguments
- ms : model spec
- inputs_i : inputs used for scenarios of interest

# Returns
Matrix of `s * 8`, where `s` is the number of scenarios and columns are:
`N_seed_TA`, `N_seed_CA`, `N_seed_CNA`, `N_seed_SM`, `N_seed_LM`, `fogging`, `SRM`,
`seed_years`, `shade_years`, `fog_years`
"""
function intervention_effort(X, ub, lb)
    return (X .- lb) ./ (ub .- lb)
end
function intervention_effort(ms::DataFrame, X::DataFrame;
    interv_cols=[
        :N_seed_TA,
        :N_seed_CA,
        :N_seed_CNA,
        :N_seed_SM,
        :N_seed_LM,
        :fogging,
        :SRM,
        :seed_years,
        :shade_years,
        :fog_years
    ]
)
    interv_s = ms[
        findall(in(interv_cols), Symbol.(ms.fieldname)),
        ["fieldname", "lower_bound", "upper_bound"]
    ]
    @assert nrow(interv_s) > 0 "No parameters for $(interv_cols) found."

    ub = interv_s[:, "upper_bound"]
    lb = interv_s[:, "lower_bound"]

    return hcat(
        [
            intervention_effort(values(X[:, interv_cols[i]]), ub[i], lb[i])
            for i in eachindex(interv_cols)
        ]...
    )
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
    return mean(gmd(intervention_effort(ms, inputs_i)))
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
    env_cols = Symbol.(ADRIA.component_params(ms, ADRIA.EnvironmentalLayer).fieldname)
    env_s = ms[
        findall(in(env_cols), Symbol.(ms.fieldname)),
        ["fieldname", "lower_bound", "upper_bound"]
    ]
    @assert nrow(env_s) > 0 "No parameters for $(env_cols) found."

    push!(env_s, ["RCP", 26, 85])  # Add lower/upper bound
    push!(env_cols, :RCP)

    ub = env_s[:, "upper_bound"]
    lb = env_s[:, "lower_bound"]

    Et = gmd.(eachcol((inputs_i[:, env_cols] .- lb') ./ (ub .- lb)'))
    if all(isnan.(Et))
        return 0.0
    elseif any(isnan.(Et))
        replace!(Et, NaN => 0.0)
    end

    return mean(mean(Et; dims=1))
end

end  # module
