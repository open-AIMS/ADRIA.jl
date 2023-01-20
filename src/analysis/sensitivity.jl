module sensitivity

using Logging
using Statistics, Distributions, HypothesisTests, Bootstrap
using NamedArrays, DataFrames

import ADRIA.analysis: normalize


"""
    ks_statistic(ks)

Calculate the Kolmogorov-Smirnov test statistic.
"""
function ks_statistic(ks)
    n = (ks.n_x * ks.n_y) / (ks.n_x + ks.n_y)

    return sqrt(n) * ks.δ
end


"""
    pawn(X::AbstractArray{T}, y::Vector{T}, dimnames::Vector{String}; S::Int64=10)::NamedArray{T} where {T<:Real}

Calculates the PAWN sensitivity index.

The PAWN method (by Pianosi and Wagener) is a moment-independent approach to Global Sensitivity Analysis.
Outputs are characterized by their Cumulative Distribution Function (CDF), quantifying the variation in
the output distribution after conditioning an input over "slices" (\$S\$) - the conditioning intervals.
If both distributions coincide at all slices (i.e., the distributions are similar or identical), then
the factor is deemed non-influential.

This implementation applies the Kolmogorov-Smirnov test as the distance measure and returns summary
statistics (min, mean, median, max, std, and cv) over the slices.

# Arguments
- `X` : Model inputs
- `y` : Model outputs
- `S` : Number of slides (default: 10)

# Returns
NamedArray, of min, mean, median, max, std, and cv summary statistics.

# References
1. Pianosi, F., Wagener, T., 2018.
   Distribution-based sensitivity analysis from a generic input-output sample.
   Environmental Modelling & Software 108, 197-207.
   https://doi.org/10.1016/j.envsoft.2018.07.019

2. Baroni, G., Francke, T., 2020.
   GSA-cvd
   Combining variance- and distribution-based global sensitivity analysis
   https://github.com/baronig/GSA-cvd
"""
function pawn(X::AbstractArray{T}, y::Vector{T}, dimnames::Vector{String}; S::Int64=10)::NamedArray{T} where {T<:Real}
    N, D = size(X)
    step = 1 / S
    seq = 0:step:1

    X_di = zeros(N)
    X_q = zeros(S + 1)
    pawn_t = zeros(S, D)
    results = zeros(D, 6)
    # Hide warnings from HypothesisTests
    with_logger(NullLogger()) do
        for d_i in 1:D
            X_di .= X[:, d_i]
            X_q .= quantile(X_di, seq)
            for s in 1:S
                Y_sel = y[(X_q[s].<X_di).&(X_di.<=X_q[s+1])]
                if length(Y_sel) == 0
                    pawn_t[s, d_i] = 0.0
                    continue  # no available samples
                end

                pawn_t[s, d_i] = ks_statistic(ApproximateTwoSampleKSTest(Y_sel, y))
            end

            p_ind = pawn_t[:, d_i]
            p_mean = mean(p_ind)
            p_sdv = std(p_ind)
            p_cv = p_sdv ./ p_mean
            results[d_i, :] .= (minimum(p_ind), p_mean, median(p_ind), maximum(p_ind), p_sdv, p_cv)
        end
    end

    replace!(results, NaN => 0.0, Inf => 0.0)

    return NamedArray(results, (dimnames, [:min, :mean, :median, :max, :std, :cv]))
end
function pawn(X::DataFrame, y::Vector{T}; S::Int64=10)::NamedArray{T} where {T<:Real}
    return pawn(Matrix(X), y, names(X); S=S)
end
function pawn(X::NamedArray, y::Vector{T}; S::Int64=10)::NamedArray{T} where {T<:Real}
    return pawn(X, y, names(X, 2); S=S)
end

"""
    tsa(X::DataFrame, y::AbstractMatrix)::NamedArray

Perform Temporal (or time-varying) Sensitivity Analysis using the PAWN sensitivity index.

The sensitivity index value for time \$t\$ is inclusive of all time steps prior to \$t\$.
Alternate approaches use a moving window, or only data for time \$t\$.

# Arguments
- `X` : Scenario specification
- `y` : scenario outcomes over time

# Returns
NamedArray, of shape \$D\$ ⋅ 6 ⋅ \$T\$, where
- \$D\$ is the number of dimensions/factors
- 6, corresponds to the min, mean, median, max, std, and cv of the PAWN indices
- \$T\$, the number of time steps

# Examples
```julia
rs = ADRIA.load_results("a ResultSet of interest")

# Get scenario outcomes over time (shape: `time ⋅ scenarios`)
y_tac = ADRIA.metrics.scenario_total_cover(rs)

# Calculate sensitivity of outcome to factors for each time step
ADRIA.sensitivity.tsa(rs.inputs, y_tac)
```
"""
function tsa(X::DataFrame, y::AbstractMatrix{T})::NamedArray{T} where {T<:Real}
    t_pawn_idx = NamedArray(
        zeros(ncol(X), 6, size(y, 1)),
        (names(X), ["min", "mean", "median", "max", "std", "cv"], string.(1:size(y, 1))),
        ("factors", "pawn", "timesteps")
    )

    for t in axes(y, 1)
        t_pawn_idx[:, :, t] .= normalize(
            pawn(X, vec(mean(y[1:t, :], dims=1)))
        )
    end

    return t_pawn_idx
end

"""
    rsa(X::DataFrame, y::Vector{<:Real}; S=20)::NamedArray

Perform Regional Sensitivity Analysis.

Regional Sensitivity Analysis is a Monte Carlo Filtering approach which aims to
identify which (group of) factors drive model outputs within or outside of a specified bound.
Outputs which fall inside the bounds are regarded as "behavioral", whereas those outside
are "non-behavioral". The distribution of behavioral/non-behavioral subsets are compared for each factor.
If the subsets are not similar, then the factor is influential. The sensitivity index is simply the
maximum distance between the two distributions, with larger values indicating greater sensitivity.

The implemented approach slices factor space into \$S\$ bins and iteratively assesses
behavioral and non-behavioral subsets with the non-parametric \$k\$-sample Anderson-Darling test.
Larger values indicate greater dissimilarity (thus, sensitivity). The Anderson-Darling test
places more weight on the tails compared to the Kolmogorov-Smirnov test.

RSA can indicate where in factor space model sensitivities may be, and contributes to a
Value-of-Information (VoI) analysis.

Increasing the value of \$S\$ increases the granularity of the analysis.

# Arguments
- `X` : scenario specification
- `y` : scenario outcomes
- `S` : number of bins to slice factor space into (default: 20)

# Returns
NamedArray, [bin values, factors]

# Examples
```julia
ADRIA.sensitivity.rsa(X, y; S=20)
```

# References
1. Pianosi, F., K. Beven, J. Freer, J. W. Hall, J. Rougier, D. B. Stephenson, and
   T. Wagener. 2016.
   Sensitivity analysis of environmental models:
   A systematic review with practical workflow.
   Environmental Modelling & Software 79:214-232.
   https://dx.doi.org/10.1016/j.envsoft.2016.02.008

2. Saltelli, A., M. Ratto, T. Andres, F. Campolongo, J. Cariboni, D. Gatelli,
   M. Saisana, and S. Tarantola. 2008.
   Global Sensitivity Analysis: The Primer.
   Wiley, West Sussex, U.K.
   https://dx.doi.org/10.1002/9780470725184
   Accessible at: http://www.andreasaltelli.eu/file/repository/Primer_Corrected_2022.pdf
"""
function rsa(X::DataFrame, y::Vector{T}; S=20)::NamedArray{T} where {T<:Real}
    factor_names = names(X)
    N, D = size(X)
    step = 1 / S
    seq = 0:step:1

    X_di = zeros(N)
    X_q = zeros(S + 1)
    r_s = NamedArray(zeros(S, D))

    setnames!(r_s, string.(collect(seq)[2:end]), 1)  # label with upper bounds
    setnames!(r_s, factor_names, 2)
    setdimnames!(r_s, [:bins, :factors])

    for d_i in 1:D
        X_di .= X[:, d_i]
        X_q .= quantile(X_di, seq)
        Threads.@threads for s in 2:S
            sel = (X_q[s-1] .< X_di) .& (X_di .<= X_q[s])
            Y_sel = y[sel]
            if length(Y_sel) == 0
                r_s[s, d_i] = 0.0
                continue  # no available samples
            end

            r_s[s, d_i] = KSampleADTest(Y_sel, y[Not(sel)]).A²k
        end
    end

    return r_s
end


"""
    outcome_map(X::DataFrame, y::AbstractVecOrMat, rule, target_factors::Vector; S::Int=20, n_boot::Int=100, conf::Float64=0.95)::NamedArray

Map normalized outcomes (defined by `rule`) to factor values discretized into `S` bins.

Produces a matrix indicating the range of (normalized) outcomes across factor space for
each dimension (the model inputs). This is similar to a Regional Sensitivity Analysis,
except that the model outputs are examined directly as opposed to a measure of sensitivity.

Note:
- `y` is normalized on a per-column basis prior to the analysis
- Empty areas of factor space (those that do not have any desired outcomes)
  will be assigned `NaN`

# Arguments
- `X` : scenario specification
- `y` : Vector or Matrix of outcomes corresponding to scenarios in `X`
- `rule` : a callable defining a "desirable" scenario outcome
- `target_factors` : list of factors of interest to perform analyses on
- `S` : number of slices of factor space. Higher values equate to finer granularity
- `n_boot` : number of bootstraps (default: 100)
- `conf` : confidence interval (default: 0.95)

# Returns
3-dimensional NamedMatrix, of shape \$S\$ ⋅ \$D\$ ⋅ 3, where:
- \$S\$ is the slices,
- \$D\$ is the number of dimensions, with
- boostrapped mean (dim 1) and the lower/upper 95% confidence interval (dims 2 and 3).

# Examples
```julia
# Factors of interest
foi = [:SRM, :fogging, :a_adapt]

# Find scenarios where all metrics are above their median
rule = y -> all(y .> 0.5)

# Map input values where to their outcomes
ADRIA.sensitivity.outcome_map(X, y, rule, foi; S=20, n_boot=100, conf=0.95)
```
"""
function outcome_map(X::DataFrame, y::AbstractVecOrMat{T}, rule, target_factors::Vector; S::Int=20, n_boot::Int=100, conf::Float64=0.95)::NamedArray{T} where {T<:Real}
    step_size = 1 / S
    steps = collect(0.0:step_size:1.0)

    p_table = NamedArray(fill(NaN, length(steps) - 1, length(target_factors), 3))
    setnames!(p_table, ["$(round(i, digits=2))" for i in steps[2:end]], 1)
    setnames!(p_table, String.(target_factors), 2)
    setnames!(p_table, ["mean", "lower", "upper"], 3)
    setdimnames!(p_table, [:bins, :factors, :CI])

    # Normalize each column in y
    y = normalize(y)

    all_p_rule = findall(rule, eachrow(y))
    num_p_rule = length(all_p_rule)
    if num_p_rule == 0
        @info "Empty result set"
        return p_table
    end

    n_scens = size(X, 1)
    behave = zeros(Bool, n_scens)
    behave[all_p_rule] .= true

    X_q = zeros(S + 1)
    for (j, fact_t) in enumerate(target_factors)
        X_q .= quantile(X[:, fact_t], steps)
        for (i, s) in enumerate(X_q[1:end-1])
            b = (X_q[i] .< X[:, fact_t]) .& (X[:, fact_t] .<= X_q[i+1]) .& behave
            if count(b) == 0
                continue
            end

            bs = bootstrap(mean, y[b], BalancedSampling(n_boot))
            ci = confint(bs, PercentileConfInt(conf))[1]

            p_table[i, j, 1] = ci[1]
            p_table[i, j, 2] = ci[2]
            p_table[i, j, 3] = ci[3]
        end
    end

    return p_table
end

end
