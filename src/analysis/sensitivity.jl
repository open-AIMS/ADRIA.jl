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
    Calculates the PAWN sensitivity index.

# Arguments
- X : Model inputs
- Y : Outputs
- S : Number of slides (default: 10)

# Returns
NamedDimsArray, of min, mean, median, max, std, and cv summary statistics.

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
function pawn(X::AbstractArray{<:Real}, Y::Vector{<:Real}, dimnames::Vector{String}; S::Int64=10)::NamedArray
    N, D = size(X)
    step = 1 / S
    seq = 0:step:1

    X_di = zeros(N)
    X_q = zeros(S + 1)
    pawn_t = zeros(S + 1, D)
    results = zeros(D, 6)
    # Hide warnings from HypothesisTests
    with_logger(NullLogger()) do
        for d_i in 1:D
            X_di .= X[:, d_i]
            X_q .= quantile(X_di, seq)
            for s in 2:S
                Y_sel = Y[(X_q[s-1].<X_di).&(X_di.<=X_q[s])]
                if length(Y_sel) == 0
                    pawn_t[s, d_i] = 0.0
                    continue  # no available samples
                end

                pawn_t[s, d_i] = ks_statistic(ApproximateTwoSampleKSTest(Y_sel, Y))
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
function pawn(X::DataFrame, Y::Vector{<:Real}; S::Int64=10)::NamedArray
    return pawn(Matrix(X), Y, names(X); S=S)
end
function pawn(X::NamedArray, Y::Vector{<:Real}; S::Int64=10)::NamedArray
    return pawn(X, Y, names(X, 2); S=S)
end

"""
    rsa(X::DataFrame, y::Vector{<:Real}; S=20)::NamedArray

Perform Regional Sensitivity Analysis.

Regional Sensitivity Analysis is a Monte Carlo Filtering approach which aims to
identify which (group of) factors drive model outputs within or outside of a specified bound.
Outputs which fall inside the bounds are regarded as "behavioral", whereas those outside
are "non-behavioral". The distribution of behavioral/non-behavioral subsets are compared for each factor.
If the subsets are not similar, then the factor is influential.

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
NamedArray, [factor names, upper bound of bins]

# Examples
```julia
ADRIA.sensitivity.rsa(X, y)
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
function rsa(X::DataFrame, y::Vector{<:Real}; S=20)::NamedArray
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

end