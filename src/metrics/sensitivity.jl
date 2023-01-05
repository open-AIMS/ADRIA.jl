module sensitivity

using Logging
using Statistics, Distributions, HypothesisTests, NamedArrays, DataFrames


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

    X_di = zeros(N)
    X_q = zeros(S + 1)
    pawn_t = zeros(S + 1, D)
    results = zeros(D, 6)
    # Hide warnings from HypothesisTests
    with_logger(NullLogger()) do
        for d_i in 1:D
            seq = 0:step:1

            X_di .= X[:, d_i]
            X_q .= quantile(X_di, seq)
            for s in 1:S
                Y_sel = Y[(X_di.>=X_q[s]).&(X_di.<X_q[s+1])]
                if length(Y_sel) == 0
                    pawn_t[s, d_i] = 0.0
                    continue  # no available samples
                end

                pawn_t[s, d_i] = ks_statistic(ApproximateTwoSampleKSTest(Y_sel, Y))
            end

            # p_ind = pawn_t[:, d_i]
            p_mean = mean(pawn_t[:, d_i])
            p_sdv = std(pawn_t[:, d_i])
            p_cv = p_sdv ./ p_mean
            results[d_i, :] .= (minimum(pawn_t[:, d_i]), p_mean, median(pawn_t[:, d_i]), maximum(pawn_t[:, d_i]), p_sdv, p_cv)
        end
    end

    replace!(results, NaN => 0.0, Inf => 0.0)

    return NamedArray(results, (dimnames, [:min, :mean, :median, :max, :std, :cv]))
end
function pawn(X::DataFrame, Y::Vector{<:Real}; S::Int64=10)::NamedArray
    return pawn(Matrix(X), Y, names(X); S=S)
end
# function pawn(X::Matrix, Y::Vector{<:Real}; S::Int64=10)::NamedArray
#     return pawn(X, Y, names(X); S=S)
# end

end