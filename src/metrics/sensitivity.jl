module sensitivity

using Logging
using Statistics, Distributions, HypothesisTests, NamedArrays


"""
    relative_importance(x)

Normalize metrics such that the values ∈ [0, 1].

# Arguments
- x : metric values

# Returns
Normalized values of same shape as x.
"""
function relative_importance(x::AbstractVector{<:Real})
    if (maximum(x) - minimum(x)) == 0.0
        return 0.0
    end

    return (x .- minimum(x)) ./ (maximum(x) - minimum(x))
end
function relative_importance(x::AbstractArray)
    S = copy(x)
    for idx in axes(x, 2)
        S[:, idx] .= relative_importance(x[:, idx])
    end

    return S
end


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
function pawn(X::AbstractArray{<:Real}, Y::Vector{<:Real}; S=10)::NamedArray
    N, D = size(X)
    step = 1 / S

    X_di = zeros(N)
    X_q = zeros(S)
    pawn_t = zeros(S, D)
    results = zeros(D, 6)
    for d_i in 1:D
        seq = 0:step:1-step  # To check

        X_di .= X[:, d_i]
        X_q .= quantile(X_di, seq)
        for s in 1:S-1
            Y_sel = Y[(X_di.>=X_q[s]).&(X_di.<X_q[s+1])]
            if length(Y_sel) == 0
                continue  # no available samples
            end

            # Hide warnings from HypothesisTests
            with_logger(NullLogger()) do
                pawn_t[s, d_i] = ks_statistic(ApproximateTwoSampleKSTest(Y_sel, Y))
            end
        end

        p_ind = pawn_t[:, d_i]
        p_min = minimum(p_ind)
        p_mean = mean(p_ind)
        p_med = median(p_ind)
        p_max = maximum(p_ind)
        p_sdv = std(p_ind)
        p_cv = p_sdv ./ p_mean
        results[d_i, :] .= [p_min, p_mean, p_med, p_max, p_sdv, p_cv]
    end

    replace!(results, NaN => 0.0, Inf => 0.0)

    return NamedArray(results, (1:D, [:min, :mean, :median, :max, :std, :cv]))
end

end