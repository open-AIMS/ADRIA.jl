"""
    rsa(X::DataFrame, y::AbstractVector{<:Real};
        top_proportion::Float64=0.9, method::Symbol=:mann_whitney) -> DataFrame
    rsa(X::DataFrame, selection_mask::Union{BitVector,AbstractVector{Bool}};
        method::Symbol=:mann_whitney) -> DataFrame

Rank-based Regional Sensitivity Analysis: score each input factor by how well it
discriminates between high- and low-outcome scenario groups using the Mann-Whitney U test.

rsa is a rank-based, scenario-conditioned sensitivity method -- it answers
"which factors' distributions differ most between high- and low-outcome scenarios?"
This complements PAWN/KS-based methods (which measure how the full output distribution
shifts across the input space) rather than replacing them. Both are forms of sensitivity
analysis; they answer related but distinct questions.

# Primary dispatch: scalar outcomes
- `X`              : Feature matrix (DataFrame, columns = factors)
- `y`              : Scalar outcome per scenario; scenarios above the `top_proportion`
                     quantile form the selected group
- `top_proportion` : Quantile threshold for the high-outcome group (default: 0.9)
- `method`         : Ranking method; only :mann_whitney is implemented

# Escape-hatch dispatch: pre-computed mask
- `X`              : Feature matrix (DataFrame, columns = factors)
- `selection_mask` : Boolean mask (true = selected/"high-outcome" group)
- `method`         : Ranking method

# Returns
DataFrame sorted descending by `prob_superiority`:
- feature          (Symbol)  : factor name
- statistic        (Float64) : raw Mann-Whitney U
- prob_superiority (Float64) : U / (n1 * n2), in [0, 1]
- effect_size      (Float64) : 1 - 2*U / (n1 * n2), in [-1, 1]

HypothesisTests.MannWhitneyUTest applies a normal approximation with tie correction.
A @warn is emitted when more than 20% of values in a feature column are tied, as the
effect_size formula becomes less reliable in that case.
"""
function rsa(
    X::DataFrame,
    y::AbstractVector{<:Real};
    top_proportion::Float64=0.9,
    method::Symbol=:mann_whitney
)::DataFrame
    threshold = quantile(y, top_proportion)
    selection_mask = y .> threshold
    return rsa(X, selection_mask; method=method)
end
function rsa(
    X::DataFrame,
    selection_mask::Union{BitVector,AbstractVector{Bool}};
    method::Symbol=:mann_whitney
)::DataFrame
    n_true = count(selection_mask)
    n_false = count(.!selection_mask)

    if n_true == 0 || n_false == 0
        throw(
            ArgumentError(
                "selection_mask must have at least one true and one false entry"
            )
        )
    end

    n_features = ncol(X)
    features = Symbol.(names(X))
    stat_vals = Vector{Float64}(undef, n_features)
    prob_sup = Vector{Float64}(undef, n_features)
    eff_size = Vector{Float64}(undef, n_features)

    for (idx, feat) in enumerate(names(X))
        col_vals = X[!, feat]

        if any(ismissing, col_vals)
            throw(
                ArgumentError(
                    "feature column " * string(feat) *
                    " contains NaN or missing values;" *
                    " preprocess before calling rsa"
                )
            )
        end

        col_f = Float64.(col_vals)

        if any(isnan, col_f)
            throw(
                ArgumentError(
                    "feature column " * string(feat) *
                    " contains NaN or missing values;" *
                    " preprocess before calling rsa"
                )
            )
        end

        if length(unique(col_f)) == 1
            @warn "Feature " * string(feat) *
                " has zero variance; returning sentinel values."
            stat_vals[idx] = 0.0
            prob_sup[idx] = 0.5
            eff_size[idx] = 0.0
            continue
        end

        n_total = length(col_f)
        val_counts = Dict{Float64,Int}()
        for v in col_f
            val_counts[v] = get(val_counts, v, 0) + 1
        end
        n_tied = sum(c for c in values(val_counts) if c > 1; init=0)
        if n_tied / n_total > 0.2
            @warn "Feature " * string(feat) *
                " has more than 20% tied values; effect size may be unreliable."
        end

        group1 = col_f[selection_mask]
        group2 = col_f[.!selection_mask]

        test = MannWhitneyUTest(group1, group2)
        U = test.U
        n1 = Float64(length(group1))
        n2 = Float64(length(group2))

        stat_vals[idx] = U
        prob_sup[idx] = U / (n1 * n2)
        eff_size[idx] = 1.0 - (2.0 * U) / (n1 * n2)
    end

    result = DataFrame(;
        feature=features,
        statistic=stat_vals,
        prob_superiority=prob_sup,
        effect_size=eff_size
    )

    sort!(result, :prob_superiority; rev=true)
    return result
end
