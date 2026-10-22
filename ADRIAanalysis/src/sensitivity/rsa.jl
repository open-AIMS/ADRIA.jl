"""
    rsa(X::DataFrame, y::AbstractVector{<:Real};
        top_proportion::Float64=0.9, method::Symbol=:auto) -> DataFrame
    rsa(X::DataFrame, selection_mask::Union{BitVector,AbstractVector{Bool}};
        method::Symbol=:auto) -> DataFrame

Rank-based Regional Sensitivity Analysis: score each input factor by how well it
discriminates between high- and low-outcome scenario groups.

rsa is a rank-based, scenario-conditioned sensitivity method -- it answers
"which factors' distributions differ most between high- and low-outcome scenarios?"
This complements PAWN/KS-based methods (which measure how the full output distribution
shifts across the input space) rather than replacing them. Both are forms of sensitivity
analysis; they answer related but distinct questions.

# Test dispatch per column
Two tests are available:
- **Mann-Whitney U** (`HypothesisTests.MannWhitneyUTest`): appropriate for continuous,
  ordered-categorical, and ordered-discrete factors, where the two groups' *ranks* are
  meaningfully comparable.
- **Cramer's V** (`HypothesisTests.ChisqTest` on a factor-level x group contingency
  table): appropriate for unordered/nominal categorical factors, where there is no
  meaningful notion of "rank".

With `method=:auto` (the default), each column of `X` is dispatched individually based on
its `DataFrames.colmetadata(X, col, "ptype", "continuous")` value: columns tagged
`"unordered categorical"` use Cramer's V; everything else (`"continuous"`,
`"ordered categorical"`, `"ordered discrete"`, `"discrete"`, or no metadata at all) uses
Mann-Whitney. Passing `method=:mann_whitney` or `method=:cramers_v` overrides this and
forces every column through that single test regardless of its `ptype` tag.

**Backwards-compatible fallback**: if `X` has no `colmetadata` attached to *any* column
(e.g. a plain, untagged `DataFrame`) and at least one column's `eltype` looks
non-numeric (a heuristic proxy for "might actually be nominal but untagged"), a single
`@warn` is emitted once before the per-feature loop, noting that no ptype metadata was
found and Mann-Whitney is being applied indiscriminately -- which is statistically
invalid for nominal categorical factors. This is informational only: every untagged
column still falls back to `"continuous"`/Mann-Whitney, matching pre-existing behaviour.

# Primary dispatch: scalar outcomes
- `X`              : Feature matrix (DataFrame, columns = factors)
- `y`              : Scalar outcome per scenario; scenarios above the `top_proportion`
                     quantile form the selected group
- `top_proportion` : Quantile threshold for the high-outcome group (default: 0.9)
- `method`         : `:auto` (default, per-column ptype dispatch), `:mann_whitney`, or
                     `:cramers_v` (the latter two force uniform application to all columns)

# Escape-hatch dispatch: pre-computed mask
- `X`              : Feature matrix (DataFrame, columns = factors)
- `selection_mask` : Boolean mask (true = selected/"high-outcome" group)
- `method`         : Ranking method, as above

# Returns
DataFrame sorted descending by `prob_superiority`:
- feature          (Symbol)  : factor name
- test             (Symbol)  : `:mann_whitney` or `:cramers_v` -- which test produced
                               this row
- statistic        (Float64) : raw Mann-Whitney U, or the chi-squared statistic for
                               Cramer's V rows
- prob_superiority (Float64) : U / (n1 * n2), in [0, 1], for Mann-Whitney rows. Cramer's
                               V has no equivalent "probability of superiority" concept
                               (it is an unsigned association strength, not a rank
                               comparison), so this is `NaN` for `:cramers_v` rows --
                               following the same "not meaningful for this row" convention
                               used for the zero-variance sentinel below, except NaN
                               rather than 0.5 since there is no neutral value to report.
- effect_size      (Float64) : 1 - 2*U / (n1 * n2), in `[-1, 1]` (signed), for Mann-Whitney
                               rows; Cramer's V itself, in `[0, 1]` (unsigned), for
                               `:cramers_v` rows. These two are NOT on a comparable scale
                               -- distinguish them via the `test` column, do not compare
                               `effect_size` values across test types directly.

Because Mann-Whitney and Cramer's V effect sizes are not comparable, a result table
containing BOTH test types does not produce a single meaningfully-ranked-together
ordering: `sort!(result, :prob_superiority; rev=true)` will place all `:cramers_v` rows
(prob_superiority = NaN) according to Julia's NaN sort ordering, not by association
strength. Callers wanting a proper within-statistic ranking should filter by `test`
first (e.g. `filter(:test => ==(:mann_whitney), result)`).

HypothesisTests.MannWhitneyUTest applies a normal approximation with tie correction.
A @warn is emitted when more than 20% of values in a feature column are tied, as the
effect_size formula becomes less reliable in that case.

**Known confound**: for factors like `mcda_method` (sentinel-zeroed/inapplicable when
`guided <= 0`), Cramer's V will partly re-detect "is guided active at all", a signal
`guided`'s own Mann-Whitney effect size already captures separately. This is a known,
expected confound -- it is documented here rather than engineered around.
"""
function rsa(
    X::DataFrame,
    y::AbstractVector{<:Real};
    top_proportion::Float64=0.9,
    method::Symbol=:auto
)::DataFrame
    threshold = quantile(y, top_proportion)
    selection_mask = y .> threshold
    return rsa(X, selection_mask; method=method)
end
function rsa(
    X::DataFrame,
    selection_mask::Union{BitVector,AbstractVector{Bool}};
    method::Symbol=:auto
)::DataFrame
    if !(method in (:auto, :mann_whitney, :cramers_v))
        throw(
            ArgumentError(
                "Unknown method :$(method); choose :auto, :mann_whitney, or :cramers_v."
            )
        )
    end

    n_true = count(selection_mask)
    n_false = count(.!selection_mask)

    if n_true == 0 || n_false == 0
        throw(
            ArgumentError(
                "selection_mask must have at least one true and one false entry"
            )
        )
    end

    if method == :auto && isempty(colmetadatakeys(X))
        looks_nominal = any(!(eltype(X[!, c]) <: Real) for c in names(X))
        if looks_nominal
            @warn "No column \"ptype\" metadata found on X; applying Mann-Whitney " *
                "indiscriminately to all columns. This is statistically invalid for " *
                "nominal/unordered categorical factors -- tag columns with " *
                "colmetadata!(df, col, \"ptype\", ...; style=:note) to enable " *
                "automatic per-column test dispatch."
        end
    end

    n_features = ncol(X)
    features = Symbol.(names(X))
    stat_vals = Vector{Float64}(undef, n_features)
    prob_sup = Vector{Float64}(undef, n_features)
    eff_size = Vector{Float64}(undef, n_features)
    test_vals = Vector{Symbol}(undef, n_features)

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

        ptype = colmetadata(X, feat, "ptype", "continuous")

        use_cramers_v = if method == :mann_whitney
            false
        elseif method == :cramers_v
            true
        else
            ptype == "unordered categorical"
        end

        if use_cramers_v
            chi2, V = _cramers_v(col_vals, selection_mask)
            stat_vals[idx] = chi2
            prob_sup[idx] = NaN
            eff_size[idx] = V
            test_vals[idx] = :cramers_v
            continue
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

        test_vals[idx] = :mann_whitney

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

        mw_test = MannWhitneyUTest(group1, group2)
        U = mw_test.U
        n1 = Float64(length(group1))
        n2 = Float64(length(group2))

        stat_vals[idx] = U
        prob_sup[idx] = U / (n1 * n2)
        eff_size[idx] = 1.0 - (2.0 * U) / (n1 * n2)
    end

    result = DataFrame(;
        feature=features,
        test=test_vals,
        statistic=stat_vals,
        prob_superiority=prob_sup,
        effect_size=eff_size
    )

    sort!(result, :prob_superiority; rev=true)
    return result
end

"""
    _cramers_v(col_vals, selection_mask) -> (statistic::Float64, effect_size::Float64)

Cramer's V association strength between a nominal/unordered-categorical feature column
`col_vals` and a binary `selection_mask`, via `HypothesisTests.ChisqTest` on the
`(n_levels x 2)` contingency table of factor level vs. group membership.

`V = sqrt(chi2 / (n * (min(n_rows, n_cols) - 1)))`, unsigned, in `[0, 1]`. `n_cols` is
always 2 (the two `selection_mask` groups); `n_rows` is the number of unique values in
`col_vals`. When `col_vals` has fewer than 2 unique values, `min(n_rows, n_cols) - 1`
would be zero (division by zero) -- mirroring the zero-variance handling in `rsa`, this
case emits a `@warn` and returns the sentinel `(0.0, 0.0)` instead.
"""
function _cramers_v(
    col_vals::AbstractVector,
    selection_mask::Union{BitVector,AbstractVector{Bool}}
)::Tuple{Float64,Float64}
    levels = unique(col_vals)
    n_rows = length(levels)
    n_cols = 2

    if n_rows < 2
        @warn "Feature has fewer than 2 unique categorical values; " *
            "returning sentinel values."
        return 0.0, 0.0
    end

    level_idx = Dict(lv => i for (i, lv) in enumerate(levels))
    table = zeros(Int, n_rows, n_cols)
    for (v, sel) in zip(col_vals, selection_mask)
        table[level_idx[v], sel ? 1 : 2] += 1
    end

    ct_test = ChisqTest(table)
    chi2 = ct_test.stat
    n = length(col_vals)
    V = sqrt(chi2 / (n * (min(n_rows, n_cols) - 1)))

    return Float64(chi2), Float64(V)
end

const DHW_STAT_COLS = Set((:dhw_mean, :dhw_stdev, :dhw_complexity))

"""
    quantile_strata_edges(vals::AbstractVector{<:Real}, n_strata::Int) -> Vector{Float64}

Equal-frequency quantile bin edges for `vals` (length `n_strata + 1`), with the
final edge nudged up by `nextfloat` so the maximum value falls inside the last
bin: stratum `s` is `edges[s] <= v < edges[s+1]`.

Used internally by `stratified_rsa` to bin `strat_col`; exposed so callers can
derive a matching bin assignment for a *different* vector against the same
edges (e.g. binning counterfactual scenarios' DHW values into the same strata
as the intervention scenarios used to derive `edges`, as `stratified_cf_mask`
does).
"""
function quantile_strata_edges(vals::AbstractVector{<:Real}, n_strata::Int)::Vector{Float64}
    edges = quantile(Float64.(vals), range(0.0, 1.0; length=n_strata + 1))
    edges[end] = nextfloat(edges[end])
    return edges
end

function _stratum_index(v::Real, edges::AbstractVector{<:Real})::Int
    n_strata = length(edges) - 1
    for s = 1:n_strata
        edges[s] <= v < edges[s + 1] && return s
    end
    return v < edges[1] ? 1 : n_strata
end

"""
    stratified_cf_mask(y_iv::AbstractVector{<:Real}, dhw_iv::AbstractVector{<:Real},
                        y_cf::AbstractVector{<:Real}, dhw_cf::AbstractVector{<:Real};
                        n_strata::Int=4, quantile_level::Float64=0.8,
                        min_stratum_n::Int=10) -> BitVector

Behavioural mask for intervention scenarios, with the "behavioural" bar set
*locally per DHW stratum*: scenario `i` is behavioural iff `y_iv[i]` beats the
`quantile_level` quantile of counterfactual outcomes (`y_cf`) drawn from the
SAME DHW stratum as scenario `i`.

Bin edges are the equal-frequency quantile edges of `dhw_iv` (see
`quantile_strata_edges`); `dhw_cf` is binned into those same edges so each
stratum's threshold is computed from counterfactual scenarios under comparable
climate conditions.

This is deliberately local rather than global: a single global
`quantile(y_cf, quantile_level)` threshold, especially at a high
`quantile_level`, is disproportionately clearable only by scenarios drawn under
mild DHW -- passed to `stratified_rsa`, that collapses the "which factors
matter within this DHW band" question back into "which DHW band is mild",
which stratification exists to avoid. Since DHW is fixed within a stratum here,
"beat the local top `1 - quantile_level` of counterfactual outcomes" isolates
intervention-parameter effects from climate severity instead.

Strata with fewer than `min_stratum_n` counterfactual scenarios fall back to
the GLOBAL `quantile_level` quantile of `y_cf`, with a `@warn`.
"""
function stratified_cf_mask(
    y_iv::AbstractVector{<:Real},
    dhw_iv::AbstractVector{<:Real},
    y_cf::AbstractVector{<:Real},
    dhw_cf::AbstractVector{<:Real};
    n_strata::Int=4,
    quantile_level::Float64=0.8,
    min_stratum_n::Int=10
)::BitVector
    if length(y_iv) != length(dhw_iv)
        throw(DimensionMismatch("y_iv and dhw_iv must have the same length."))
    end
    if length(y_cf) != length(dhw_cf)
        throw(DimensionMismatch("y_cf and dhw_cf must have the same length."))
    end

    dhw_iv_f = Float64.(dhw_iv)
    dhw_cf_f = Float64.(dhw_cf)
    y_cf_f = Float64.(y_cf)

    edges = quantile_strata_edges(dhw_iv_f, n_strata)
    global_threshold = quantile(y_cf_f, quantile_level)

    thresholds = Vector{Float64}(undef, n_strata)
    for s = 1:n_strata
        in_stratum = (dhw_cf_f .>= edges[s]) .& (dhw_cf_f .< edges[s + 1])
        n_cf_s = count(in_stratum)
        if n_cf_s < min_stratum_n
            @warn "DHW stratum $s has only $n_cf_s counterfactual scenarios " *
                "(< $min_stratum_n); falling back to the global CF quantile for this stratum."
            thresholds[s] = global_threshold
        else
            thresholds[s] = quantile(y_cf_f[in_stratum], quantile_level)
        end
    end

    mask = falses(length(y_iv))
    for i in eachindex(y_iv)
        s = _stratum_index(dhw_iv_f[i], edges)
        mask[i] = y_iv[i] > thresholds[s]
    end
    return mask
end

"""
    stratified_rsa(fs::DataFrame, y::AbstractVector{<:Real};
                   strat_col::Symbol=:dhw_mean,
                   n_strata::Int=4,
                   top_proportion::Float64=0.9) -> DataFrame
    stratified_rsa(fs::DataFrame, selection_mask::Union{BitVector,AbstractVector{Bool}};
                   strat_col::Symbol=:dhw_mean,
                   n_strata::Int=4) -> DataFrame

Run `rsa` independently within each DHW quantile stratum and return a
long-format DataFrame summarising factor importance across strata.

Scenarios are binned into `n_strata` equal-frequency quantile groups on
`strat_col` (default `:dhw_mean`).  DHW stat columns
(`dhw_mean`, `dhw_stdev`, `dhw_complexity`) are dropped from the feature
matrix within each stratum before calling `rsa`, so they do not confound
the within-stratum ranking.

# Primary dispatch: scalar outcomes
- `fs`             : Feature matrix as returned by `feature_set(rs)`, which
                     must contain `strat_col` as a column.
- `y`              : Scalar outcome per scenario (same row order as `fs`).
- `strat_col`      : Column of `fs` used to form strata (default `:dhw_mean`).
- `n_strata`       : Number of equal-frequency quantile bins (default 4).
- `top_proportion` : Passed through to `rsa`; quantile threshold for the
                     high-outcome group within each stratum (default 0.9).

# Escape-hatch dispatch: pre-computed mask
- `fs`             : Feature matrix, as above.
- `selection_mask` : Boolean mask (true = selected/"behavioural" group), same
                     row order as `fs`. Sliced per stratum and passed straight
                     to `rsa`'s mask dispatch -- the caller decides how the
                     behavioural split is defined (e.g. against a counterfactual
                     threshold) rather than `stratified_rsa` re-deriving a
                     per-stratum quantile from `y` itself.
- `strat_col`, `n_strata` : As above.

# Returns
Long-format `DataFrame` with columns:
- `feature`          (Symbol)  : factor name
- `stratum`          (Int)     : stratum index 1..n_strata (low->high DHW)
- `prob_superiority` (Float64) : Mann-Whitney P(X1 > X2) for this factor in this stratum
- `effect_size`      (Float64) : 1 - 2U/(n1*n2)
- `mean_importance`  (Float64) : mean of abs(prob_superiority - 0.5) across all strata
                                  (same value repeated per row). Measures average deviation
                                  from neutrality; a factor with reversed importance across
                                  strata (e.g. [0.9, 0.1]) scores the same as a consistently
                                  important one ([0.9, 0.9]), both higher than a neutral
                                  factor ([0.5, 0.5]).

Strata that contain fewer than 10 scenarios are skipped with a `@warn`. For the
scalar-outcome dispatch, strata where all outcome values are identical are
also skipped; for the mask dispatch, strata where the sliced mask is all-true
or all-false are skipped instead. Skipped strata are absent from the returned
DataFrame.
"""
function stratified_rsa(
    fs::DataFrame,
    y::AbstractVector{<:Real};
    strat_col::Symbol=:dhw_mean,
    n_strata::Int=4,
    top_proportion::Float64=0.9
)::DataFrame
    return _stratified_rsa(fs, y, strat_col, n_strata) do X_s, y_s
        if length(unique(y_s)) == 1
            @warn "All outcome values are identical in this stratum; skipping."
            return nothing
        end
        return rsa(X_s, y_s; top_proportion=top_proportion)
    end
end
function stratified_rsa(
    fs::DataFrame,
    selection_mask::Union{BitVector,AbstractVector{Bool}};
    strat_col::Symbol=:dhw_mean,
    n_strata::Int=4
)::DataFrame
    return _stratified_rsa(fs, selection_mask, strat_col, n_strata) do X_s, mask_s
        n_true = count(mask_s)
        n_false = count(.!mask_s)
        if n_true == 0 || n_false == 0
            @warn "This stratum's selection_mask is all-true or all-false; skipping."
            return nothing
        end
        return rsa(X_s, mask_s)
    end
end

"""
    _stratified_rsa(compute_si, fs, labels, strat_col, n_strata) -> DataFrame

Shared stratification/aggregation core for `stratified_rsa`. Bins `fs` into
`n_strata` equal-frequency quantile groups on `strat_col`, drops DHW stat
columns and zero-variance columns within each stratum, and calls
`compute_si(X_s, labels_s)` per stratum -- `labels` is either the scalar
outcome vector `y` or a `selection_mask`, sliced to the stratum's rows.
`compute_si` returns a per-stratum `rsa` result DataFrame, or `nothing` to
skip the stratum (after emitting its own `@warn`).
"""
function _stratified_rsa(
    compute_si,
    fs::DataFrame,
    labels::AbstractVector,
    strat_col::Symbol,
    n_strata::Int
)::DataFrame
    if !(strat_col in Symbol.(names(fs)))
        throw(
            ArgumentError(
                "strat_col :$(strat_col) not found in fs; run feature_set(rs) first."
            )
        )
    end
    if length(labels) != nrow(fs)
        throw(
            DimensionMismatch(
                "fs has $(nrow(fs)) rows but labels has $(length(labels)) elements."
            )
        )
    end

    # Drop DHW stat columns from the feature matrix (they stratify, not discriminate)
    drop_cols = [c for c in names(fs) if Symbol(c) in DHW_STAT_COLS]
    X = fs[:, Not(drop_cols)]

    strat_vals = Float64.(fs[!, strat_col])
    edges = quantile_strata_edges(strat_vals, n_strata)

    rows = DataFrame[]
    for s = 1:n_strata
        lo, hi = edges[s], edges[s + 1]
        mask = (strat_vals .>= lo) .& (strat_vals .< hi)
        n_in = count(mask)
        if n_in < 10
            @warn "Stratum $s has only $n_in scenarios (< 10); skipping."
            continue
        end
        labels_s = labels[mask]
        X_s = X[mask, :]
        # Drop zero-variance columns within this stratum to avoid sentinel spam.
        # NOTE: column-selection-by-name (X_s[:, varying]) is expected to preserve
        # :note-style colmetadata (as `Not()` selection is confirmed to), the same as
        # the row-mask X[mask, :] above -- but this has not been runtime-verified
        # against real tagged data; spot-check once real data is available.
        varying = [col for col in names(X_s) if length(unique(X_s[!, col])) > 1]
        si = compute_si(X_s[:, varying], labels_s)
        si === nothing && continue
        si[!, :stratum] .= s
        push!(rows, si)
    end

    isempty(rows) && return DataFrame(;
        feature=Symbol[], stratum=Int[], prob_superiority=Float64[],
        effect_size=Float64[], mean_importance=Float64[]
    )

    combined = vcat(rows...)

    # Compute mean_importance = mean(abs(prob_superiority - 0.5)) across strata per feature.
    # This correctly identifies both consistently-important and DHW-reversed factors,
    # unlike mean(prob_superiority) which would score reversed factors as neutral.
    imp = combine(
        groupby(combined, :feature),
        :prob_superiority => (ps -> mean(abs.(ps .- 0.5))) => :mean_importance
    )
    combined = leftjoin(combined, imp; on=:feature)

    return combined[
        :, [:feature, :stratum, :prob_superiority, :effect_size, :mean_importance]
    ]
end

"""
    _rank_aligned_delta(y_iv, y_cf_raw) -> (Vector{Float64}, Vector{Int})

Sort both outcome vectors ascending and compute element-wise delta at matched
rank positions.  When lengths differ, the CF empirical quantile function is
evaluated at the midpoint probability for each of the `n_iv` rank positions
(`p_r = (r - 0.5) / n_iv`) to avoid boundary extrapolation.

Returns `(y_delta, perm)` where `perm` is the sort permutation applied to
`y_iv`; apply `fs_iv[perm, :]` to keep the feature matrix aligned with
`y_delta`.
"""
function _rank_aligned_delta(
    y_iv::Vector{Float64},
    y_cf_raw::Vector{Float64}
)::Tuple{Vector{Float64},Vector{Int}}
    n_iv = length(y_iv)
    n_cf = length(y_cf_raw)
    perm = sortperm(y_iv)
    y_iv_sorted = y_iv[perm]
    y_cf_sorted = if n_iv == n_cf
        sort(y_cf_raw)
    else
        # Midpoint-rule probabilities: avoids p=0 and p=1 boundary issues.
        probs = ((1:n_iv) .- 0.5) ./ n_iv
        quantile(y_cf_raw, probs)
    end
    return y_iv_sorted .- y_cf_sorted, perm
end

"""
    counterfactual_delta(
        rs_intervention::ResultSet,
        rs_counterfactual::ResultSet,
        metric_fn;
        bootstrap_n::Int=0
    ) -> NamedTuple

Compute per-scenario intervention lift (delta) and return the intervention
feature matrix with DHW columns dropped, ready for `rsa`.

Both `rs_intervention` and `rs_counterfactual` are assumed to be independently
Sobol'-sampled over the same parameter bounds.

Two estimators are available via `delta_method`:

**`:rank_aligned`** (default): Both outcome vectors are sorted ascending; the
CF vector is interpolated to `n_iv` quantile positions if lengths differ.
The lift for rank position `r` is:

    y_delta[r] = y_iv^(r) - y_cf^(r)

`fs_intervention` rows are permuted to the same sort order so that
`y_delta[r]` and `fs_intervention[r, :]` always correspond to the same
intervention scenario.  RSA on this delta asks: "which factors characterise
scenarios that beat the same-rank counterfactual outcome?"

**`:mean_difference`**: The lift for intervention scenario `i` is:

    y_delta[i] = y_iv[i] - mean(y_cf)

Note: because `mean(y_cf)` is a constant offset, RSA on this delta is
mathematically equivalent to RSA on `y_iv` directly. `fs_intervention` row
order is unchanged.

# Arguments
- `rs_intervention`   : ResultSet from an intervention run.
- `rs_counterfactual` : ResultSet from a no-intervention (counterfactual) run.
- `metric_fn`         : Function `rs -> AbstractVector{<:Real}` mapping a
                        ResultSet to a scalar outcome per scenario.
- `delta_method`      : `:rank_aligned` (default) or `:mean_difference`; see above.
- `bootstrap_n`       : Number of bootstrap resamples for a 95% CI on
                        `mean(y_delta)`.  0 (default) skips bootstrapping.
- `fs`                : Precomputed `feature_set(rs_intervention)` to reuse, for
                        callers that invoke `counterfactual_delta` repeatedly
                        against the same `rs_intervention` with different
                        `metric_fn`/location restrictions (e.g. per-period or
                        per-location-set repeats). `feature_set` is a pure
                        function of `rs_intervention` alone -- independent of
                        `metric_fn` -- so recomputing it on every call
                        (deployment-log summaries, DHW stats) is pure overhead
                        once a caller needs more than one delta from the same
                        `rs_intervention`. Must have exactly `length(metric_fn(
                        rs_intervention))` rows, in `rs_intervention`'s original
                        scenario order (i.e. `feature_set(rs_intervention)`
                        itself, unpermuted). `nothing` (default) computes it
                        internally, as before.

# Returns
`NamedTuple` with fields:
- `y_delta`        : `Vector{Float64}` -- per-scenario lift (length = n intervention
                     scenarios).  With `:rank_aligned`, sorted ascending by outcome.
- `fs_intervention`: `DataFrame` -- `feature_set(rs_intervention)` with DHW stat columns
                     (`:dhw_mean`, `:dhw_stdev`, `:dhw_complexity`) removed.
                     With `:rank_aligned`, rows are permuted to match `y_delta` order.
- `iv_perm`        : `Vector{Int}` -- permutation applied to the intervention rows,
                     i.e. `fs_intervention`/`y_delta` row `i` corresponds to the `i`-th
                     original intervention scenario (in `metric_fn(rs_intervention)`
                     order) at index `iv_perm[i]`. Identity for `:mean_difference`.
                     Callers holding other per-original-intervention-row arrays (e.g.
                     a `guided`/`unguided` mask) must index them with `iv_perm` before
                     lining them up against `fs_intervention`/`y_delta`.
- `bootstrap_ci`   : `Union{Nothing, Tuple{Float64,Float64}}` -- bootstrapped 95% CI on
                     `mean(y_delta)`, or `nothing` if `bootstrap_n == 0`

# Edge cases
- If `all(y_delta .== 0)`, a `@warn` is emitted (ATE is zero; RSA will return sentinels).
- If either ResultSet has zero scenarios, throws `ArgumentError`.
"""
function counterfactual_delta(
    rs_intervention::ResultSet,
    rs_counterfactual::ResultSet,
    metric_fn;
    delta_method::Symbol=:rank_aligned,
    bootstrap_n::Int=0,
    fs::Union{Nothing,DataFrame}=nothing
)
    y_iv = Float64.(metric_fn(rs_intervention))
    y_cf_raw = Float64.(metric_fn(rs_counterfactual))

    n_iv = length(y_iv)
    n_cf = length(y_cf_raw)

    if n_iv == 0
        throw(ArgumentError("rs_intervention has zero scenarios."))
    end
    if n_cf == 0
        throw(ArgumentError("rs_counterfactual has zero scenarios."))
    end

    y_delta, sort_perm = if delta_method === :rank_aligned
        _rank_aligned_delta(y_iv, y_cf_raw)
    elseif delta_method === :mean_difference
        y_iv .- mean(y_cf_raw), collect(1:n_iv)
    else
        throw(
            ArgumentError(
                "Unknown delta_method :$(delta_method); " *
                "choose :rank_aligned or :mean_difference."
            )
        )
    end

    if all(y_delta .== 0.0)
        @warn "All y_delta values are zero; intervention has no measurable effect. " *
            "RSA will return sentinel values."
    end

    # Build feature set (or reuse a precomputed one) and drop DHW stat columns
    fs_iv = fs === nothing ? feature_set(rs_intervention) : fs
    if nrow(fs_iv) != n_iv
        throw(
            ArgumentError(
                "Precomputed `fs` has $(nrow(fs_iv)) rows; expected $n_iv " *
                "(length of metric_fn(rs_intervention))."
            )
        )
    end
    drop_cols = [c for c in names(fs_iv) if Symbol(c) in DHW_STAT_COLS]
    fs_iv = fs_iv[sort_perm, Not(drop_cols)]

    # Optional bootstrap CI on mean(y_delta)
    ci = if bootstrap_n > 0
        bs = bootstrap(mean, y_delta, BasicSampling(bootstrap_n))
        _, lo, hi = confint(bs, BasicConfInt(0.95))[1]
        (lo, hi)
    else
        nothing
    end

    return (; y_delta=y_delta, fs_intervention=fs_iv, iv_perm=sort_perm, bootstrap_ci=ci)
end

"""
    counterfactual_delta(
        rs::ResultSet,
        cf_mask::AbstractVector{Bool},
        metric_fn;
        bootstrap_n::Int=0
    ) -> NamedTuple

Variant of `counterfactual_delta` for result sets where intervention and
counterfactual scenarios are stored together (e.g. `rs.inputs.guided .== -1`
marks the counterfactual group).

`metric_fn` is called once with the full `rs`; its output (one value per scenario)
is then split by `cf_mask`.  Intervention scenarios are all rows where `cf_mask`
is `false`.

All other behaviour (`delta_method`, DHW column removal, optional bootstrap CI,
the `iv_perm` return field) is identical to the two-ResultSet overload. Here
`iv_perm` indexes into `findall(.!cf_mask)` order, i.e. the intervention rows
of `rs` in their original order.

`fs` (optional): precomputed `feature_set(rs)` to reuse across repeated calls
against the same `rs` with different `metric_fn`/location restrictions --
see the two-ResultSet overload's docstring for the rationale. Must have
`length(cf_mask)` rows in `rs`'s original scenario order (i.e.
`feature_set(rs)` itself, unpermuted and NOT pre-filtered to `iv_mask`).
`nothing` (default) computes it internally, as before.
"""
function counterfactual_delta(
    rs::ResultSet,
    cf_mask::AbstractVector{Bool},
    metric_fn;
    delta_method::Symbol=:rank_aligned,
    bootstrap_n::Int=0,
    fs::Union{Nothing,DataFrame}=nothing
)
    iv_mask = .!cf_mask

    if !any(iv_mask)
        throw(
            ArgumentError(
                "cf_mask selects all scenarios; no intervention scenarios remain."
            )
        )
    end
    if !any(cf_mask)
        throw(ArgumentError("cf_mask selects no counterfactual scenarios."))
    end

    y_all = Float64.(metric_fn(rs))
    y_iv = y_all[iv_mask]
    y_cf_raw = y_all[cf_mask]

    n_iv = length(y_iv)
    n_cf = length(y_cf_raw)

    y_delta, sort_perm = if delta_method === :rank_aligned
        _rank_aligned_delta(y_iv, y_cf_raw)
    elseif delta_method === :mean_difference
        y_iv .- mean(y_cf_raw), collect(1:n_iv)
    else
        throw(
            ArgumentError(
                "Unknown delta_method :$(delta_method); " *
                "choose :rank_aligned or :mean_difference."
            )
        )
    end

    if all(y_delta .== 0.0)
        @warn "All y_delta values are zero; intervention has no measurable effect. " *
            "RSA will return sentinel values."
    end

    fs_all = fs === nothing ? feature_set(rs) : fs
    if nrow(fs_all) != length(cf_mask)
        throw(
            ArgumentError(
                "Precomputed `fs` has $(nrow(fs_all)) rows; expected $(length(cf_mask)) " *
                "(length of cf_mask)."
            )
        )
    end
    fs_iv = fs_all[iv_mask, :]
    drop_cols = [c for c in names(fs_iv) if Symbol(c) in DHW_STAT_COLS]
    fs_iv = fs_iv[sort_perm, Not(drop_cols)]

    ci = if bootstrap_n > 0
        bs = bootstrap(mean, y_delta, BasicSampling(bootstrap_n))
        _, lo, hi = confint(bs, BasicConfInt(0.95))[1]
        (lo, hi)
    else
        nothing
    end

    return (; y_delta=y_delta, fs_intervention=fs_iv, iv_perm=sort_perm, bootstrap_ci=ci)
end
