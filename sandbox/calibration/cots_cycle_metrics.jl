using Statistics

# Cycle-aware validation metrics for COTS calibration.
#
# These functions are intentionally dependency-light so they can move into the
# standalone calibration study repo. The objective is a loss: lower is better.

struct CotsPeak
    year::Int
    value::Float64
    prominence::Float64
end

struct LagCorrelation
    pearson::Float64
    spearman::Float64
    lag::Int
    n::Int
end

struct CotsCycleScore
    total_loss::Float64
    matched_rmse::Float64
    percent_bias_abs::Float64
    peak_count_penalty::Float64
    peak_timing_penalty::Float64
    period_penalty::Float64
    amplitude_penalty::Float64
    flatline_penalty::Float64
    lag_correlation_penalty::Float64
    best_lag_pearson::Float64
    best_lag_spearman::Float64
    best_lag_years::Int
    n_obs::Int
    n_sim_peaks::Int
    n_obs_peaks::Int
    sim_peak_years::Vector{Int}
    obs_peak_years::Vector{Int}
end

function _finite_nonnegative(v::AbstractVector{<:Real})::Vector{Float64}
    return [isfinite(Float64(x)) ? max(0.0, Float64(x)) : 0.0 for x in v]
end

function normalize_series(v::AbstractVector{<:Real})::Vector{Float64}
    x = _finite_nonnegative(v)
    isempty(x) && return Float64[]
    mx = maximum(x)
    mx <= 0.0 && return zeros(Float64, length(x))
    return x ./ mx
end

function moving_average(v::AbstractVector{<:Real}, window::Int=3)::Vector{Float64}
    x = _finite_nonnegative(v)
    n = length(x)
    n == 0 && return Float64[]
    window <= 1 && return x
    radius = max(0, div(window, 2))
    out = similar(x)
    for i in 1:n
        lo = max(1, i - radius)
        hi = min(n, i + radius)
        out[i] = mean(@view x[lo:hi])
    end
    return out
end

function rank_average(v::AbstractVector{<:Real})::Vector{Float64}
    order = sortperm(v)
    ranks = zeros(Float64, length(v))
    i = 1
    while i <= length(v)
        j = i
        while j < length(v) && v[order[j + 1]] == v[order[i]]
            j += 1
        end
        avg_rank = (i + j) / 2
        for k in i:j
            ranks[order[k]] = avg_rank
        end
        i = j + 1
    end
    return ranks
end

function safe_cor(x::AbstractVector{<:Real}, y::AbstractVector{<:Real})::Float64
    length(x) < 2 && return NaN
    std(x) == 0.0 && return NaN
    std(y) == 0.0 && return NaN
    return cor(Float64.(x), Float64.(y))
end

function detect_cots_peaks(
    years::AbstractVector{<:Integer},
    values::AbstractVector{<:Real};
    min_height::Float64=0.25,
    min_prominence::Float64=0.12,
    min_separation::Int=4,
    smooth_window::Int=3
)::Vector{CotsPeak}
    @assert length(years) == length(values) "years and values must have same length"
    n = length(values)
    n < 3 && return CotsPeak[]

    y = moving_average(normalize_series(values), smooth_window)
    candidates = CotsPeak[]
    for i in 2:(n - 1)
        if y[i] >= y[i - 1] && y[i] > y[i + 1] && y[i] >= min_height
            left_min = minimum(@view y[1:i])
            right_min = minimum(@view y[i:n])
            prominence = y[i] - max(left_min, right_min)
            if prominence >= min_prominence
                push!(candidates, CotsPeak(Int(years[i]), y[i], prominence))
            end
        end
    end

    sorted = sort(candidates; by=p -> (p.value, p.prominence), rev=true)
    kept = CotsPeak[]
    for p in sorted
        if all(abs(p.year - q.year) >= min_separation for q in kept)
            push!(kept, p)
        end
    end
    return sort(kept; by=p -> p.year)
end

function _nearest_year_losses(sim_years::Vector{Int}, obs_years::Vector{Int}; tolerance::Int=6)::Tuple{Float64,Int}
    isempty(obs_years) && return (0.0, 0)
    isempty(sim_years) && return (Float64(tolerance), 0)

    total = 0.0
    matched = 0
    used = falses(length(sim_years))
    for oy in obs_years
        best_idx = 0
        best_dist = typemax(Int)
        for (i, sy) in pairs(sim_years)
            used[i] && continue
            d = abs(sy - oy)
            if d < best_dist
                best_dist = d
                best_idx = i
            end
        end
        if best_idx > 0 && best_dist <= tolerance
            used[best_idx] = true
            matched += 1
            total += best_dist / tolerance
        else
            total += 1.0
        end
    end
    return (total / length(obs_years), matched)
end

function _period_penalty(years::Vector{Int}; expected_period::Float64=15.0, tolerance::Float64=6.0)::Float64
    length(years) < 2 && return 1.0
    periods = diff(years)
    isempty(periods) && return 1.0
    return mean(min(abs(Float64(p) - expected_period) / tolerance, 1.0) for p in periods)
end

function _matched_series(
    sim_years::AbstractVector{<:Integer},
    sim::AbstractVector{<:Real},
    obs_years::AbstractVector{<:Integer},
    obs::AbstractVector{<:Real}
)::Tuple{Vector{Float64},Vector{Float64}}
    sim_lookup = Dict(Int(y) => Float64(v) for (y, v) in zip(sim_years, sim))
    xs = Float64[]
    ys = Float64[]
    for (y, v) in zip(obs_years, obs)
        key = Int(y)
        if haskey(sim_lookup, key)
            push!(xs, sim_lookup[key])
            push!(ys, Float64(v))
        end
    end
    return xs, ys
end

function lagged_correlation(
    sim_years::AbstractVector{<:Integer},
    sim_values::AbstractVector{<:Real},
    obs_years::AbstractVector{<:Integer},
    obs_values::AbstractVector{<:Real};
    max_lag::Int=6,
    min_overlap::Int=4
)::LagCorrelation
    best = LagCorrelation(NaN, NaN, 0, 0)
    best_score = -Inf
    sim_norm = normalize_series(sim_values)
    obs_norm = normalize_series(obs_values)

    for lag in -max_lag:max_lag
        # Positive lag means the simulated trajectory is shifted later before
        # comparison. A best positive lag therefore suggests the raw simulation
        # is too early relative to observations.
        shifted_years = Int.(sim_years) .+ lag
        x, y = _matched_series(shifted_years, sim_norm, obs_years, obs_norm)
        length(x) < min_overlap && continue
        r = safe_cor(x, y)
        rho = safe_cor(rank_average(x), rank_average(y))
        score = coalesce(isnan(r) ? missing : r, -1.0) + 0.5 * coalesce(isnan(rho) ? missing : rho, -1.0)
        score -= 0.02 * abs(lag)
        if score > best_score
            best = LagCorrelation(r, rho, lag, length(x))
            best_score = score
        end
    end
    return best
end

function cots_cycle_score(
    sim_years::AbstractVector{<:Integer},
    sim_values::AbstractVector{<:Real},
    obs_years::AbstractVector{<:Integer},
    obs_values::AbstractVector{<:Real};
    expected_period::Float64=15.0,
    peak_tolerance::Int=6,
    max_correlation_lag::Int=6,
    min_peak_height::Float64=0.25,
    min_peak_prominence::Float64=0.12,
    rmse_weight::Float64=0.35,
    bias_weight::Float64=0.002,
    peak_count_weight::Float64=1.0,
    peak_timing_weight::Float64=1.5,
    period_weight::Float64=0.75,
    amplitude_weight::Float64=0.75,
    flatline_weight::Float64=2.0,
    lag_correlation_weight::Float64=0.6
)::CotsCycleScore
    @assert length(sim_years) == length(sim_values) "sim years and values must have same length"
    @assert length(obs_years) == length(obs_values) "obs years and values must have same length"

    sim_norm = normalize_series(sim_values)
    obs_norm = normalize_series(obs_values)
    sim_peaks = detect_cots_peaks(sim_years, sim_norm; min_height=min_peak_height, min_prominence=min_peak_prominence)
    obs_peaks = detect_cots_peaks(obs_years, obs_norm; min_height=min_peak_height, min_prominence=min_peak_prominence)
    sim_peak_years = [p.year for p in sim_peaks]
    obs_peak_years = [p.year for p in obs_peaks]

    sim_matched, obs_matched = _matched_series(sim_years, sim_norm, obs_years, obs_norm)
    n_obs = length(obs_matched)
    matched_rmse = n_obs == 0 ? 1.0 : sqrt(mean((sim_matched .- obs_matched) .^ 2))
    percent_bias_abs = if n_obs == 0 || sum(obs_matched) == 0.0
        100.0
    else
        abs(100.0 * (sum(sim_matched) - sum(obs_matched)) / sum(obs_matched))
    end

    n_obs_peaks = length(obs_peaks)
    n_sim_peaks = length(sim_peaks)
    peak_count_penalty = n_obs_peaks == 0 ? (n_sim_peaks == 0 ? 0.0 : 1.0) : abs(n_sim_peaks - n_obs_peaks) / max(n_obs_peaks, 1)
    peak_count_penalty = min(peak_count_penalty, 1.0)

    peak_timing_penalty, matched_peaks = _nearest_year_losses(sim_peak_years, obs_peak_years; tolerance=peak_tolerance)
    period_penalty = n_obs_peaks >= 2 ? _period_penalty(sim_peak_years; expected_period=expected_period) : 0.0

    sim_amp = isempty(sim_norm) ? 0.0 : maximum(sim_norm) - minimum(sim_norm)
    obs_amp = isempty(obs_norm) ? 0.0 : maximum(obs_norm) - minimum(obs_norm)
    amplitude_penalty = obs_amp <= 0.0 ? 0.0 : min(abs(sim_amp - obs_amp) / obs_amp, 1.0)

    flatline_penalty = sim_amp < 0.15 ? 1.0 : 0.0
    if n_obs_peaks > 0 && matched_peaks == 0
        flatline_penalty = max(flatline_penalty, 0.5)
    end

    lag_corr = lagged_correlation(sim_years, sim_norm, obs_years, obs_norm; max_lag=max_correlation_lag)
    lag_r = isnan(lag_corr.pearson) ? 0.0 : lag_corr.pearson
    lag_rho = isnan(lag_corr.spearman) ? 0.0 : lag_corr.spearman
    lag_correlation_penalty = 1.0 - clamp(0.7 * max(lag_r, 0.0) + 0.3 * max(lag_rho, 0.0), 0.0, 1.0)

    total_loss =
        rmse_weight * matched_rmse +
        bias_weight * percent_bias_abs +
        peak_count_weight * peak_count_penalty +
        peak_timing_weight * peak_timing_penalty +
        period_weight * period_penalty +
        amplitude_weight * amplitude_penalty +
        flatline_weight * flatline_penalty +
        lag_correlation_weight * lag_correlation_penalty

    return CotsCycleScore(
        total_loss,
        matched_rmse,
        percent_bias_abs,
        peak_count_penalty,
        peak_timing_penalty,
        period_penalty,
        amplitude_penalty,
        flatline_penalty,
        lag_correlation_penalty,
        lag_corr.pearson,
        lag_corr.spearman,
        lag_corr.lag,
        n_obs,
        n_sim_peaks,
        n_obs_peaks,
        sim_peak_years,
        obs_peak_years
    )
end