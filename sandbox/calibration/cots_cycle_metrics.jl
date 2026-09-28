using Statistics

# Calendar-aware validation metrics for COTS calibration. ADRIA adult density
# and AIMS COTS-per-tow are linked with a fitted per-reef observation scale.

struct CotsPeak
    year::Int
    value::Float64
    prominence::Float64
    normalized_value::Float64
    normalized_prominence::Float64
end

# Calendar-year helpers follow.

function yearly_series(years::AbstractVector{<:Integer}, values::AbstractVector{<:Real})
    @assert length(years) == length(values) "years and values must have same length"
    buckets = Dict{Int,Vector{Float64}}()
    for (year, value) in zip(years, _finite_nonnegative(values))
        push!(get!(buckets, Int(year), Float64[]), value)
    end
    sorted_years = sort!(collect(keys(buckets)))
    return sorted_years, [mean(buckets[year]) for year in sorted_years]
end

function calendar_moving_average(
    years::AbstractVector{<:Integer}, values::AbstractVector{<:Real}; window_years::Int=3
)::Vector{Float64}
    @assert length(years) == length(values)
    x = _finite_nonnegative(values)
    isempty(x) && return Float64[]
    window_years <= 1 && return x
    radius = max(0, div(window_years, 2))
    integer_years = Int.(years)
    return [mean(x[abs.(integer_years .- year) .<= radius]) for year in integer_years]
end

function moving_average(v::AbstractVector{<:Real}, window::Int=3)::Vector{Float64}
    x = _finite_nonnegative(v)
    n = length(x)
    n == 0 && return Float64[]
    window <= 1 && return x
    radius = max(0, div(window, 2))
    return [mean(@view x[max(1, i - radius):min(n, i + radius)]) for i in 1:n]
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
    years::AbstractVector{<:Integer}, values::AbstractVector{<:Real};
    min_height::Float64=0.25,
    min_prominence::Float64=0.12,
    min_separation::Int=10,
    smooth_window_years::Int=1,
    prominence_window_years::Int=10,
    max_neighbor_gap_years::Int=6
)::Vector{CotsPeak}
    sorted_years, raw = yearly_series(years, values)
    n = length(raw)
    n < 3 && return CotsPeak[]
    smoothed = calendar_moving_average(sorted_years, raw; window_years=smooth_window_years)
    normalized = normalize_series(smoothed)
    candidates = CotsPeak[]

    # Boundary values lack observations on both sides and are not confirmed peaks.
    for i in 2:(n - 1)
        left_gap = sorted_years[i] - sorted_years[i - 1]
        right_gap = sorted_years[i + 1] - sorted_years[i]
        (left_gap > max_neighbor_gap_years || right_gap > max_neighbor_gap_years) && continue
        is_local_peak = normalized[i] >= normalized[i - 1] && normalized[i] > normalized[i + 1]
        (!is_local_peak || normalized[i] < min_height) && continue

        left_idx = findall(y -> sorted_years[i] - prominence_window_years <= y < sorted_years[i], sorted_years)
        right_idx = findall(y -> sorted_years[i] < y <= sorted_years[i] + prominence_window_years, sorted_years)
        (isempty(left_idx) || isempty(right_idx)) && continue
        norm_prominence = normalized[i] - max(minimum(normalized[left_idx]), minimum(normalized[right_idx]))
        norm_prominence < min_prominence && continue
        raw_prominence = smoothed[i] - max(minimum(smoothed[left_idx]), minimum(smoothed[right_idx]))
        push!(candidates, CotsPeak(sorted_years[i], smoothed[i], max(0.0, raw_prominence), normalized[i], norm_prominence))
    end

    sorted = sort(candidates; by=p -> (p.normalized_value, p.normalized_prominence), rev=true)
    kept = CotsPeak[]
    for peak in sorted
        all(abs(peak.year - other.year) >= min_separation for other in kept) && push!(kept, peak)
    end
    return sort(kept; by=p -> p.year)
end

function _matched_series(sim_years, sim, obs_years, obs)
    sim_lookup = Dict(Int(y) => Float64(v) for (y, v) in zip(sim_years, sim))
    xs = Float64[]
    ys = Float64[]
    for (year, value) in zip(obs_years, obs)
        if haskey(sim_lookup, Int(year))
            push!(xs, sim_lookup[Int(year)])
            push!(ys, Float64(value))
        end
    end
    return xs, ys
end

function _common_window(sim_years, sim_values, obs_years, obs_values)
    if isempty(sim_years) || isempty(obs_years)
        return Int[], Float64[], Int[], Float64[], 0, -1
    end
    start_year = max(first(sim_years), first(obs_years))
    end_year = min(last(sim_years), last(obs_years))
    end_year < start_year && return Int[], Float64[], Int[], Float64[], 0, -1
    sim_keep = (sim_years .>= start_year) .& (sim_years .<= end_year)
    obs_keep = (obs_years .>= start_year) .& (obs_years .<= end_year)
    return sim_years[sim_keep], sim_values[sim_keep], obs_years[obs_keep], obs_values[obs_keep], start_year, end_year
end

function fit_observation_scale(
    sim::AbstractVector{<:Real}, obs::AbstractVector{<:Real};
    mode::Symbol=:least_squares, fixed_scale::Union{Nothing,Float64}=nothing
)::Float64
    fixed_scale !== nothing && return max(0.0, fixed_scale)
    x = _finite_nonnegative(sim)
    y = _finite_nonnegative(obs)
    isempty(x) && return 0.0
    if mode == :least_squares
        denominator = sum(abs2, x)
        return denominator <= eps(Float64) ? 0.0 : max(0.0, sum(x .* y) / denominator)
    elseif mode == :max
        return maximum(x) <= 0.0 ? 0.0 : maximum(y) / maximum(x)
    elseif mode == :none
        return 1.0
    end
    throw(ArgumentError(string(:unknown_observation_scale_mode, :_, mode)))
end

function _match_peaks(sim_peaks::Vector{CotsPeak}, obs_peaks::Vector{CotsPeak}; tolerance::Int=6)
    used_sim = falses(length(sim_peaks))
    matches = Tuple{Int,Int}[]
    for (obs_idx, obs_peak) in enumerate(obs_peaks)
        best_sim_idx = 0
        best_distance = typemax(Int)
        for (sim_idx, sim_peak) in enumerate(sim_peaks)
            used_sim[sim_idx] && continue
            distance = abs(sim_peak.year - obs_peak.year)
            if distance < best_distance
                best_distance = distance
                best_sim_idx = sim_idx
            end
        end
        if best_sim_idx > 0 && best_distance <= tolerance
            used_sim[best_sim_idx] = true
            push!(matches, (best_sim_idx, obs_idx))
        end
    end
    return matches
end

function _period_penalty(sim_peaks, obs_peaks, matches; tolerance::Float64=6.0)
    length(obs_peaks) < 2 && return 0.0
    length(matches) < 2 && return 1.0
    ordered = sort(matches; by=match -> obs_peaks[match[2]].year)
    penalties = Float64[]
    for i in 1:(length(ordered) - 1)
        sim_period = sim_peaks[ordered[i + 1][1]].year - sim_peaks[ordered[i][1]].year
        obs_period = obs_peaks[ordered[i + 1][2]].year - obs_peaks[ordered[i][2]].year
        push!(penalties, min(abs(Float64(sim_period - obs_period)) / tolerance, 1.0))
    end
    return isempty(penalties) ? 1.0 : mean(penalties)
end

function lagged_correlation(
    sim_years::AbstractVector{<:Integer}, sim_values::AbstractVector{<:Real},
    obs_years::AbstractVector{<:Integer}, obs_values::AbstractVector{<:Real};
    max_lag::Int=6, min_overlap::Int=4
)::LagCorrelation
    best = LagCorrelation(NaN, NaN, 0, 0)
    best_score = -Inf
    sim_norm = normalize_series(sim_values)
    obs_norm = normalize_series(obs_values)
    for lag in -max_lag:max_lag
        shifted_years = Int.(sim_years) .+ lag
        x, y = _matched_series(shifted_years, sim_norm, obs_years, obs_norm)
        length(x) < min_overlap && continue
        r = safe_cor(x, y)
        rho = safe_cor(rank_average(x), rank_average(y))
        score = (isnan(r) ? -1.0 : r) + 0.5 * (isnan(rho) ? -1.0 : rho) - 0.02 * abs(lag)
        if score > best_score
            best = LagCorrelation(r, rho, lag, length(x))
            best_score = score
        end
    end
    return best
end

function cots_cycle_score(
    sim_years::AbstractVector{<:Integer}, sim_values::AbstractVector{<:Real},
    obs_years::AbstractVector{<:Integer}, obs_values::AbstractVector{<:Real};
    peak_tolerance::Int=6,
    period_tolerance::Float64=6.0,
    max_correlation_lag::Int=6,
    min_peak_height::Float64=0.25,
    min_peak_prominence::Float64=0.12,
    observation_scale_mode::Symbol=:least_squares,
    fixed_observation_scale::Union{Nothing,Float64}=nothing,
    rmse_weight::Float64=0.25,
    bias_weight::Float64=0.001,
    peak_count_weight::Float64=1.25,
    peak_timing_weight::Float64=2.0,
    period_weight::Float64=0.75,
    amplitude_weight::Float64=2.0,
    flatline_weight::Float64=2.0,
    lag_correlation_weight::Float64=0.3
)::CotsCycleScore
    sim_yearly, sim_raw = yearly_series(sim_years, sim_values)
    obs_yearly, obs_raw = yearly_series(obs_years, obs_values)
    sim_yearly, sim_raw, obs_yearly, obs_raw, start_year, end_year = _common_window(
        sim_yearly, sim_raw, obs_yearly, obs_raw
    )

    sim_matched_raw, obs_matched = _matched_series(sim_yearly, sim_raw, obs_yearly, obs_raw)
    observation_scale = fit_observation_scale(
        sim_matched_raw, obs_matched;
        mode=observation_scale_mode, fixed_scale=fixed_observation_scale
    )
    sim_scaled = observation_scale .* sim_raw
    sim_matched, obs_matched = _matched_series(sim_yearly, sim_scaled, obs_yearly, obs_raw)
    sim_peaks = detect_cots_peaks(
        sim_yearly, sim_scaled;
        min_height=min_peak_height, min_prominence=min_peak_prominence,
        smooth_window_years=3
    )
    obs_peaks = detect_cots_peaks(
        obs_yearly, obs_raw;
        min_height=min_peak_height, min_prominence=min_peak_prominence,
        smooth_window_years=1
    )
    matches = _match_peaks(sim_peaks, obs_peaks; tolerance=peak_tolerance)

    n_obs = length(obs_matched)
    obs_reference = isempty(obs_raw) ? 1.0 : max(maximum(obs_raw), eps(Float64))
    matched_rmse = n_obs == 0 ? 1.0 : sqrt(mean((sim_matched .- obs_matched) .^ 2)) / obs_reference
    percent_bias_abs = if n_obs == 0 || sum(obs_matched) <= eps(Float64)
        100.0
    else
        abs(100.0 * (sum(sim_matched) - sum(obs_matched)) / sum(obs_matched))
    end

    n_obs_peaks = length(obs_peaks)
    n_sim_peaks = length(sim_peaks)
    peak_count_penalty = if n_obs_peaks == 0
        n_sim_peaks == 0 ? 0.0 : 1.0
    else
        min(abs(n_sim_peaks - n_obs_peaks) / n_obs_peaks, 1.0)
    end

    if n_obs_peaks == 0
        peak_timing_penalty = n_sim_peaks == 0 ? 0.0 : 1.0
        peak_height_penalty = n_sim_peaks == 0 ? 0.0 : 1.0
        peak_prominence_penalty = peak_height_penalty
    else
        timing_losses = fill(1.0, n_obs_peaks)
        height_losses = fill(1.0, n_obs_peaks)
        prominence_losses = fill(1.0, n_obs_peaks)
        for (sim_idx, obs_idx) in matches
            sim_peak = sim_peaks[sim_idx]
            obs_peak = obs_peaks[obs_idx]
            timing_losses[obs_idx] = abs(sim_peak.year - obs_peak.year) / peak_tolerance
            height_losses[obs_idx] = min(abs(sim_peak.value - obs_peak.value) / obs_reference, 2.0)
            prominence_reference = max(obs_peak.prominence, 0.1 * obs_reference)
            prominence_losses[obs_idx] = min(
                abs(sim_peak.prominence - obs_peak.prominence) / prominence_reference,
                2.0
            )
        end
        peak_timing_penalty = mean(timing_losses)
        peak_height_penalty = mean(height_losses)
        peak_prominence_penalty = mean(prominence_losses)
    end

    amplitude_penalty = 0.6 * peak_height_penalty + 0.4 * peak_prominence_penalty
    period_penalty = _period_penalty(sim_peaks, obs_peaks, matches; tolerance=period_tolerance)
    sim_relative_amplitude = if isempty(sim_scaled) || maximum(sim_scaled) <= 0.0
        0.0
    else
        (maximum(sim_scaled) - minimum(sim_scaled)) / maximum(sim_scaled)
    end
    flatline_penalty = sim_relative_amplitude < 0.15 ? 1.0 : 0.0
    n_obs_peaks > 0 && isempty(matches) && (flatline_penalty = max(flatline_penalty, 0.5))

    lag_corr = lagged_correlation(
        sim_yearly, sim_scaled, obs_yearly, obs_raw; max_lag=max_correlation_lag
    )
    lag_r = isnan(lag_corr.pearson) ? 0.0 : lag_corr.pearson
    lag_rho = isnan(lag_corr.spearman) ? 0.0 : lag_corr.spearman
    lag_correlation_penalty = 1.0 - clamp(
        0.7 * max(lag_r, 0.0) + 0.3 * max(lag_rho, 0.0), 0.0, 1.0
    )

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
        peak_height_penalty,
        peak_prominence_penalty,
        flatline_penalty,
        lag_correlation_penalty,
        lag_corr.pearson,
        lag_corr.spearman,
        lag_corr.lag,
        observation_scale,
        start_year,
        end_year,
        n_obs,
        n_sim_peaks,
        n_obs_peaks,
        length(matches),
        [peak.year for peak in sim_peaks],
        [peak.year for peak in obs_peaks]
    )
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
    peak_height_penalty::Float64
    peak_prominence_penalty::Float64
    flatline_penalty::Float64
    lag_correlation_penalty::Float64
    best_lag_pearson::Float64
    best_lag_spearman::Float64
    best_lag_years::Int
    observation_scale::Float64
    scoring_start_year::Int
    scoring_end_year::Int
    n_obs::Int
    n_sim_peaks::Int
    n_obs_peaks::Int
    n_matched_peaks::Int
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
