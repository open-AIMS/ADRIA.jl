"""Frozen score v2 (crash_mechanism_protocol.md, "Preregistered frozen score v2", 2026-10-14).

A wrapper around the unchanged v1 `cots_cycle_score` (cots_cycle_metrics.jl):

1. per-reef observation treatment: North Direction is scored against the raw observed
   series (its 1995 peak); the other reefs keep `smoothed_3y`;
2. Eyrie's first wave is unobserved (triennial surveys): model peaks earlier than
   (earliest observed peak - peak tolerance) are dropped from the peak count only.

Pure functions; the caller supplies the v1 score.
"""

const SCORE_V2_TREATMENT = Dict(
    "Lizard Island Reef" => "smoothed_3y",
    "MacGillivray Reef" => "smoothed_3y",
    "North Direction Reef" => "raw",
    "Eyrie Reef" => "smoothed_3y",
)
const SCORE_V2_UNOBSERVED_FIRST_WAVE = ("Eyrie Reef",)
const SCORE_V2_PEAK_TOLERANCE = 6         # v1 default `peak_tolerance`
const SCORE_V2_PEAK_COUNT_WEIGHT = 1.25   # v1 default `peak_count_weight`

score_v2_primary(reef::AbstractString, treatment::AbstractString) =
    get(SCORE_V2_TREATMENT, reef, "smoothed_3y") == treatment

"""Peak-count penalty with model peaks before the observable window dropped."""
function score_v2_peak_count_penalty(sim_peak_years::AbstractVector{<:Integer},
        obs_peak_years::AbstractVector{<:Integer}; tolerance::Int=SCORE_V2_PEAK_TOLERANCE)
    n_obs = length(obs_peak_years)
    cutoff = n_obs == 0 ? typemin(Int) : minimum(obs_peak_years) - tolerance
    n_sim = count(>=(cutoff), sim_peak_years)
    n_obs == 0 && return (n_sim == 0 ? 0.0 : 1.0), n_sim
    return min(abs(n_sim - n_obs) / n_obs, 1.0), n_sim
end

"""v2 loss for one reef-treatment row, from the v1 loss and its peak-count terms."""
function score_v2_loss(reef::AbstractString, loss_v1::Real, peak_count_penalty_v1::Real,
        sim_peak_years::AbstractVector{<:Integer}, obs_peak_years::AbstractVector{<:Integer})
    reef in SCORE_V2_UNOBSERVED_FIRST_WAVE || return Float64(loss_v1)
    penalty_v2, _ = score_v2_peak_count_penalty(sim_peak_years, obs_peak_years)
    return Float64(loss_v1) - SCORE_V2_PEAK_COUNT_WEIGHT * (peak_count_penalty_v1 - penalty_v2)
end
