"""Frozen score v3 (crash_mechanism_protocol.md, "Preregistered frozen score v3", 2026-10-14).

v3 = score v2 on 1992-2021 (manta) + 1.5 x the diver-anchored 2020s error per reef:

    D_reef = mean over target years of |log((A_model + 0.5) / (T + 0.5))|

with A_model adults/ha and T the diver-anchored target (build_2020s_density_targets.py).
The manta cut is 2021, not the preregistered 2019 (correction, 2026-10-14): Eyrie has no
manta survey between 2015 and 2021, so a 2019 cut shortened its scoring window to 2015 and
dropped its 2014 model peak. 2021 is the first post-2019 survey year on every scored reef.
Diver targets are scored from 2020. Pure functions; the caller supplies the v2 part.
"""

const SCORE_V3_LAST_MANTA_YEAR = 2021
const SCORE_V3_FIRST_DIVER_YEAR = 2020
const SCORE_V3_DIVER_WEIGHT = 1.5
const SCORE_V3_OFFSET_HA = 0.5

"""Mean absolute log error of modelled adults against diver targets over the target years
present in the model; returns (error, n). n = 0 gives NaN."""
function diver_log_error(model_years::AbstractVector{<:Integer}, adults::AbstractVector{<:Real},
        target_years::AbstractVector{<:Integer}, targets::AbstractVector{<:Real};
        offset::Float64=SCORE_V3_OFFSET_HA)
    lookup = Dict(Int(y) => Float64(a) for (y, a) in zip(model_years, adults))
    errors = Float64[]
    for (year, target) in zip(target_years, targets)
        haskey(lookup, Int(year)) || continue
        isfinite(target) && target >= 0.0 || throw(ArgumentError("diver targets must be finite and non-negative"))
        push!(errors, abs(log((lookup[Int(year)] + offset) / (Float64(target) + offset))))
    end
    return isempty(errors) ? (NaN, 0) : (sum(errors) / length(errors), length(errors))
end

score_v3_loss(loss_v2_to2019::Real, diver_error::Real) =
    Float64(loss_v2_to2019) + SCORE_V3_DIVER_WEIGHT * Float64(diver_error)
