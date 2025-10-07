using FLoops, DataStructures

"""
    reef_condition_index(rc::AbstractArray, sv::AbstractArray, juves::AbstractArray,)::AbstractArray
    reef_condition_index(rc::AbstractArray, juves::AbstractArray, sv::AbstractArray, rubble::AbstractArray)::AbstractArray
    reef_condition_index(rs::ResultSet)::AbstractArray{<:Real}
    reef_condition_index(rs::ResultSet, rubble::AbstractArray)::AbstractArray{<:Real}

Estimates a Reef Condition Index (RCI) using either the 3-metric version using relative
cover, juveniles, shelter volume or the 4-metric versions with rubble added.

The RCI is a single value that indicates the condition of a reef.

# Notes
Juveniles are made relative to maximum observed juvenile density (51.8/m²)
See notes for `juvenile_indicator()`

# Arguments
- `rc` : Relative coral cover across all groups
- `juves` : Abundance of coral juveniles < 5 cm diameter
- `sv` : Shelter volume based on coral sizes and abundances
- `rubble` : Cover of rubble (optional)
- `rs` : A ResultSet object.

# Returns
YAXArray[timesteps ⋅ locations ⋅ scenarios]
"""
function _reef_condition_index(
    rc::AbstractArray{<:Real,3},
    sv::AbstractArray{<:Real,3},
    juves::AbstractArray{<:Real,3}
)::AbstractArray{<:Real}
    out_rrci = zeros(eltype(rc), size(rc)...)

    juves_rel_baseline = juves ./ maximum(juves)
    for scen_idx in axes(rc, axis_index(rc, :scenarios))
        @views ADRIAIndicators.reef_condition_index!(
            rc.data[:, :, scen_idx],
            sv.data[:, :, scen_idx],
            juves_rel_baseline.data[:, :, scen_idx],
            out_rrci[:, :, scen_idx]
        )
    end

    return DataCube(out_rrci, (:timesteps, :locations, :scenarios))
end
function _reef_condition_index(rs::ResultSet)::AbstractArray{<:Real}
    return _reef_condition_index(
        relative_cover(rs),
        relative_shelter_volume(rs),
        juvenile_indicator(rs)
    )
end
function _reef_condition_index(
    rc::AbstractArray{<:Real,3},
    sv::AbstractArray{<:Real,3},
    juves::AbstractArray{<:Real,3},
    rubble::AbstractArray{<:Real,3}
)::AbstractArray{<:Real}
    juves = collect(juves ./ maximum(juves))

    out_rci = zeros(eltype(rc), size(rc)...)
    for scen_idx in axes(rc, axis_index(rc, :scenarios))
        @views ADRIAIndicators.reef_condition_index!(
            rc.data[:, :, scen_idx],
            sv.data[:, :, scen_idx],
            juves[:, :, scen_idx],
            rubble.data[:, :, scen_idx],
            out_rci[:, :, scen_idx];
        )
    end

    return DataCube(out_rci, (:timesteps, :locations, :scenarios))
end
function _reef_condition_index(
    rs::ResultSet, rubble::YAXArray
)::AbstractArray{<:Real}
    return _reef_condition_index(
        relative_cover(rs),
        relative_shelter_volume(rs),
        juvenile_indicator(rs),
        rubble
    )
end
reef_condition_index = Metric(
    _reef_condition_index,
    (:timesteps, :locations, :scenarios),
    (:timesteps, :locations, :scenarios),
    "RCI",
    IS_NOT_RELATIVE
)

"""
    scenario_rci(rci::YAXArray, tac::YAXArray; kwargs...)
    scenario_rci(rci::YAXArray, rubble::YAXArray; kwargs...)
    scenario_rci(rs::ResultSet; kwargs...)

Extract the total populated area of locations with Reef Condition Index of "Good" or higher
for each scenario for the entire domain.
"""
function _scenario_rci(rci::YAXArray, tac::YAXArray; kwargs...)
    rci_sliced = slice_results(rci; kwargs...)
    tac_sliced = slice_results(tac; kwargs...)

    # We want sum of populated area >= "good" condition
    return scenario_trajectory(tac_sliced .* (rci_sliced .> 0.35); metric=sum)
end
function _scenario_rci(rs::ResultSet, rubble::YAXArray; kwargs...)
    rci = reef_condition_index(rs, rubble)
    tac = total_absolute_cover(rs)

    return _scenario_rci(rci, tac; kwargs...)
end
function _scenario_rci(rs::ResultSet; kwargs...)
    rci = reef_condition_index(rs)
    tac = total_absolute_cover(rs)

    return _scenario_rci(rci, tac; kwargs...)
end
scenario_rci = Metric(
    _scenario_rci,
    (:timesteps, :locations, :scenarios),
    (:timesteps, :scenarios),
    "RCI",
    IS_NOT_RELATIVE
)

"""
    reef_tourism_index(rc::AbstractArray{<:Real,3}, sv::AbstractArray{<:Real,3}, juves::AbstractArray{<:Real,3}, cots::AbstractArray{<:Real,3}, rubble::AbstractArray{<:Real,3})::AbstractArray
    reef_tourism_index(rc::AbstractArray{<:Real,3}, ce::AbstractArray{<:Real,3}, sv::AbstractArray{<:Real,3}, juves::AbstractArray{<:Real,3})::AbstractArray
    reef_tourism_index(rs::ResultSet, cots::YAXArray, rubble::YAXArray)::AbstractArray
    reef_tourism_index(rs::ResultSet)::AbstractArray

Estimate tourism index.

This metric is a variation of the Reef Condition Index, but weighted by metrics known to be
of importance to tourists. This version uses 5 metrics: relative cover, shelter volume,
juvenile abundance, CoTS, and rubble.

# Arguments
- `rs` : ResultSet
- `cots` : Outbreak status of Crown-of-Thorns Starfish
- `rubble` : Cover of rubble

"""
function _reef_tourism_index(
    rc::AbstractArray{<:Real,3},
    ce::AbstractArray{<:Real,3},
    sv::AbstractArray{<:Real,3},
    juves::AbstractArray{<:Real,3}
)::AbstractArray
    out_rrti = zeros(eltype(rc), size(rc)...)
    for scen_idx in axes(rc, axis_index(rc, :scenarios))
        @views ADRIAIndicators.reef_tourism_index!(
            rc.data[:, :, scen_idx],
            ce.data[:, :, scen_idx],
            sv.data[:, :, scen_idx],
            juves.data[:, :, scen_idx],
            out_rrti[:, :, scen_idx]
        )
    end

    return DataCube(out_rrti, (:timesteps, :locations, :scenarios))
end
function _reef_tourism_index(rs::ResultSet)::AbstractArray
    rc = relative_cover(rs)
    ce = coral_evenness(rs)
    sv = relative_shelter_volume(rs)
    juves = juvenile_indicator(rs)

    return _reef_tourism_index(rc, ce, sv, juves)
end
function _reef_tourism_index(
    rc::AbstractArray{<:Real,3},
    sv::AbstractArray{<:Real,3},
    juves::AbstractArray{<:Real,3},
    cots::AbstractArray{<:Real,3},
    rubble::AbstractArray{<:Real,3}
)::AbstractArray
    out_rti = zeros(eltype(rc), size(rc)...)
    for scen_idx in axes(rc, axis_index(rc, :scenarios))
        @views ADRIAIndicators.reef_tourism_index!(
            rc.data[:, :, scen_idx],
            sv.data[:, :, scen_idx],
            juves.data[:, :, scen_idx],
            cots.data[:, :, scen_idx],
            rubble.data[:, :, scen_idx],
            out_rti[:, :, scen_idx]
        )
    end

    return DataCube(out_rti, (:timesteps, :locations, :scenarios))
end
function _reef_tourism_index(rs::ResultSet, cots::YAXArray, rubble::YAXArray)::AbstractArray
    rc = relative_cover(rs)
    sv = relative_shelter_volume(rs)
    juves = juvenile_indicator(rs)

    return _reef_tourism_index(rc, sv, juves, cots, rubble)
end
reef_tourism_index = Metric(
    _reef_tourism_index,
    (:timesteps, :locations, :scenarios),
    (:timesteps, :locations, :scenarios),
    "RTI",
    IS_NOT_RELATIVE
)

"""
    scenario_rti(rs::ResultSet, cots::YAXArray, rubble::YAXArray; kwargs...)
    scenario_rti(rs::ResultSet; kwargs)

Calculate the mean Reef Tourism Index (RTI) for each scenario for the entire domain.
"""
function _scenario_rti(rti::YAXArray; kwargs...)
    rti_sliced = slice_results(rti; kwargs...)
    return scenario_trajectory(rti_sliced)
end
function _scenario_rti(rs::ResultSet; kwargs...)
    rti = reef_tourism_index(rs)
    return _scenario_rti(rti; kwargs...)
end
function _scenario_rti(rs::ResultSet, cots::YAXArray, rubble::YAXArray; kwargs...)
    return _scenario_rti(reef_tourism_index(rs, cots, rubble); kwargs...)
end
scenario_rti = Metric(
    _scenario_rti,
    (:timesteps, :locations, :scenarios),
    (:timesteps, :scenarios),
    "RTI",
    IS_NOT_RELATIVE
)

"""
    reef_fish_index(rc::AbstractArray)
    reef_fish_index(rs::ResultSet)

The Reef Fish Index (RFI) estimates fish biomass from relative coral cover.

A linear regression (developed by Dr. R. Heneghan, Queensland University of Technology)
is used to indicate the relationship between coral cover and fish biomass.
The regression was developed with digitized data from Figures 4a and 6b in
Graham & Nash (2013; see [1]).

Values are provided ∈ [0, 1], where 1 indicates maximum fish biomass.

Note: Coral cover here is relative to coral habitable area (\$k\$ area).

# Arguments
- `rc` : Relative cover

# Returns
YAXArray[timesteps ⋅ locations ⋅ scenarios], values in kg/km²

# References
1. Graham, N.A.J., Nash, K.L., 2013.
The importance of structural complexity in coral reef ecosystems.
Coral Reefs 32, 315–326.
https://doi.org/10.1007/s00338-012-0984-y
"""
function _reef_fish_index(rc::AbstractArray{T})::AbstractArray{T} where {T<:Real}
    n_scenarios = size(rc, axis_index(rc, :scenarios))
    out_rfi = zeros(eltype(rc), size(rc)...)

    for scen_idx in 1:n_scenarios
        @views ADRIAIndicators.reef_fish_index!(
            rc.data[:, :, scen_idx],
            out_rfi[:, :, scen_idx]
        )
    end

    dims = (
        rc.timesteps,
        rc.locations,
        rc.scenarios
    )

    return YAXArray(dims, out_rfi, rc.properties)
end
function _reef_fish_index(rs::ResultSet)
    rc = relative_cover(rs)

    return _reef_fish_index(rc)
end
reef_fish_index = Metric(
    _reef_fish_index,
    (:timesteps, :locations, :scenarios),
    (:timesteps, :locations, :scenarios),
    "RFI",
    IS_NOT_RELATIVE
)

"""
scenario_rfi(rfi::YAXArray; kwargs...)
scenario_rfi(rs::ResultSet; kwargs...)

Calculate the mean Reef Fish Index (RFI) for each scenario for the entire domain.
"""
function _scenario_rfi(rfi::YAXArray; kwargs...)
    rfi_sliced = slice_results(rfi; kwargs...)
    return scenario_trajectory(rfi_sliced)
end
function _scenario_rfi(rs::ResultSet; kwargs...)
    return _scenario_rfi(reef_fish_index(rs); kwargs...)
end
scenario_rfi = Metric(
    _scenario_rfi,
    (:timesteps, :locations, :scenarios),
    (:timesteps, :scenarios),
    "RFI",
    IS_NOT_RELATIVE
)

"""
    _reef_biodiviersity_condition_index(rs::ResultSet)
"""
function _reef_biodiversity_condition_index(
    rc::AbstractArray{<:Real,3},
    sv::AbstractArray{<:Real,3},
    ce::AbstractArray{<:Real,3}
)::AbstractArray{<:Real,3}
    cd = collect(1 .- 1 ./ ce)
    cd = replace!(
        cd, NaN => 0.0, Inf => 0.0
    )

    rbci::Array = zeros(eltype(rc), size(rc)...)
    @views for scen_idx in axes(rc, axis_index(rc, :scenarios))
        ADRIAIndicators.reef_biodiversity_condition_index!(
            rc[:, :, scen_idx],
            cd[:, :, scen_idx],
            sv[:, :, scen_idx],
            rbci[:, :, scen_idx]
        )
    end
    dims = (
        :timesteps,
        :locations,
        :scenarios
    )

    return DataCube(rbci, dims, rc.properties)
end
function _reef_biodiversity_condition_index(rs::ResultSet)
    rc = relative_cover(rs)
    sv = relative_shelter_volume(rs)
    ce = coral_evenness(rs)

    return _reef_biodiversity_condition_index(rc, sv, ce)
end
reef_biodiversity_condition_index = Metric(
    _reef_biodiversity_condition_index,
    (:timesteps, :locations, :scenarios),
    (:timesteps, :locations, :scenarios),
    "RBCI",
    IS_NOT_RELATIVE
)

function _scenario_rbci(rbci::YAXArray; kwargs...)
    rbci_sliced = slice_results(rbci; kwargs...)
    return scenario_trajectory(rbci_sliced)
end
function _scenario_rbci(rs::ResultSet; kwargs...)
    rbci = _reef_biodiversity_condition_index(rs)
    return _scenario_rbci(rbci)
end
scenario_rbci = Metric(
    _scenario_rbci,
    (:timesteps, :locations, :scenarios),
    (:timesteps, :scenarios),
    "RBCI",
    IS_NOT_RELATIVE
)
