using FLoops, DataStructures

"""
    reef_condition_index(rc::AbstractArray, juves::AbstractArray, sv::AbstractArray)::AbstractArray
    reef_condition_index(rs::ResultSet)::AbstractArray{<:Real}

Estimates a Reef Condition Index (RCI) providing a single value that indicates the condition
of a reef across four metrics:
- coral cover
- abundance of juveniles
- shelter volume

# Notes
Juveniles are made relative to maximum observed juvenile density (51.8/m²)
See notes for `juvenile_indicator()`

# Arguments
- `rc` : Relative coral cover across all groups
- `evenness` : Evenness across all coral groups
- `sv` : Shelter volume based on coral sizes and abundances
- `juves` : Abundance of coral juveniles < 5 cm diameter

# Returns
YAXArray[timesteps ⋅ locations ⋅ scenarios]
"""
function _reef_condition_index(
    rc::AbstractArray{<:Real, 3},
    juves::AbstractArray{<:Real,3},
    sv::AbstractArray{<:Real,3}
)::AbstractArray{<:Real}
    juves = collect(juves ./ maximum(juves))

    n_timesteps, n_locations, _= size(rc)
    dummy_metric = ones(eltype(rc), n_timesteps, n_locations)

    out_rci = zeros(eltype(rc), size(rc)...)
    for scen_idx in axes(rc, axis_index(rc, :scenarios))
        @views ADRIAIndicators.reef_condition_index!(
            rc.data[:, :, scen_idx],
            sv.data[:, :, scen_idx],
            juves[:, :, scen_idx],
            dummy_metric,
            dummy_metric,
            out_rci[:, :, scen_idx];
            n_metrics=3
        )
    end

    return DataCube(out_rci, (:timesteps, :locations, :scenarios))
end
function _reef_condition_index(rs::ResultSet)::AbstractArray{<:Real}
    return _reef_condition_index(
        relative_cover(rs),
        juvenile_indicator(rs),
        relative_shelter_volume(rs)
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
    reef_tourism_index(rc::AbstractArray, evenness::AbstractArray, sv::AbstractArray, juves::AbstractArray, intcp_u::Vector)::AbstractArray
    reef_tourism_index(rs::ResultSet; intcp_u::Bool=false)::AbstractArray

Estimate tourism index.

Note:
This metric assumes all inputs (relative cover, evenness, shelter volume, coral juveniles)
are scaled between 0 and 1. For evenness, shelter volume and coral juveniles, a value of 1
may represent a theoretical maximum.

# Arguments
- `rc` : Relative coral cover across all groups
- `evenness` : Evenness across all coral groups
- `sv` : Shelter volume based on coral sizes and abundances
- `juves` : Abundance of coral juveniles < 5 cm diameter
- `intcp_u` : ?
"""
function _reef_tourism_index(
    rc::AbstractArray{<:Real,3},
    evenness::AbstractArray{<:Real,3},
    sv::AbstractArray{<:Real,3},
    juves::AbstractArray{<:Real,3},
    intcp_u::Vector{<:Real}
)::AbstractArray
    intcp = 0.47947 .+ intcp_u

    # TODO: Ryan to reinterpolate to account for no CoTS and no rubble
    # Apply unique intercepts for each scenario
    rti = cat(
        map(
            axe -> (
                intcp[axe] .+ (0.12764 .* rc[:, :, axe]) .+
                (0.31946 .* evenness[:, :, axe]) .+
                (0.11676 .* sv[:, :, axe]) .+
                (-0.0036065 .* juves[:, :, axe])
            ),
            axes(rc, 3)
        )...;
        dims=3
    )

    timesteps = rti.timesteps.val.data
    locations = rti.locations.val.data
    n_scenarios = size(rti, 3)
    rti = DataCube(
        rti.data; timesteps=timesteps, locations=locations, scenarios=1:n_scenarios
    )

    return round.(clamp.(rti, 0.1, 0.9), digits=2)
end
function _reef_tourism_index(rs::ResultSet; intcp_u::Bool=false)::AbstractArray
    rc::AbstractArray{<:Real} = relative_cover(rs)
    juves::AbstractArray{<:Real} = juvenile_indicator(rs)
    evenness::AbstractArray{<:Real} = coral_evenness(rs)
    sv::AbstractArray{<:Real} = relative_shelter_volume(rs)

    n_scens = size(rc, :scenarios)
    intcp = intcp_u ? rand(Normal(0.0, 0.163), n_scens) : zeros(n_scens)

    return _reef_tourism_index(rc, evenness, sv, juves, intcp)
end
reef_tourism_index = Metric(
    _reef_tourism_index,
    (:timesteps, :locations, :scenarios),
    (:timesteps, :locations, :scenarios),
    "RTI",
    IS_NOT_RELATIVE
)

"""
    scenario_rti(rti::YAXArray; kwargs...)
    scenario_rti(rs::ResultSet; kwargs...)

Calculate the mean Reef Tourism Index (RTI) for each scenario for the entire domain.
"""
function _scenario_rti(rti::YAXArray; kwargs...)
    rti_sliced = slice_results(rti; kwargs...)
    return scenario_trajectory(rti_sliced)
end
function _scenario_rti(rs::ResultSet; kwargs...)
    return _scenario_rti(reef_tourism_index(rs); kwargs...)
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
