using DataStructures

"""
    reef_condition_index(rc::AbstractArray, evenness::AbstractArray, sv::AbstractArray, juves::AbstractArray)::AbstractArray
    reef_condition_index(rs)

Estimates a Reef Condition Index (RCI) providing a single value that indicates the condition
of a reef across four metrics: 
- coral cover
- evenness (coral diversity)
- shelter volume, and
- abundance of juveniles

# Notes
Juveniles are made relative to maximum observed juvenile density (51.8/m²)
See notes for `juvenile_indicator()`

# Arguments
- `rc` : Relative coral cover across all groups
- `evenness` : Evenness across all coral groups
- `sv` : Shelter volume based on coral sizes and abundances
- `juves` : Abundance of coral juveniles < 5 cm diameter

# Returns
NamedArray[timesteps ⋅ locations ⋅ scenarios]
"""
function _reef_condition_index(rc::AbstractArray, evenness::AbstractArray, sv::AbstractArray, juves::AbstractArray; threshold=2)::AbstractArray
    # Compare outputs against reef condition criteria provided by experts
    # These are median values for 8 experts.
    criteria = NamedDimsArray([
            0.0 0.0 0.0 0.0
            0.05 0.15 0.175 0.15  # Very Poor
            0.15 0.25 0.3 0.25    # Poor
            0.25 0.25 0.3 0.25    # Fair
            0.35 0.35 0.35 0.25   # Good
            0.45 0.45 0.45 0.35   # Very Good
            Inf Inf Inf Inf       # Note: some metrics might return a value > 1.0
        ],
        condition=[:lower, :very_poor, :poor, :fair, :good, :very_good, :upper],
        metric=[:RC, :E, :SV, :Juv]
    )

    index_metrics = zeros(size(rc)..., 4)
    for (idx, met) in enumerate([rc, evenness, sv, juves])
        lower = collect(criteria[1:end-1, idx])
        upper = collect(criteria[2:end, idx])
        met_cp = map(x -> criteria[2:end, idx][lower.<=x.<=upper][1], met)
        replace!(met_cp, Inf => 0.9)
        index_metrics[:, :, :, idx] .= met_cp
    end

    rci = similar(rc)
    for (ts, loc, scen) in Iterators.product(axes(index_metrics, 1), axes(index_metrics, 2), axes(index_metrics, 3))
        c = counter(index_metrics[ts, loc, scen, :])

        # RCI is assigned the minimunm score of the greatest number of metrics that meet 
        # `threshold`.
        # e.g., if RC and evenness are good, and sv and juves are "poor" (both meet threshold)
        #       the score is "poor".
        #       If scores are spread out, we return the 3rd highest score.
        if any(values(c) .>= threshold)
            rci[ts, loc, scen] = minimum(sort(collect(keys(c))[values(c).>=threshold]))
        else
            rci[ts, loc, scen] = minimum(sort(collect(keys(c))[1:3]))
        end
    end

    # TODO: bootstrap to cover expert uncertainty
    return rci
end
function _reef_condition_index(rs::ResultSet)::AbstractArray{<:Real}
    rc::AbstractArray{<:Real} = relative_cover(rs)
    juves::AbstractArray{<:Real} = juvenile_indicator(rs)
    evenness::AbstractArray{<:Real} = coral_evenness(rs)
    sv::AbstractArray{<:Real} = relative_shelter_volume(rs)

    return _reef_condition_index(rc, evenness, sv, juves)
end
reef_condition_index = Metric(_reef_condition_index, (:timesteps, :sites, :scenarios))

"""
    scenario_rci(rci::NamedDimsArray; kwargs...)
    scenario_rci(rs::ResultSet; kwargs...)

Calculate the mean Reef Fish Index for each scenario for the entire domain.
"""
function _scenario_rci(rti::NamedDimsArray; kwargs...)
    rti_sliced = slice_results(rti; kwargs...)
    return scenario_trajectory(rti_sliced)
end
function _scenario_rci(rs::ResultSet; kwargs...)
    return _scenario_rci(reef_condition_index(rs); kwargs...)
end
scenario_rci = Metric(_scenario_rci, (:timesteps, :scenario))

"""
    reef_tourism_index(rc::AbstractArray, evenness::AbstractArray, sv::AbstractArray, juves::AbstractArray)::AbstractArray
    reef_tourism_index(rs::ResultSet)::AbstractArray

Estimate tourism index.

# Arguments
- `rc` : Relative coral cover across all groups
- `evenness` : Evenness across all coral groups
- `sv` : Shelter volume based on coral sizes and abundances
- `juves` : Abundance of coral juveniles < 5 cm diameter
"""
function _reef_tourism_index(rc::AbstractArray, evenness::AbstractArray, sv::AbstractArray, juves::AbstractArray, intcp_u::Vector)::AbstractArray
    intcp = -0.498 .+ intcp_u

    # TODO: Ryan to reinterpolate to account for no CoTS and no rubble
    # Apply unique intercepts for each scenario
    rti = cat(map(axe -> (intcp[axe] .+ (0.291 .* rc[:, :, axe]) .+
                          (0.056 .* evenness[:, :, axe]) .+
                          (0.628 .* sv[:, :, axe]) .+
                          (1.335 .* juves[:, :, axe])
            ), axes(rc, 3))..., dims=3)

    return round.(clamp.(rti, 0.1, 0.9), digits=2)
end
function _reef_tourism_index(rs::ResultSet; intcp_u::Bool=true)::AbstractArray
    rc::AbstractArray{<:Real} = relative_cover(rs)
    juves::AbstractArray{<:Real} = juvenile_indicator(rs)
    evenness::AbstractArray{<:Real} = coral_evenness(rs)
    sv::AbstractArray{<:Real} = relative_shelter_volume(rs)

    n_scens = size(rc, :scenarios)
    intcp = intcp_u ? rand(Normal(0.0, 0.163), n_scens) : zeros(n_scens)

    return _reef_tourism_index(rc, evenness, sv, juves, intcp)
end
reef_tourism_index = Metric(_reef_tourism_index, (:timesteps, :sites, :scenarios))

"""
    scenario_rti(rti::NamedDimsArray; kwargs...)
    scenario_rti(rs::ResultSet; kwargs...)

Calculate the mean Reef Fish Index for each scenario for the entire domain.
"""
function _scenario_rti(rti::NamedDimsArray; kwargs...)
    rti_sliced = slice_results(rti; kwargs...)
    return scenario_trajectory(rti_sliced)
end
function _scenario_rti(rs::ResultSet; kwargs...)
    return _scenario_rti(reef_tourism_index(rs); kwargs...)
end
scenario_rti = Metric(_scenario_rti, (:timesteps, :scenario))


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
NamedArray[timesteps ⋅ locations ⋅ scenarios], values in kg/km²

# References
1. Graham, N.A.J., Nash, K.L., 2013. 
The importance of structural complexity in coral reef ecosystems. 
Coral Reefs 32, 315–326. 
https://doi.org/10.1007/s00338-012-0984-y
"""
function _reef_fish_index(rc::AbstractArray, intcp_u1, intcp_u2)
    intcp1 = 1.232 .+ intcp_u1
    intcp2 = -1623.6 .+ intcp_u2

    slope1 = 0.007476
    slope2 = 1883.3

    # Apply unique intercepts for each scenario
    rfi = cat(map(axe -> 0.01 .* (
                intcp2[axe] .+ slope2 .*
                               (intcp1[axe] .+ slope1 .* (rc[:, :, axe] .* 100.0))
            ),
            axes(rc, 3))...,
        dims=3)

    # Calculate total fish biomass, kg km2
    # 0.01 coefficient is to convert from kg ha to kg km2
    return round.(rfi, digits=2)
end
function _reef_fish_index(rs::ResultSet; intcp_u1::Bool=true, intcp_u2::Bool=true)
    rc = relative_cover(rs)
    n_scens = size(rc, :scenarios)
    icp1 = intcp_u1 ? rand(Normal(0.0, 0.195), n_scens) : zeros(n_scens)
    icp2 = intcp_u2 ? rand(Normal(0.0, 533), n_scens) : zeros(n_scens)

    return _reef_fish_index(rc, icp1, icp2)
end
reef_fish_index = Metric(_reef_fish_index, (:timesteps, :sites, :scenarios))

"""
    scenario_rfi(rfi::NamedDimsArray; kwargs...)
    scenario_rfi(rs::ResultSet; kwargs...)

Calculate the mean Reef Fish Index for each scenario for the entire domain.
"""
function _scenario_rfi(rfi::NamedDimsArray; kwargs...)
    rfi_sliced = slice_results(rfi; kwargs...)
    return scenario_trajectory(rfi_sliced)
end
function _scenario_rfi(rs::ResultSet; kwargs...)
    return _scenario_rfi(reef_fish_index(rs); kwargs...)
end
scenario_rfi = Metric(_scenario_rfi, (:timesteps, :scenario))
