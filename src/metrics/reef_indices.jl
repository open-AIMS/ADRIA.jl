using FLoops, DataStructures

"""
    reef_condition_index(rc::AbstractArray, evenness::AbstractArray, sv::AbstractArray, juves::AbstractArray; threshold=2)::AbstractArray
    reef_condition_index(rs::ResultSet)::AbstractArray{<:Real}

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
YAXArray[timesteps ⋅ locations ⋅ scenarios]
"""
function _reef_condition_index(
    rc::AbstractArray,
    evenness::AbstractArray,
    sv::AbstractArray,
    juves::AbstractArray;
    threshold=2
)::AbstractArray
    # Compare outputs against reef condition criteria provided by experts
    # These are median values for 8 experts.
    condition_category = [:lower, :very_poor, :poor, :fair, :good, :very_good, :very_good]
    criteria = DataCube(
        [
            0.0 0.0 0.0 0.0
            0.05 0.15 0.175 0.15  # Very Poor [0 to these values]
            0.15 0.25 0.3 0.25    # Poor
            0.25 0.25 0.3 0.25    # Fair
            0.35 0.35 0.35 0.25   # Good
            0.45 0.45 0.45 0.35   # Very Good [these values and above]
            Inf Inf Inf Inf       # Note: some metrics might return a value > 1.0
        ];
        condition=condition_category,
        metric=[:RC, :E, :SV, :Juv]
    )
    index_metrics = zeros(size(rc)..., 4)
    @floop for (idx, met) in enumerate([rc, evenness, sv, juves])
        lower = collect(criteria[1:(end - 1), idx])
        upper = collect(criteria[2:end, idx])
        met_cp = map(x -> criteria[2:end, idx][lower .<= x .< upper][1], met.data)
        replace!(met_cp, Inf => 0.9)
        index_metrics[:, :, :, idx] .= met_cp
    end

    replace!(criteria, Inf => 0.9)

    rci = similar(rc)
    @floop for (ts, loc, scen) in Iterators.product(
        axes(index_metrics, 1),
        axes(index_metrics, 2),
        axes(index_metrics, 3))

        # Find the category name condition score relates to for each metric
        rci_cats = [
            condition_category[argmax(index_metrics[ts, loc, scen, i] .== criteria[:, i])]
            for i in axes(index_metrics, 4)
        ]
        c = counter(rci_cats)

        # Get the corresponding score values by index (based on RC only)
        # Because the index is used and not the score value (0.15, 0.25, etc), it will
        # always align with the correct category ("poor", "good", etc)
        scores = dropdims(
            criteria[condition=At(collect(keys(c))), metric=At([:RC])]; dims=2
        )

        # RCI is assigned the minimunm score of the greatest number of metrics that meet
        # `threshold`.
        # e.g., if RC and evenness are good, and sv and juves are "poor" (both meet threshold)
        #       the score is "poor".
        #       If scores are spread out (e.g., unique values), we return the 2nd lowest score
        if any(values(c) .>= threshold)
            rci[ts, loc, scen] = minimum(scores[values(c) .>= threshold])
        else
            rci[ts, loc, scen] = sort(scores)[2]
        end
    end

    clamp!(rci, 0.1, 0.9)

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
reef_condition_index = Metric(_reef_condition_index, (:timesteps, :locations, :scenarios))

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
    tac = total_absolute_cover(rs)
    rci = reef_condition_index(rs)

    return _scenario_rci(rci, tac; kwargs...)
end
scenario_rci = Metric(_scenario_rci, (:timesteps, :scenarios))

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
    rc::AbstractArray,
    evenness::AbstractArray,
    sv::AbstractArray,
    juves::AbstractArray,
    intcp_u::Vector
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
reef_tourism_index = Metric(_reef_tourism_index, (:timesteps, :locations, :scenarios))

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
scenario_rti = Metric(_scenario_rti, (:timesteps, :scenarios))

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
function _reef_fish_index(rc::AbstractArray, intcp_u1, intcp_u2)
    intcp1 = 1.232 .+ intcp_u1
    intcp2 = -1623.6 .+ intcp_u2

    slope1 = 0.007476
    slope2 = 1883.3

    # Apply unique intercepts for each scenario
    rfi = cat(
        map(
            axe ->
                0.01 .* (
                    intcp2[axe] .+
                    slope2 .*
                    (intcp1[axe] .+ slope1 .* (rc[:, :, axe] .* 100.0))
                ),
            axes(rc, 3))...;
        dims=3)

    timesteps = rfi.timesteps.val.data
    locations = rfi.locations.val.data
    n_scenarios = size(rfi, 3)
    rfi = DataCube(
        rfi.data; timesteps=timesteps, locations=locations, scenarios=1:n_scenarios
    )

    # Calculate total fish biomass, kg km2
    # 0.01 coefficient is to convert from kg ha to kg km2
    return round.(rfi, digits=2)
end
function _reef_fish_index(rs::ResultSet; intcp_u1::Bool=false, intcp_u2::Bool=false)
    rc = relative_cover(rs)
    n_scens = size(rc, :scenarios)
    icp1 = intcp_u1 ? rand(Normal(0.0, 0.195), n_scens) : zeros(n_scens)
    icp2 = intcp_u2 ? rand(Normal(0.0, 533), n_scens) : zeros(n_scens)

    return _reef_fish_index(rc, icp1, icp2)
end
reef_fish_index = Metric(_reef_fish_index, (:timesteps, :locations, :scenarios))

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
scenario_rfi = Metric(_scenario_rfi, (:timesteps, :scenarios))
