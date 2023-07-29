"""
    reef_condition_index(rc, E, SV, juveniles)
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
function _reef_condition_index(rc::AbstractArray, evenness::AbstractArray, sv::AbstractArray, juves::AbstractArray)::AbstractArray
    # Compare outputs against reef condition criteria provided by experts

    # These are median values for 7 experts. TODO: draw from distributions
    #  Condition        RC       E       SV      Juv
    # 'VeryGood'        0.45     0.45    0.45    0.35
    # 'Good'            0.35     0.35    0.35    0.25
    # 'Fair'            0.25     0.25    0.30    0.25
    # 'Poor'            0.15     0.25    0.30    0.25
    # 'VeryPoor'        0.05     0.15    0.18    0.15

    # Note that the scores for evenness and juveniles are slightly different
    lin_grid::Gridded{Linear{Throw{OnGrid}}} = Gridded(Linear())
    TC_func::GriddedInterpolation{Float64,1,Vector{Float64},Gridded{Linear{Throw{OnGrid}}},Tuple{Vector{Float64}}} = interpolate((Float64[0, 0.05, 0.15, 0.25, 0.35, 0.45, 1.0],), Float64[0, 0.1, 0.3, 0.5, 0.7, 0.9, 1.0], lin_grid)
    E_func::GriddedInterpolation{Float64,1,Vector{Float64},Gridded{Linear{Throw{OnGrid}}},Tuple{Vector{Float64}}} = interpolate((Float64[0, 0.15, 0.25, 0.35, 0.45, 1.0],), Float64[0, 0.1, 0.5, 0.7, 0.9, 1.0], lin_grid)
    SV_func::GriddedInterpolation{Float64,1,Vector{Float64},Gridded{Linear{Throw{OnGrid}}},Tuple{Vector{Float64}}} = interpolate((Float64[0, 0.18, 0.30, 0.35, 0.45, 1.0],), Float64[0, 0.1, 0.3, 0.5, 0.9, 1.0], lin_grid)
    juv_func::GriddedInterpolation{Float64,1,Vector{Float64},Gridded{Linear{Throw{OnGrid}}},Tuple{Vector{Float64}}} = interpolate((Float64[0, 0.15, 0.25, 0.35, 1.0],), Float64[0, 0.1, 0.5, 0.9, 1.0], lin_grid)

    rc_i::AbstractArray{<:Real} = TC_func.(rc)
    E_i::AbstractArray{<:Real} = E_func.(evenness)
    SV_i::AbstractArray{<:Real} = SV_func.(sv)
    juv_i::AbstractArray{<:Real} = juv_func.(juves)

    return mean([rc_i, E_i, SV_i, juv_i])
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
    reef_tourism_index(rc::AbstractArray, evenness::AbstractArray, sv::AbstractArray, juves::AbstractArray)::AbstractArray
    reef_tourism_index(rs::ResultSet)::AbstractArray

Estimate tourism index.

# Arguments
- `rc` : Relative coral cover across all groups
- `evenness` : Evenness across all coral groups
- `sv` : Shelter volume based on coral sizes and abundances
- `juves` : Abundance of coral juveniles < 5 cm diameter
"""
function _reef_tourism_index(rc::AbstractArray, evenness::AbstractArray, sv::AbstractArray, juves::AbstractArray; intcp_u=0.0)::AbstractArray
    intcp = -0.498 + intcp_u
    rti = (intcp .+ (0.291 .* rc) .+
           (0.056 .* evenness) .+
           (0.628 .* sv) .+
           (1.335 .* juves)
    )

    return clamp.(rti, 0.1, 0.9)
end
function _reef_tourism_index(rs::ResultSet; intcp_u::Bool=true)::AbstractArray
    rc::AbstractArray{<:Real} = relative_cover(rs)
    juves::AbstractArray{<:Real} = juvenile_indicator(rs)
    evenness::AbstractArray{<:Real} = coral_evenness(rs)
    sv::AbstractArray{<:Real} = relative_shelter_volume(rs)

    if intcp_u
        intcp = rand(Normal(0.0, 0.163))
    end

    return _reef_tourism_index(rc, evenness, sv, juves; intcp_u=intcp)
end
reef_tourism_index = Metric(_reef_tourism_index, (:timesteps, :sites, :scenarios))

"""
    reef_fish_index(rc::AbstractArray)
    reef_fish_index(rs::ResultSet)

The Reef Fish Index (RFI) estimates fish biomass from relative coral cover.
The relationship is developed by R. Heneghan (QUT) by digitizing Figures 4a and 6b
from Graham & Nash (2013; see [1])

A linear regression (developed by R. Heneghan, Queensland University of Technology) 
is used to indicate the relationship between coral cover and fish biomass.

Values are provided ∈ [0, 1], where 1 indicates maximum fish biomass.

Note: Coral cover here is relative to coral habitable area (\$k\$ area).

# Arguments
- `rc` : Relative cover

# Returns
NamedArray[timesteps ⋅ locations ⋅ scenarios]

# References
1. Graham, N.A.J., Nash, K.L., 2013. 
The importance of structural complexity in coral reef ecosystems. 
Coral Reefs 32, 315–326. 
https://doi.org/10.1007/s00338-012-0984-y
"""
function _reef_fish_index(rc::AbstractArray; intcp_u1=0.0, intcp_u2=0.0)
    intcp1 = 1.232 + intcp_u1
    intcp2 = -1623.6 + intcp_u2

    slope1 = 0.007476
    slope2 = 1883.3

    return 0.01 * (intcp2 .+ slope2 .* (intcp1 .+ slope1 .* (rc .* 100.0))) ./ 100.0
end
function _reef_fish_index(rs::ResultSet; intcp_u1::Bool=true, intcp_u2::Bool=true)
    icp1 = intcp_u1 ? rand(Normal(0.0, 0.195)) : 0.0
    icp2 = intcp_u2 ? rand(Normal(0.0, 533)) : 0.0

    return _reef_fish_index(_relative_cover(rs); intcp_u1=icp1, intcp_u2=icp2)
end
reef_fish_index = Metric(_reef_fish_index, (:timesteps, :sites, :scenarios))
