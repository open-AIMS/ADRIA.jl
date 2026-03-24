# This represents how many corals survive after one year for each deployment
# TODO Use a realistic value for this
const MC_CORALS_PER_POOL = 1

"""
    distribute_moving_corals(
        loc_k_m²::Union{Vector{Float64},SubArray{Float64,1}},
        available_space_per_loc_m2::Union{Vector{Float64},SubArray{Float64,1}},
        n_mc_settlers::Float64,
        colony_areas::Union{Vector{Float64},SubArray{Float64,1}},
        prop_fecundity::Union{Matrix{Float64},SubArray{Float64,2}}
    )::Tuple{Matrix{Float64},Matrix{Float64}}

Distribute pools for moving corals intervention. Number of corals per pool is determined by
the constant `MC_CORALS_PER_POOL`. Number of pools follows available space relative to the
overall available space across all selected locations.

# Arguments
- `loc_k_m²` : Carrying capacity area of locations to seed in m².
- `available_space_per_loc_m2` : Currently available space at each seed location in m².
- `n_mc_settlers` : Number of yo settlers to add
- `colony_areas` : Area of one coral of each functional group and size class
- `prop_fecundity` : Proportional fecundity to infer
"""
function distribute_moving_corals(
    loc_k_m²::Union{Vector{Float64},SubArray{Float64,1}},
    available_space_per_loc_m2::Union{Vector{Float64},SubArray{Float64,1}},
    n_mc_settlers::Float64,
    colony_areas::Union{Vector{Float64},SubArray{Float64,1}},
    prop_fecundity::Union{Matrix{Float64},SubArray{Float64,2}}
)::Tuple{Matrix{Float64},Matrix{Float64}}
    # Proportion of available space on each location
    total_available_space::Float64 = sum(available_space_per_loc_m2)
    prop_available_space = available_space_per_loc_m2 ./ total_available_space

    n_pools = (n_mc_settlers / MC_CORALS_PER_POOL)

    # Use prop_area_avail to determine proportion of corals for each location
    n_pools_per_loc = prop_available_space .* n_pools

    pool_surplus = 0.0
    for i in eachindex(n_pools_per_loc)
        _nppl = floor(n_pools_per_loc[i])
        pool_surplus += (n_pools_per_loc[i] - _nppl)
        n_pools_per_loc[i] = _nppl
    end

    # Designate pool_surplus to location with more pools (because it has more available space)
    n_pools_per_loc[findmax(n_pools_per_loc)[2]] += round(pool_surplus)
    n_settlers_per_loc = n_pools_per_loc .* MC_CORALS_PER_POOL

    # Use prop_fecundity do determine how many settlers of each group per location
    n_deployed_corals = zeros(size(prop_fecundity))
    prop_increase = zeros(size(prop_fecundity))
    for i in eachindex(n_settlers_per_loc)
        # ? Is it ok to just round here?
        # ? The implication is: if you sum this it won't match exactly n_settlers_per_loc
        n_deployed_corals[:, i] .= round.(prop_fecundity[:, i] .* n_settlers_per_loc[i])

        area_increase_per_coral_m2 = n_deployed_corals[:, i] .* colony_areas

        # If increase in area will exceed available space, decrease the
        if sum(area_increase_per_coral_m2) > available_space_per_loc_m2[i]
            @warn "MC area increase exceeds available space for some locations. " .*
                "Constraining area and number of larvae."

            # Cap increase in cover for each group proportional to the increase in each one
            area_increase_per_coral_m2 -=
                available_space_per_loc_m2[i] .*
                (area_increase_per_coral_m2 / sum(area_increase_per_coral_m2))

            # Adjust number of deployed corals accordingly
            n_deployed_corals[:, i] .=
                floor.(
                    area_increase_per_coral_m2 ./ colony_areas
                )
        end

        prop_increase[:, i] .= (area_increase_per_coral_m2 ./ loc_k_m²[i])
    end

    return prop_increase, n_deployed_corals
end
