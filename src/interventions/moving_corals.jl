# TODO turn into constant
# This represents how many adult settlers survive after one year for each deployed box
MC_SETTLERS_PER_BOX = 1

# TODO Update with realistic values
# basal_area_per_settler is the area in m^2 of a size class one coral
# _basal_area_per_settler = colony_mean_area(
#     corals.mean_colony_diameter_m[corals.class_id .== 1]
# )

# TODO use MCDA
# TODO filter locations with no available space
function mc_location_idx(domain::Domain, param_set::YAXArray)::Vector{Int64}
    min_mc_locations::Int64 = param_set[At("min_mc_locations")]

    # ! Just so I can implement the rest. Should be changed before merging !!
    return 1:min_mc_locations
end

function distribute_moving_corals_settlers(
    loc_k_m²::Union{Vector{Float64},SubArray{Float64,1}},
    available_space_per_loc_m2::Union{Vector{Float64},SubArray{Float64,1}},
    n_mc_settlers::Float64,
    colony_areas::Union{Vector{Float64},SubArray{Float64,1}},
    prop_fecundity::Union{Matrix{Float64},SubArray{Float64,2}}
)#::Tuple{YAXArray,Matrix{Float64}}
    # Proportion of available space on each site relative to available space at these
    # locations
    total_available_space::Float64 = sum(available_space_per_loc_m2)
    prop_available_space = available_space_per_loc_m2 ./ total_available_space

    # Intervention uses boxes that hold a certain number of settlers
    n_boxes = (n_mc_settlers / MC_SETTLERS_PER_BOX)

    # - Use prop_area_avail to determine proportion of seeds for each location
    n_boxes_per_loc = prop_available_space .* n_boxes

    box_surplus = 0.0
    for i in eachindex(n_boxes_per_loc)
        _nbpl = floor(n_boxes_per_loc[i])
        box_surplus += (n_boxes_per_loc[i] - _nbpl)
        n_boxes_per_loc[i] = _nbpl
    end

    # Designate box_surplus to location with more boxes (because it has more available space)
    n_boxes_per_loc[findmax(n_boxes_per_loc)[2]] += round(box_surplus)

    n_settlers_per_loc = n_boxes_per_loc .* MC_SETTLERS_PER_BOX

    # - Use prop_fecundity do determine how many settlers of each group per location
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

        prop_increase[:, i] .= (area_increase_per_coral_m2 ./ loc_k_m²)
    end

    return prop_increase, n_deployed_corals
end
