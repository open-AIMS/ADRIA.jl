using ADRIA: mcda_normalize

"""
    zones_criteria(zones::Vector{String}, priority_zones::Vector{String},
        strong_pred::Vector{Int64}, site_ids::Vector{Int64})

Calculates values of each site/reef for priority zones. Higher value if within or strongest 
connector for a priority zone.

# Arguments
- `zones` : Zones types for each reef/site (of length nsites).
- `priority_zones` : Priority zones for a decision instance in order of priority (entries unique).
- `strong_pred` : strongest connectivity predecessors for each reef/site (of length nsites).
- `site_ids` : full list of site ids as integers.

# Returns
- `zone_criteria` : Vector of floats indicating value of each reef/site for priority zones.

"""
function zones_criteria(zones::Vector{String}, priority_zones::Vector{String},
    strong_pred::Vector{Int64}, site_ids::Vector{Int64})
    n_sites = size(zones, 1)
    zone_ids = intersect(priority_zones, unique(zones))
    zone_weights = mcda_normalize(collect(length(zone_ids):-1:1))
    zone_preds = zeros(n_sites, 1)
    zone_sites = zeros(n_sites, 1)

    for k in axes(zone_ids, 1)
        # find sites which are strongest predecessors of sites in the zone
        zone_preds_temp = strong_pred[zones.==zone_ids[k]]
        for s in unique(zone_preds_temp)
            # for each predecessor site, add zone_weights* (no. of zone sites the site is a strongest predecessor for)
            zone_preds[site_ids.==s] .= zone_preds[site_ids.==s] .+ (zone_weights[k]) .* sum(zone_preds_temp .== s)
        end
        #     # add zone_weights for sites in the zone (whether a strongest predecessor of a zone or not)
        zone_sites[zones.==zone_ids[k]] .= zone_weights[k]
    end

    # # add weights for strongest predecessors and zones to get zone criteria
    zones_criteria = zone_preds .+ zone_sites
    return zones_criteria
end

"""
    priority_predecessor_criteria(strong_pred::Vector{Int64}, priority_sites::Vector{Any})

Calculates values of each site/reef for priority sites. Higher value if within or strongest 
connector for a priority site.

# Arguments
- `strong_pred` : strongest connectivity predecessors for each reef/site (of length nsites).
- `priority_sites` : Sites/reefs to prioritise for a management goal.

# Returns
- `predec` : Vector of floats indicating value of each reef/site for priority sites.

"""
function priority_predecessor_criteria(strong_pred::Vector{Int64}, priority_sites::Vector{Any})
    # work out which priority predecessors are connected to priority sites
    predec::Array{Float64} = zeros(length(strong_pred), 3)
    predec[:, 1:2] .= strong_pred
    predprior = predec[in.(predec[:, 1], [priority_sites']), 2]
    predprior = [x for x in predprior if !isnan(x)]

    predec[predprior, 3] .= 1.0
    return predec
end

"""
    coral_cover_criteria(site_data::DataFrame, coral_cover::Matrix{Float64})

Calculates current coral cover area and space available for coral to grow for 
each reef/site.

# Arguments
- `site_data` : DataFrame containing site/reef k values and area.
- `coral_cover` : Proportional cover at each site/reef.

"""
function coral_cover_criteria(site_data::DataFrame, coral_cover::AbstractArray)
    max_area = site_data.k .* site_data.area
    coral_cover_area = site_data.area .* coral_cover
    return max.(coral_cover_area, 0.0), max.(max_area .- coral_cover_area, 0.0)
end

"""
    env_stress_criteria(env_stress::AbstractArray)

Calculates environmental stress as a proportion of the maximum and minimum 
stress across reefs/sites.

# Arguments
- `env_stress` : e.g. heat_stress as dhws, wave_stress as probabilities etc.

"""
function env_stress_criteria(env_stress::AbstractArray)
    return maximum(env_stress) != 0.0 ? (env_stress .- minimum(env_stress)) ./ (maximum(env_stress) - minimum(env_stress)) : zeros(Float64, size(env_stress))
end

"""
    connectivity_criteria(conn::Vector{Float64}, sum_cover::Matrix,
        area::Matrix{Float64})

Calculates connectivity criterium for each reef/site as connectivity*(area of coral at site).

# Arguments
- `conn` : In-coming or out-going connectivity for each reef/site.
- `sum_cover` : Proportional coral cover for each reef/site.
- `area` : Area of each site (m^2).

"""
function connectivity_criteria(conn::Vector{Float64}, sum_cover::AbstractArray,
    area::Array{Float64})
    cov_area = conn .* sum_cover .* area
    return maximum(cov_area) != 0.0 ? cov_area / maximum(cov_area) : zeros(Float64, size(cov_area))
end

"""
    initialize_mcda(domain::Domain, param_set::NamedDimsArray, sim_params::SimConstants,
        site_data::DataFrame, depth_priority::Vector{Int64}, init_sum_cover::Array{Float64},
        area_to_seed::Float64)

Initialises variable strucutres required for dynamic site selection in ADRIA.

# Arguments
- `domain` : Domain object for site selection problem.
- `param_set` : Set of parameters for a single scenario (contains criteria weightings and tolerances).
- `site_data` : Containing site area, k values etc.
- `depth_priority` : Depth filtered set of site ids as integers.
- `init_sum_cover` : Initial proportional coral cover.
- `area_to_seed` : Area in m^2 to be covered by seeding coral at a single time step.

"""
function initialize_mcda(domain::Domain, param_set::NamedDimsArray, site_ids::Vector{Int64},
    tolerances::NamedTuple)

    rankings = [site_ids zeros(Int, length(site_ids)) zeros(Int, length(site_ids))]

    # initialize thresholds
    thresholds = create_tolerances_store(tolerances)

    # calculate values for criteria which do not change over time
    min_distance = domain.median_site_distance .* param_set("dist_thresh")

    # initialize criteria
    criteria_store = create_criteria_store(site_ids, domain.mcda_criteria)

    return rankings, criteria_store, thresholds, min_distance

end

"""
    update_criteria_store!(criteria_store::NamedDimsArray, wave_stress::AbstractArray,
        heat_stress::AbstractArray, in_conn::AbstractArray, out_conn::AbstractArray,
        site_area::AbstractArray, site_coral_cover::AbstractArray, site_data::DataFrame,
        depth_priority::AbstractArray)


Updates NamedDimsArray of criteria values required for dynamic site selection in ADRIA.

# Arguments
- `criteria_store` : Containing values of selection criteria, as created by 
    `create_criteria_store()`.
- `wave_stress` : Current wave stress values (not normalized).
- `heat_stress` : Current heat stress values (not normalized).
- `in_conn` : In-coming connectivity for each site/reef.
- `out_conn` : Out-going connectivity for each site/reef.
- `site_area` : Area in m^2 of each site/reef.
- `site_coral_cover` : Current proportional coral cover at each site/reef.
- `site_data` : Containing site area, k values etc.
- `depth_priority` : Depth filtered set of site ids as integers.

"""
function update_criteria_store!(criteria_store::NamedDimsArray, wave_stress::AbstractArray,
    heat_stress::AbstractArray, in_conn::AbstractArray, out_conn::AbstractArray,
    site_area::AbstractArray, site_coral_cover::AbstractArray, site_data::DataFrame,
    depth_priority::AbstractArray)

    criteria_store(:iv__wave_stress) .= env_stress_criteria(wave_stress)[depth_priority]
    criteria_store(:iv__heat_stress) .= env_stress_criteria(heat_stress)[depth_priority]
    criteria_store(:iv__in_connectivity) .= connectivity_criteria(in_conn, site_coral_cover, site_area)[depth_priority]
    criteria_store(:iv__in_connectivity) .= connectivity_criteria(out_conn, site_coral_cover, site_area)[depth_priority]
    coral_cover, coral_space = coral_cover_criteria(site_data, site_coral_cover')
    criteria_store(:iv__coral_cover) .= coral_cover[depth_priority]
    criteria_store(:iv__coral_space) .= coral_space[depth_priority]

    return criteria_store
end

"""
    ranks_to_location_order(ranks::NamedDimsArray, int_type::String)


Post-processing function for location ranks output of `run_site_selection()`. Gives the order 
of location preference for each scenario as location ids.

# Arguments
- `ranks` : Contains location ranks for each scenario of location selection, as created by 
    `run_site_selection()`.
- `int_type` : String indicating the intervention type to perform aggregation on.

"""
function ranks_to_location_order(ranks::NamedDimsArray, int_type::String)
    ranks_set = ranks(:, :, string(int_type, "_rank"))
    location_orders = NamedDimsArray(repeat([""], size(ranks, 1), size(ranks, 2)), scenarios=1:size(ranks, 1), ranks=1:size(ranks, 2))

    for scen in collect(1:size(ranks, 1))
        location_orders[scenarios=scen, ranks=1:sum(ranks_set[scenarios=scen] .!= 0.0)] .= sort(Int.(ranks_set[scenarios=scen][ranks_set[scenarios=scen].!=0.0])).locations
    end
    return location_orders
end

"""
    ranks_to_frequencies(ranks::NamedDimsArray, int_type::String)


Post-processing function for location ranks output of `run_site_selection()`. Gives the frequency 
with which each location was selected at each rank across the location selection scenarios.

# Arguments
- `ranks` : Contains location ranks for each scenario of location selection, as created by 
    `run_site_selection()`.
- `int_type` : String indicating the intervention type to perform aggregation on.

"""
function ranks_to_frequencies(ranks::NamedDimsArray, int_type::String)
    rank_frequencies = NamedDimsArray(zeros(size(ranks, 2), size(ranks, 2)), locations=ranks.locations, ranks=1:size(ranks, 2))

    for rank in collect(range(1, size(ranks, 2), size(ranks, 2)))
        rank_frequencies[ranks=Int(rank)] .= sum(ranks(:, :, string(int_type, "_rank")) .== rank, dims=:scenarios)[scenarios=1]
    end
    return rank_frequencies
end