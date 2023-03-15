using ADRIA: mcda_normalize, DMCDA_vars


function zones_criteria(zones, priority_zones, strong_pred, site_ids)
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

function priority_predecessor_criteria(strong_pred::AbstractArray, priority_sites::AbstractArray)::AbstractArray
    # # work out which priority predecessors are connected to priority sites
    predec::Array{Float64} = zeros(n_sites, 3)
    predec[:, 1:2] .= strong_pred
    predprior = predec[in.(predec[:, 1], [priority_sites']), 2]
    predprior = [x for x in predprior if !isnan(x)]

    predec[predprior, 3] .= 1.0
    return predec
end

function coral_cover_criteria(site_data, coral_cover)
    max_area = site_data.k .* site_data.area
    return max.(max_area .* coral_cover, 0.0), max.(max_area .- coral_cover_area, 0.0)
end

function env_stress_criteria(env_stress)
    return maximum(env_stress) != 0.0 ? (env_stress .- minimum(env_stress)) ./ (maximum(env_stress) - minimum(env_stress)) : 0.0
end

function connectivity_criteria(conn)
    cov_area = conn .* sum_cover .* area
    return maximum(cov_area) != 0.0 ? cov_area / maximum(cov_area) : 0.0
end