using ADRIA: mcda_normalize


function zones_criteria(zones::Vector{String}, priority_zone::Vector{String})::AbstractArray
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
    coral_cover_area = max_area .* coral_cover
    coral_cover_space = max_area .- coral_cover_area
    return coral_cover_area, coral_cover_space
end