
"""
    site_selection(domain::Domain, scenario::DataFrameRow{DataFrame,DataFrames.Index}, w_scens::NamedDimsArray, dhw_scens::NamedDimsArray, site_ids::Vector{Int64}, sum_cover::AbstractArray, area_to_seed::Float64)

Perform site selection using a chosen mcda aggregation method, domain, initial cover, criteria weightings and thresholds.

# Arguments
- `domain` : ADRIA Domain type, indicating geographical domain to perform site selection over.
- `scenario` : contains criteria weightings and thresholds for a single scenario.
- `w_scens` : array of length nsites containing wave scenario.
- `dhw_scens` : array of length nsites containing dhw scenario.
- `site_ids` : locations to consider
- `sum_cover` : summed cover (over species) for single scenario being run, for each site.
- `area_to_seed` : area of coral to be seeded at each time step in km^2

# Returns
- `ranks` : n_reps * sites * 3 (last dimension indicates: site_id, seeding rank, shading rank)
    containing ranks for single scenario.
"""
function site_selection(domain::Domain, scenario::DataFrameRow, w_scens::NamedDimsArray, dhw_scens::NamedDimsArray,
    site_ids::Vector{Int64}, sum_cover::NamedDimsArray, area_to_seed::Float64)::Matrix{Int64}

    mcda_vars = DMCDA_vars(domain, scenario, site_ids, sum_cover, area_to_seed, w_scens, dhw_scens)
    n_sites = length(site_ids)

    # site_id, seeding rank, shading rank
    rankingsin = [mcda_vars.site_ids zeros(Int64, n_sites) zeros(Int64, n_sites)]

    prefseedsites::Vector{Int64} = zeros(Int64, mcda_vars.n_site_int)
    prefshadesites::Vector{Int64} = zeros(Int64, mcda_vars.n_site_int)

    # Determine connectivity strength
    # Account for cases where no coral cover
    in_conn, out_conn, strong_pred = connectivity_strength(domain.TP_data .* site_k_area(domain), Array(sum_cover))
    in_conn = in_conn[site_ids]
    out_conn = out_conn[site_ids]
    strong_pred = strong_pred[site_ids]

    (_, _, ranks) = guided_site_selection(mcda_vars, scenario.guided, true, true, prefseedsites, prefshadesites, rankingsin, in_conn, out_conn, strong_pred)

    return ranks
end

"""
    run_site_selection(domain::Domain, scenarios::DataFrame, sum_cover::AbstractArray, area_to_seed::Float64; aggregate_function=nothing,
        target_seed_sites=nothing, target_shade_sites=nothing)

Perform site selection for a given domain for multiple scenarios defined in a dataframe.

# Arguments
- `domain` : ADRIA Domain type, indicating geographical domain to perform site selection over.
- `scenarios` : DataFrame of criteria weightings and thresholds for each scenario.
- `sum_cover` : array of size (number of scenarios * number of sites) containing the summed coral cover for each site selection scenario.
- `area_to_seed` : area of coral to be seeded at each time step in km^2
- `aggregate_function` : function which aggregates ranks output into preferred form, e.g `ranks_to_frequencies`, `ranks_to_location_order`.
- `target_seed_sites` : list of candidate locations for seeding (indices)
- `target_shade_sites` : list of candidate location to shade (indices)

# Returns
- `ranks_store` : number of scenarios * sites * 3 (last dimension indicates: site_id, seed rank, shade rank)
    containing ranks for each scenario run.
"""
function run_site_selection(domain::Domain, scenarios::DataFrame, sum_cover::AbstractArray, area_to_seed::Float64; aggregate_function=nothing,
    target_seed_sites=nothing, target_shade_sites=nothing)
    ranks_store = NamedDimsArray(
        zeros(nrow(scenarios), length(domain.site_ids), 3),
        scenarios=1:nrow(scenarios),
        sites=domain.site_ids,
        ranks=["site_id", "seed_rank", "shade_rank"],
    )

    dhw_scens = domain.dhw_scens
    wave_scens = domain.wave_scens

    # Pre-calculate maximum depth to consider
    scenarios[:, "max_depth"] .= scenarios.depth_min .+ scenarios.depth_offset
    target_dhw_scens = unique(scenarios[:, "dhw_scenario"])
    target_wave_scens = unique(scenarios[:, "wave_scenario"])

    target_site_ids = Int64[]
    if !isnothing(target_seed_sites)
        append!(target_site_ids, target_seed_sites)
    end

    if !isnothing(target_shade_sites)
        append!(target_site_ids, target_shade_sites)
    end

    n_sites = length(dom.site_ids)
    for (scen_idx, scen) in enumerate(eachrow(scenarios))
        depth_criteria = (dom.site_data.depth_med .<= scen.max_depth) .& (dom.site_data.depth_med .>= scen.depth_min)
        depth_priority = findall(depth_criteria)

        considered_sites = target_site_ids[findall(in(depth_priority), target_site_ids)]
        ranks_store(scenarios=scen_idx, sites=dom.site_ids[considered_sites]) .= site_selection(
            dom,
            scen,
            (mean(wave_scens[:, :, target_wave_scens], dims=(:timesteps, :scenarios)) .+ std(wave_scens[:, :, target_wave_scens], dims=(:timesteps, :scenarios))) .* 0.5,
            (mean(dhw_scens[:, :, target_dhw_scens], dims=(:timesteps, :scenarios)) .+ std(dhw_scens[:, :, target_dhw_scens], dims=(:timesteps, :scenarios))) .* 0.5,
            considered_sites,
            sum_cover[scen_idx, :],
            area_to_seed
        )
    end
    if !isnothing(aggregate_function)
        return ranks_store, aggregate_function(ranks_store, "seed"), aggregate_function(ranks_store, "shade")
    else
        return ranks_store
    end
end

"""
    ranks_to_location_order(ranks::NamedDimsArray, interv_type::String)


Post-processing function for location ranks output of `run_location_selection()`. Gives the order 
of location preference for each scenario as location ids.

# Arguments
- `ranks` : Contains location ranks for each scenario of location selection, as created by 
    `run_location_selection()`.
- `interv_type` : String indicating the intervention type to perform aggregation on.
"""
function ranks_to_location_order(ranks::NamedDimsArray, interv_type::String)
    ranks_set = ranks(:, :, string(interv_type, "_rank"))
    location_orders = NamedDimsArray(repeat([""], size(ranks, 1), size(ranks, 2)), scenarios=1:size(ranks, 1), ranks=1:size(ranks, 2))

    for scen in 1:size(ranks, 1)
        location_orders[scenarios=scen, ranks=1:sum(ranks_set[scenarios=scen] .!= 0.0)] .= sort(Int.(ranks_set[scenarios=scen][ranks_set[scenarios=scen].!=0.0])).sites
    end
    return location_orders
end

"""
    ranks_to_frequencies(ranks::NamedDimsArray, interv_type::String)


Post-processing function for location ranks output of `run_location_selection()`. Gives the frequency 
with which each location was selected at each rank across the location selection scenarios.

# Arguments
- `ranks` : Contains location ranks for each scenario of location selection, as created by 
    `run_location_selection()`.
- `interv_type` : String indicating the intervention type to perform aggregation on.
"""
function ranks_to_frequencies(ranks::NamedDimsArray, interv_type::String)
    rank_frequencies = NamedDimsArray(zeros(size(ranks, 2), size(ranks, 2)), sites=ranks.sites, ranks=1:size(ranks, 2))

    for rank in range(1, size(ranks, 2), size(ranks, 2))
        rank_frequencies[ranks=Int64(rank)] .= sum(ranks(:, :, string(interv_type, "_rank")) .== rank, dims=:scenarios)[scenarios=1]
    end
    return rank_frequencies
end
