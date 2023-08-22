
"""
    site_selection(domain::Domain, scenario::DataFrameRow{DataFrame,DataFrames.Index}, w_scens::AbstractArray, dhw_scens::AbstractArray, site_ids::Vector{Int64}, sum_cover::AbstractArray, area_to_seed::Float64)

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
function site_selection(domain::Domain, scenario::DataFrameRow, w_scens::AbstractArray, dhw_scens::AbstractArray,
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

    return ranks[:, 2:3]
end

"""
    run_site_selection(domain::Domain, scenarios::DataFrame, sum_cover::AbstractArray, area_to_seed::Float64; aggregate_function=nothing,
        target_seed_sites=nothing, target_shade_sites=nothing)
    run_site_selection(domain::Domain, scenarios::DataFrame, sum_cover::AbstractArray, area_to_seed::Float64, aggregation_function::Function;
        target_seed_sites=nothing, target_shade_sites=nothing)

Perform site selection for a given domain for multiple scenarios defined in a dataframe.

# Arguments
- `domain` : ADRIA Domain type, indicating geographical domain to perform site selection over.
- `scenarios` : DataFrame of criteria weightings and thresholds for each scenario.
- `sum_cover` : array of size (number of scenarios * number of sites) containing the summed coral cover for each site selection scenario.
- `area_to_seed` : area of coral to be seeded at each time step in km^2
- `aggregation_function` : function which aggregates ranks output into preferred form, e.g `ranks_to_frequencies`, `ranks_to_location_order`.
- `target_seed_sites` : list of candidate locations for seeding (indices)
- `target_shade_sites` : list of candidate location to shade (indices)

# Returns
- `ranks_store` : number of scenarios * sites * 3 (last dimension indicates: site_id, seed rank, shade rank)
    containing ranks for each scenario run.
"""
function run_site_selection(domain::Domain, scenarios::DataFrame, sum_cover::AbstractArray, area_to_seed::Float64;
    target_seed_sites=nothing, target_shade_sites=nothing)
    ranks_store = NamedDimsArray(
        zeros(length(domain.site_ids), 2, nrow(scenarios)),
        sites=domain.site_ids,
        intervention=["seed", "shade"],
        scenarios=1:nrow(scenarios),
    )

    dhw_scens = domain.dhw_scens
    wave_scens = domain.wave_scens

    # Pre-calculate maximum depth to consider
    target_dhw_scens = unique(scenarios[:, "dhw_scenario"])
    target_wave_scens = unique(scenarios[:, "wave_scenario"])

    target_site_ids = Int64[]
    if !isnothing(target_seed_sites)
        append!(target_site_ids, target_seed_sites)
    end

    if !isnothing(target_shade_sites)
        append!(target_site_ids, target_shade_sites)
    end

    if isnothing(target_seed_sites) && isnothing(target_shade_sites)
        target_site_ids = collect(1:length(domain.site_ids))
    end

    n_sites = length(domain.site_ids)
    for (scen_idx, scen) in enumerate(eachrow(scenarios))

        depth_criteria = (domain.site_data.depth_med .<= (scen.depth_min .+ scen.depth_offset)) .& (domain.site_data.depth_med .>= scen.depth_min)
        depth_priority = findall(depth_criteria)

        considered_sites = target_site_ids[findall(in(depth_priority), target_site_ids)]

        ranks_store(scenarios=scen_idx, sites=domain.site_ids[considered_sites]) .= site_selection(
            domain,
            scen,
            vec((mean(wave_scens[:, :, target_wave_scens], dims=(:timesteps, :scenarios)) .+ std(wave_scens[:, :, target_wave_scens], dims=(:timesteps, :scenarios))) .* 0.5),
            vec((mean(dhw_scens[:, :, target_dhw_scens], dims=(:timesteps, :scenarios)) .+ std(dhw_scens[:, :, target_dhw_scens], dims=(:timesteps, :scenarios))) .* 0.5),
            considered_sites,
            sum_cover[scen_idx, :],
            area_to_seed
        )
    end

    return ranks_store

end
function run_site_selection(domain::Domain, scenarios::DataFrame, sum_cover::AbstractArray, area_to_seed::Float64, aggregation_function::Function, iv_type::String;
    target_seed_sites=nothing, target_shade_sites=nothing)

    ranks = run_site_selection(domain, scenarios, sum_cover, area_to_seed; target_seed_sites=nothing,
        target_shade_sites=nothing)
    return aggregation_function(ranks, iv_type)

end

"""
    ranks_to_location_order(ranks::NamedDimsArray)
    ranks_to_location_order(ranks::NamedDimsArray, iv_type::String)

Post-processing function for location ranks output of `run_location_selection()`. Gives the order 
of location preference for each scenario as location ids.

# Arguments
- `ranks` : Contains location ranks for each scenario of location selection, as created by 
    `run_location_selection()`.
- `iv_type` : String indicating the intervention type to perform aggregation on.
"""
function ranks_to_location_order(ranks::NamedDimsArray)

    n_scens = length(ranks.scenarios)
    location_orders = NamedDimsArray(repeat([""], length(ranks.scenarios), length(ranks.sites)), scenarios=ranks.scenarios, ranks=1:length(ranks.sites))

    for scen in 1:n_scens
        location_orders[scenarios=scen, ranks=1:sum(ranks[scenarios=scen] .!= 0.0)] .= sort(Int.(ranks[scenarios=scen][ranks[scenarios=scen].!=0.0])).sites
    end
    return location_orders
end
function ranks_to_location_order(ranks::NamedDimsArray, iv_type::String)
    ranks_set = ranks[intervention=iv_type]
    return ranks_to_location_order(ranks_set)
end

"""
    ranks_to_frequencies(ranks::NamedDimsArray, iv_type::String, rank_frequencies::NamedDimsArray, agg_dims::Union{Symbol,Tuple{Symbol,Symbol}})
    ranks_to_frequencies(ranks::NamedDimsArray, iv_type::String)

Post-processing function for location ranks output of `run_location_selection()`. Gives the frequency 
with which each location was selected at each rank across the location selection scenarios.

# Arguments
- `ranks` : Contains location ranks for each scenario of location selection, as created by 
    `run_location_selection()`.
- `iv_type` : String indicating the intervention type to perform aggregation on.
- `rank_frequencies` : Storage container for rank frequencies.
-`agg_dims` : dimensions to aggregate over when calculating rank frequencies.

"""
function ranks_to_frequencies(ranks::NamedDimsArray, iv_type::String, rank_frequencies::NamedDimsArray, agg_dims::Union{Symbol,Tuple{Symbol,Symbol}}; n_ranks=length(ranks.sites))
    selected_ranks = ranks[intervention=iv_type]
    for rank in range(1, n_ranks, n_ranks)
        rank_frequencies[ranks=Int64(rank)] .= sum(selected_ranks .== rank, dims=agg_dims)[scenarios=1]
    end

    return rank_frequencies
end
function ranks_to_frequencies(ranks::NamedDimsArray, iv_type::String; agg_dims=:scenarios)
    if ndims(ranks) == 4
        agg_dims = (:timesteps, :scenarios)
    end
    rank_frequencies = NamedDimsArray(zeros(length(ranks.sites), length(ranks.sites)), sites=ranks.sites, ranks=1:length(ranks.sites))
    return ranks_to_frequencies(ranks, iv_type, rank_frequencies, agg_dims)
end

"""
    ranks_to_frequencies_ts(ranks::NamedDimsArray, iv_type::String; n_ranks=length(ranks.sites))

Post-processing function for location ranks output of `run_location_selection()`. Gives the frequency 
with which each location was selected at each rank across the location selection scenarios and over time.

# Arguments
- `ranks` : Contains location ranks for each scenario of location selection, as created by 
    `run_location_selection()`.
- `iv_type` : String indicating the intervention type to perform aggregation on.
- `n_ranks` : Consider first n_ranks, default is all ranks (n_locs).

"""
function ranks_to_frequencies_ts(ranks::NamedDimsArray, iv_type::String; n_ranks=length(ranks.sites))
    rank_frequencies = NamedDimsArray(zeros(length(ranks.timesteps), length(ranks.sites), length(ranks.sites)), timesteps=ranks.timesteps, sites=ranks.sites, ranks=1:length(ranks.sites))
    return ranks_to_frequencies(ranks, iv_type, rank_frequencies, :scenarios; n_ranks=n_ranks)
end

"""
    location_selection_frequencies(ranks::NamedDimsArray, ind_metrics::Vector{Int64}, iv_type::String, n_loc_int::Int64)

Post-processing function for intervention logs. Calculates the frequencies with which locations were selected for a particular intervention,
for a selection of scenarios (e.g. selected robust scenarios)
# Arguments
- `ranks` : Contains location ranks for each scenario of location selection, as created by 
    `run_location_selection()`.
- `ind_metrics` : Indices for selected scenarios (such as robust scenarios).
- `iv_type` : indicates intervention log to use ("seed", "shade" or "fog").
- `n_loc_int` : number of locations which are intervened at for each intervention decision.
- `ind_metrics` : Optional vector indicating selection of scenarios to calculate for.

"""
function location_selection_frequencies(ranks::NamedDimsArray, iv_type::String, n_loc_int::Int64; ind_metrics=collect(1:length(ranks.scenarios)))

    ranks_frequencies = ranks_to_frequencies_ts(ranks[scenarios=ind_metrics], iv_type; n_ranks=n_loc_int)
    loc_count = dropdims(sum(ranks_frequencies[ranks=1:n_loc_int], dims=[1, 3]), dims=3)[timesteps=1]

    return loc_count
end
