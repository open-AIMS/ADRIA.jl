
"""
    _site_selection(domain::Domain, mcda_vars::DMCDA_vars, guided::Int64)

Perform site selection using a chosen aggregation method, domain, initial cover, criteria weightings and thresholds.

# Arguments
- `domain` : ADRIA Domain type, indicating geographical domain to perform site selection over.
- `mcda_vars` : Contains relevant parameters for performing guided location selection.
- `guided` : Integer indicating aggegation algorithm to use for guided location selection.

# Returns
- `ranks` : n_reps * sites * 3 (last dimension indicates: site_id, seeding rank, shading rank)
    containing ranks for single scenario.
"""
function _site_selection(domain::Domain, mcda_vars::DMCDA_vars, guided::Int64)
    site_ids = mcda_vars.site_ids
    n_sites = length(site_ids)

    # site_id, seeding rank, shading rank
    rankingsin = [mcda_vars.site_ids zeros(Int64, n_sites) zeros(Int64, n_sites)]

    prefseedsites::Vector{Int64} = zeros(Int64, mcda_vars.n_site_int)
    prefshadesites::Vector{Int64} = zeros(Int64, mcda_vars.n_site_int)

    # Determine connectivity strength
    # Account for cases where no coral cover
    in_conn, out_conn, strong_pred = connectivity_strength(domain.TP_data .* site_k_area(domain), Array(mcda_vars.sum_cover))

    # Perform location selection for seeding and shading.
    seed_true, shade_true = [true, true]

    (_, _, ranks) = guided_site_selection(mcda_vars, guided, seed_true, shade_true, prefseedsites, prefshadesites, rankingsin, in_conn[site_ids], out_conn[site_ids], strong_pred[site_ids])

    return ranks[:, 2:3]
end

"""
    run_site_selection(domain::Domain, scenarios::DataFrame, sum_cover::NamedDimsArray, area_to_seed::Float64;
        target_seed_sites=nothing, target_shade_sites=nothing)
    run_site_selection(domain::Domain, scenarios::DataFrame, sum_cover::NamedDimsArray, area_to_seed::Float64, 
        aggregation_function::Function; target_seed_sites=nothing, target_shade_sites=nothing)

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
function run_site_selection(domain::Domain, scenarios::DataFrame, sum_cover::NamedDimsArray, area_to_seed::Float64;
    target_seed_sites=nothing, target_shade_sites=nothing)
    ranks_store = NamedDimsArray(
        zeros(length(domain.site_ids), 2, nrow(scenarios)),
        sites=domain.site_ids,
        intervention=1:2,
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
        mcda_vars_temp = DMCDA_vars(domain, scen, considered_sites,  sum_cover[scen_idx, :], area_to_seed, 
            vec((mean(wave_scens[:, :, target_wave_scens], dims=(:timesteps, :scenarios)) .+ std(wave_scens[:, :, target_wave_scens], dims=(:timesteps, :scenarios))) .* 0.5),
            vec((mean(dhw_scens[:, :, target_dhw_scens], dims=(:timesteps, :scenarios)) .+ std(dhw_scens[:, :, target_dhw_scens], dims=(:timesteps, :scenarios))) .* 0.5),
            )
    
        ranks_store(scenarios=scen_idx, sites=domain.site_ids[considered_sites]) .= _site_selection(domain, mcda_vars_temp, scen.guided)
    end

    return ranks_store

end
function run_site_selection(domain::Domain, scenarios::DataFrame, sum_cover::NamedDimsArray, area_to_seed::Float64, aggregation_function::Function, iv_type::String;
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

# Returns
- Location order after ranking as string IDs for each scenario.
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
    ranks_set = _get_iv_type(ranks, iv_type)
    return ranks_to_location_order(ranks_set)
end

"""
    ranks_to_frequencies(ranks::NamedDimsArray, iv_type::String; n_ranks=length(ranks.sites))
    ranks_to_frequencies(rs::ResultSet, iv_type::String; n_ranks=length(ranks.sites))

Post-processing function for location ranks output of `run_location_selection()`. Gives the frequency 
with which each location was selected at each rank across the location selection scenarios.

# Arguments
- `ranks` : Contains location ranks for each scenario of location selection, as created by 
    `run_location_selection()`.
- `rs` : ADRIA result set.
- `iv_type` : String indicating the intervention type to perform aggregation on.

# Returns
- Frequency with which each location was selected for each rank.
"""
function ranks_to_frequencies(ranks::NamedDimsArray, rank_frequencies::NamedDimsArray; n_ranks::Int64=length(ranks.sites))
    for rank in range(1, n_ranks, n_ranks)
        rank_frequencies[ranks=Int64(rank)] .= sum(ranks .== rank, dims=:scenarios)[scenarios=1]
    end
    return rank_frequencies
end
function ranks_to_frequencies(ranks::NamedDimsArray, iv_type::String; n_ranks::Int64=length(ranks.sites))
    selected_ranks = _get_iv_type(ranks, iv_type)
    rank_frequencies = NamedDimsArray(zeros(length(ranks.sites), length(ranks.sites)), sites=ranks.sites, ranks=1:length(ranks.sites))
    return ranks_to_frequencies(selected_ranks, rank_frequencies; n_ranks=n_ranks)
end
function ranks_to_frequencies(rs::ResultSet, iv_type::String; n_ranks=length(rs.ranks.sites))
    return sum(ranks_to_frequencies_ts(rs, iv_type; n_ranks=n_ranks), dims=:timesteps)[timesteps=1]
end

"""
    ranks_to_frequencies_ts(ranks::NamedDimsArray, iv_type::String; n_ranks=length(ranks.sites))
    ranks_to_frequencies_ts(rs::ResultSet, iv_type::String; n_ranks=length(ranks.sites))

Post-processing function for location ranks output of `run_location_selection()`. Gives the frequency 
with which each location was selected at each rank across the location selection scenarios and over time.

# Arguments
- `ranks` : Contains location ranks for each scenario of location selection, as created by 
    `run_location_selection()`.
- `rs` : ADRIA result set. 
- `iv_type` : String indicating the intervention type to perform aggregation on.
- `n_ranks` : Consider first n_ranks, default is all ranks (n_locs).

# Returns 
- Frequency with which each location was selected for each rank over time.
"""
function ranks_to_frequencies_ts(ranks::NamedDimsArray; n_ranks::Int64=length(ranks.sites))
    rank_frequencies = NamedDimsArray(zeros(length(ranks.timesteps), length(ranks.sites), length(ranks.sites)), timesteps=ranks.timesteps, sites=ranks.sites, ranks=1:length(ranks.sites))
    return ranks_to_frequencies(ranks, rank_frequencies; n_ranks=n_ranks)
end
function ranks_to_frequencies_ts(ranks::NamedDimsArray, iv_type::String; n_ranks=length(ranks.sites))
    selected_ranks = _get_iv_type(ranks, iv_type)
    return ranks_to_frequencies_ts(selected_ranks; n_ranks=n_ranks)
end
function ranks_to_frequencies_ts(rs::ResultSet, iv_type::String; n_ranks=length(rs.ranks.sites))
    selected_ranks = _get_iv_type(rs.ranks, iv_type)
    return ranks_to_frequencies_ts(selected_ranks; n_ranks=n_ranks)
end

"""
    location_selection_frequencies(ranks::NamedDimsArray, iv_type::String; n_loc_int=5, ind_metrics=collect(1:length(ranks.scenarios)))
    location_selection_frequencies(rs::ResultSet, iv_type::String; n_loc_int=5, ind_metrics=collect(1:length(ranks.scenarios)))
 
Post-processing function for intervention logs. Calculates the frequencies with which locations were selected for a particular intervention,
for a selection of scenarios (e.g. selected robust scenarios).

# Arguments
- `ranks` : Contains location ranks for each scenario of location selection, as created by 
    `run_location_selection()`.
- `rs` : ADRIA result set.
- `ind_metrics` : Indices for selected scenarios (such as robust scenarios).
- `iv_type` : indicates intervention log to use ("seed", "shade" or "fog").
- `n_loc_int` : number of locations which are intervened at for each intervention decision.

# Returns 
- Counts for location selection at each location in the domain.
"""
function location_selection_frequencies(ranks::NamedDimsArray, iv_type::String; n_loc_int::Int64=5, ind_metrics::Vector{Int64}=ranks.scenarios)
    ranks_frequencies = ranks_to_frequencies(ranks[scenarios=ind_metrics], iv_type; n_ranks=n_loc_int)
    loc_count = sum(ranks_frequencies[ranks=1:n_loc_int], dims=2)[ranks=1]

    return loc_count
end
function location_selection_frequencies(rs::ResultSet, iv_type::String; n_loc_int::Int64=5, ind_metrics::Vector{Int64}=rs.ranks.scenarios)
    selected_ranks = _get_iv_type(rs.ranks[scenarios=ind_metrics], iv_type)
    ranks_frequencies = ranks_to_frequencies_ts(selected_ranks; n_ranks=n_loc_int)
    loc_count = dropdims(sum(ranks_frequencies[ranks=1:n_loc_int], dims=[1, 3]), dims=3)[timesteps=1]

    return loc_count
end

"""
    summed_inverse_rank(ranks::NamedDimsArray, iv_type::String; dims=:scenarios,agg_func=x->x)
    summed_inverse_rank(rs::ResultSet, iv_type::String; dims=:scenarios,agg_func=x->x)
    summed_inverse_rank(ranks::NamedDimsArray; dims=:scenarios,agg_func=x->x)
 
Calculates (number of sites) .- ranks summed over the dimension dims and transformed using agg_func 
    (default no transformation).

# Arguments
- `ranks` : Contains location ranks for each scenario of location selection, as created by 
    `run_location_selection()`.
- `rs` : ADRIA result set.
- `iv_type` : indicates intervention log to use ("seed", "shade" or "fog").
- `dims` : Dimensions to sum over.
- `agg_func` : function to transform result.

# Returns 
- Inverse rankings (i.e. the greater the number the higher ranked the site).
"""
function summed_inverse_rank(ranks::NamedDimsArray, iv_type::String; dims::Symbol=:scenarios, agg_func::Function=x->x)
    selected_ranks = _get_iv_type(ranks, iv_type)
    return summed_inverse_rank(selected_ranks; dims=dims,agg_func=agg_func)

end
function summed_inverse_rank(rs::ResultSet, iv_type::String; dims::Symbol=:scenarios, agg_func=x->x)
    selected_ranks = _get_iv_type(rs.ranks, iv_type)
    return summed_inverse_rank(selected_ranks; dims=dims,agg_func=agg_func)

end
function summed_inverse_rank(ranks::NamedDimsArray; dims::Symbol=:scenarios,agg_func::Function=x->x)
    n_locs = size(ranks,1)
    inv_ranks = agg_func(sum(ranks .- n_locs,dims=dims))
    
    return inv_ranks
end


"""
    _get_iv_type(ranks, iv_type)

Get ranks for intervention based on name.
    
# Arguments
- `ranks` : Contains location ranks for each scenario of location selection, as created by 
    `run_location_selection()`.
- `iv_type` : indicates intervention log to use ("seed", "shade" or "fog").

# Returns
- Ranks for specified iv_type.
"""
function _get_iv_type(ranks, iv_type)
    iv_dict = Dict([("seed", 1), ("shade", 2)])
    return ranks[intervention=iv_dict[iv_type]]
end

