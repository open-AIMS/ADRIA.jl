using NamedDims, AxisKeys


"""
    _site_selection(domain::Domain, mcda_vars::DMCDA_vars, guided::Int64)

Perform site selection using a chosen aggregation method, domain, initial cover, criteria weightings and thresholds.

# Arguments
- `domain` : ADRIA Domain type, indicating geographical domain to perform site selection over.
- `mcda_vars` : Contains relevant parameters for performing guided location selection.
- `guided` : Integer indicating aggegation algorithm to use for guided location selection.

# Returns
`ranks` : n_reps * sites * 3 (last dimension indicates: site_id, seeding rank, shading rank)
    containing ranks for single scenario.
"""
function _site_selection(domain::Domain, 
    mcda_vars::DMCDA_vars, 
    guided::Int64)
    
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
    run_site_selection(domain::Domain,scenarios::DataFrame, sum_cover::NamedDimsArray, area_to_seed::Float64, aggregation_function::Function, 
        iv_type::Union{String,Int64};target_seed_sites=nothing, target_shade_sites=nothing)

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
-`ranks_store` : number of scenarios * sites * 3 (last dimension indicates: site_id, seed rank, shade rank)
    containing ranks for each scenario run.
"""
function run_site_selection(domain::Domain, 
    scenarios::DataFrame, 
    sum_cover::NamedDimsArray, 
    area_to_seed::Float64;
    target_seed_sites=nothing, 
    target_shade_sites=nothing)

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
    # Set filtered locations as n_locs+1 for consistency with time dependent ranks
    ranks_store[ranks_store.==0.0] .= length(domain.site_ids)+1

    return ranks_store

end
function run_site_selection(domain::Domain, 
    scenarios::DataFrame, 
    sum_cover::NamedDimsArray, 
    area_to_seed::Float64, 
    aggregation_function::Function, 
    iv_type::String;
    target_seed_sites=nothing, 
    target_shade_sites=nothing)

    ranks = run_site_selection(domain, 
        scenarios, 
        sum_cover, 
        area_to_seed; 
        target_seed_sites=target_seed_sites,
        target_shade_sites=target_shade_sites)
    return aggregation_function(ranks, iv_type)

end

"""
    ranks_to_location_order(ranks::NamedDimsArray)

Post-processing function for location ranks output of `run_location_selection()`. Gives the order 
of location preference for each scenario as location ids.

# Arguments
- `ranks` : Contains location ranks for each scenario of location selection, as created by 
    `run_location_selection()`.

# Returns
Location order after ranking as string IDs for each scenario.
"""
function ranks_to_location_order(ranks::NamedDimsArray)

    n_scens = length(ranks.scenarios)
    location_orders = NamedDimsArray(repeat([""], length(ranks.scenarios), length(ranks.sites)), scenarios=ranks.scenarios, ranks=1:length(ranks.sites))

    for scen in 1:n_scens
        location_orders[scenarios=scen, ranks=1:sum(ranks[scenarios=scen] .!= 0.0)] .= sort(Int.(ranks[scenarios=scen][ranks[scenarios=scen].!=0.0])).sites
    end
    return location_orders
end


"""
    ranks_to_frequencies(ranks::NamedDimsArray, n_ranks::Int64)ranks_to_frequencies(rs::ResultSet, iv_type::String; n_ranks=length(ranks.sites))
    ranks_to_frequencies(ranks::NamedDimsArray{D,T,3,A};n_ranks=length(ranks.sites),agg_func=x -> dropdims(sum(x; dims=:timesteps); dims=:timesteps),) where {D,T,A}
    ranks_to_frequencies(ranks::NamedDimsArray{D,T,2,A};n_ranks=length(ranks.sites),agg_func=nothing) where {D,T,A}

Post-processing function for location ranks output of `run_location_selection()`. Gives the frequency 
with which each location was selected at each rank across the location selection scenarios.

# Arguments
- `ranks` : Contains location ranks for each scenario of location selection, as created by 
    `run_location_selection()`.
- `n_ranks` : number of rankings (defualt is number of locations).
- `agg_func` : Aggregation function to appy after frequencies are calculated.

# Returns
Frequency with which each location was selected for each rank.
"""
function ranks_to_frequencies(ranks::NamedDimsArray, n_ranks::Int64)
    dn = NamedDims.dimnames(ranks)
    freq_dims = [n for n in dn if n != :scenarios]
    dn_subset = vcat(freq_dims, [:ranks])
    freq_elements = vcat(
        [1:size(ranks, n) for n in dn if n != :scenarios],
        [1:size(ranks, :sites)],
    )
    mn = ([size(ranks, k) for k in freq_dims]..., size(ranks, :sites))

    rank_frequencies = NamedDimsArray(zeros(mn...); zip(dn_subset, freq_elements)...)

    for rank in 1:n_ranks
        rank_frequencies[ranks=Int64(rank)] .= sum(ranks .== rank; dims=:scenarios)[scenarios=1]
    end

    return rank_frequencies
end
function ranks_to_frequencies(
    ranks::NamedDimsArray{D,T,3,A};
    n_ranks=length(ranks.sites),
    agg_func=x -> dropdims(sum(x; dims=:timesteps); dims=:timesteps),
) where {D,T,A}
    return agg_func(ranks_to_frequencies(ranks, n_ranks))
end
function ranks_to_frequencies(
    ranks::NamedDimsArray{D,T,2,A};
    n_ranks=length(ranks.sites),
    agg_func=nothing,
) where {D,T,A}
    if !isnothing(agg_func)
        return agg_func(ranks_to_frequencies(ranks, n_ranks))
    end

    return ranks_to_frequencies(ranks, n_ranks)
end

"""
    location_selection_frequencies(ranks::NamedDimsArray;n_loc_int::Int64=5)
    location_selection_frequencies(inv_log::NamedDimsArray{D,T,4,A};dims::Union{Symbol,Vector{Symbol}}=:coral_id) where {D,T,A}

Post-processing function for intervention logs. Calculates the frequencies with which locations were selected for a particular intervention,
for a selection of scenarios (e.g. selected robust scenarios).

# Arguments
- `ranks` : Contains location ranks for each scenario of location selection, as created by 
    `run_location_selection()`.
- `n_loc_int` : number of locations which are intervened at for each intervention decision.
- `dims` : dimensions to sum selection frequencies over.

# Returns 
Counts for location selection at each location in the domain.
"""
function location_selection_frequencies(
    ranks::NamedDimsArray;
    n_loc_int::Int64=5,
)
    ranks_frequencies = ranks_to_frequencies(ranks; n_ranks=n_loc_int)
    loc_count = sum(ranks_frequencies[ranks=1:n_loc_int], dims=2)[ranks=1]

    return loc_count
end
function location_selection_frequencies(
    inv_log::NamedDimsArray{D,T,4,A};
    dims::Union{Symbol,Vector{Symbol}}=:coral_id,
) where {D,T,A}
    loc_count = dropdims(
        sum(dropdims(sum(inv_log; dims=dims); dims=dims) .> 0; dims=:scenarios);
        dims=:scenarios,
    )
    return loc_count
end

"""
    summed_inverse_rank(ranks::NamedDimsArray{D,T,3,A};dims::Union{Symbol,Vector{Symbol}}=[:scenarios, :timsteps],
            ) where {D,T,A}
    summed_inverse_rank(ranks::NamedDimsArray{D,T,2,A};) where {D,T,A}
    summed_inverse_rank(ranks::NamedDimsArray,dims::Union{Symbol,Vector{Symbol}},)

Calculates (number of sites) .- ranks summed over the dimension dims and transformed using agg_func 
    (default no transformation).

# Arguments
- `ranks` : Contains location ranks for each scenario of location selection, as created by 
    `run_location_selection()`.
- `dims` : Dimensions to sum over.

# Returns 
Inverse rankings (i.e. the greater the number the higher ranked the site).
"""
function summed_inverse_rank(
    ranks::NamedDimsArray{D,T,3,A};
    dims::Union{Symbol,Vector{Symbol}}=[:scenarios, :timsteps],
) where {D,T,A}
    return summed_inverse_rank(ranks, dims)
end
function summed_inverse_rank(
    ranks::NamedDimsArray{D,T,2,A};
) where {D,T,A}
    return summed_inverse_rank(ranks, dims)
end
function summed_inverse_rank(
    ranks::NamedDimsArray,
    dims::Union{Symbol,Vector{Symbol}},
)
    n_ranks = maximum(ranks)
    inv_ranks = dropdims(sum(n_ranks .- ranks; dims=dims); dims=dims)
    return inv_ranks ./ n_ranks
end
