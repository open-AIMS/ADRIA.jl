using NamedDims, AxisKeys

using ADRIA: connectivity_strength, relative_leftover_space, site_k_area

"""
	_location_selection(domain::Domain, sum_cover::AbstractArray, mcda_vars::DMCDA_vars, guided::Int64)::Matrix

Select locations for a given domain and criteria/weightings/thresholds, using a chosen
MCDA method.

# Arguments
- `domain` : The geospatial domain to assess
- `sum_cover` :  Relative coral cover at each location
- `mcda_vars` : Parameters for MCDA
- `guided` : ID of MCDA algorithm to apply

# Returns
Matrix[n_sites ⋅ 2], where columns hold seeding and shading ranks for each location.
"""
function _location_selection(
	domain::Domain, sum_cover::AbstractArray, mcda_vars::DMCDA_vars, guided::Int64
)::Matrix
	site_ids = mcda_vars.site_ids
	n_sites = length(site_ids)

	# site_id, seeding rank, shading rank
	rankingsin = [mcda_vars.site_ids zeros(Int64, n_sites) zeros(Int64, n_sites)]

	pref_seed_sites::Vector{Int64} = zeros(Int64, mcda_vars.n_site_int)
	pref_fog_sites::Vector{Int64} = zeros(Int64, mcda_vars.n_site_int)

	# Determine connectivity strength
	# Account for cases where no coral cover
	in_conn, out_conn, strong_pred = connectivity_strength(
		domain.TP_data .* site_k_area(domain),
		vec(sum_cover),
		similar(domain.TP_data)  # tmp cache matrix
	)

	# Perform location selection for seeding and shading.
	seed_true, fog_true = true, true

	(_, _, ranks) = guided_site_selection(
		mcda_vars,
		guided,
		seed_true,
		fog_true,
		pref_seed_sites,
		pref_fog_sites,
		rankingsin,
		in_conn[site_ids],
		out_conn[site_ids],
		strong_pred[site_ids],
	)

	return ranks[:, 2:3]
end

"""
	rank_locations(domain::Domain, scenarios::DataFrame, sum_cover::NamedDimsArray, area_to_seed::Float64; target_seed_sites=nothing, target_fog_sites=nothing)::NamedDimsArray
	rank_locations(domain::Domain,scenarios::DataFrame, sum_cover::NamedDimsArray, area_to_seed::Float64, agg_func::Function,
		iv_type::Union{String,Int64}; target_seed_sites=nothing, target_fog_sites=nothing)::AbstractArray

Return location ranks for a given domain and scenarios.

# Arguments
- `domain` : The geospatial domain locations were selected from
- `scenarios` : Scenario specification
- `sum_cover` : Matrix[n_scenarios ⋅ n_sites] containing the total relative coral cover at each
	location, for each scenario
- `area_to_seed` : Area of coral to be seeded at each time step in km²
- `agg_func` : Aggregation function to apply, e.g `ranks_to_frequencies` or
	`ranks_to_location_order`
- `iv_type` : ID of intervention (1 = seeding, 2 = fogging)
- `target_seed_sites` : list of candidate locations for seeding (indices)
- `target_fog_sites` : list of candidate location to fog (indices)

# Returns
Array[n_locations ⋅ 2 ⋅ n_scenarios], where columns hold seeding and shading ranks.
"""
function rank_locations(
	domain::Domain,
	scenarios::DataFrame,
	sum_cover::NamedDimsArray,
	area_to_seed::Float64;
	target_seed_sites = nothing,
	target_fog_sites = nothing,
)::NamedDimsArray
	n_locs = n_locations(domain)
	k_area_locs = site_k_area(domain)

	ranks_store = NamedDimsArray(
		zeros(n_locs, 2, nrow(scenarios));
		sites = 1:n_locs,
		intervention = ["seed", "fog"],
		scenarios = 1:nrow(scenarios),
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

	if !isnothing(target_fog_sites)
		append!(target_site_ids, target_fog_sites)
	end

	if isnothing(target_seed_sites) && isnothing(target_fog_sites)
		target_site_ids = collect(1:length(domain.site_ids))
	end

	leftover_space_scens = relative_leftover_space(sum_cover.data) .* k_area_locs'

	for (scen_idx, scen) in enumerate(eachrow(scenarios))
		depth_criteria = within_depth_bounds(
			domain.site_data.depth_med,
			scen.depth_min .+ scen.depth_offset,
			scen.depth_min
		)
		depth_priority = findall(depth_criteria)

		considered_sites = target_site_ids[findall(in(depth_priority), target_site_ids)]
		mcda_vars_temp = DMCDA_vars(
			domain,
			scen,
			considered_sites,
			leftover_space_scens[scen_idx, :],
			area_to_seed,
			summary_stat_env(wave_scens[:, :, target_wave_scens], (:timesteps, :scenarios)),
			summary_stat_env(dhw_scens[:, :, target_dhw_scens], (:timesteps, :scenarios)),
		)

		ranks_store(; scenarios = scen_idx, sites = considered_sites) .= _location_selection(
			domain,
			collect(sum_cover[scen_idx, :]),
			mcda_vars_temp,
			scen.guided
		)
	end

	# Set filtered locations as n_locs+1 for consistency with time dependent ranks
	ranks_store[ranks_store .== 0.0] .= length(domain.site_ids) + 1

	return ranks_store
end
function rank_locations(
	domain::Domain,
	scenarios::DataFrame,
	sum_cover::NamedDimsArray,
	area_to_seed::Float64,
	agg_func::Function,
	iv_type::Union{Int64, Symbol, String};
	target_seed_sites = nothing,
	target_fog_sites = nothing,
)::AbstractArray
	ranks = rank_locations(
		domain,
		scenarios,
		sum_cover,
		area_to_seed;
		target_seed_sites = target_seed_sites,
		target_fog_sites = target_fog_sites,
	)
	local iv_id
	try
		iv_id = axiskeys(ranks, :intervention)[iv_type]
	catch err
		if !(err isa ArgumentError)
			rethrow(err)
		end

		iv_id = iv_type
	end

	return agg_func(ranks(; intervention = iv_id))
end

"""
	ranks_to_frequencies(ranks::NamedDimsArray, n_ranks::Int64)
	ranks_to_frequencies(ranks::NamedDimsArray{D,T,3,A}; n_ranks=length(ranks.sites), agg_func=x -> dropdims(sum(x; dims=:timesteps); dims=:timesteps),) where {D,T,A}
	ranks_to_frequencies(ranks::NamedDimsArray{D,T,2,A}; n_ranks=length(ranks.sites), agg_func=nothing) where {D,T,A}

Returns the frequency with which each location was ranked across scenarios.
Uses the results from `rank_locations()`.

Note: By default, ranks are aggregated over time where `ranks` is a 3-dimensional array.

# Arguments
- `ranks` : Location ranks from `rank_locations()`
- `n_ranks` : Number of rankings (default is number of locations).
- `agg_func` : Aggregation function to appy after frequencies are calculated.

# Returns
Frequency with which each location was selected for each rank.
"""
function ranks_to_frequencies(ranks::NamedDimsArray, n_ranks::Int64)::NamedDimsArray
	dn = NamedDims.dimnames(ranks)
	freq_dims = [n for n in dn if n != :scenarios]
	dn_subset = vcat(freq_dims, [:ranks])
	freq_elements = vcat(
		[1:size(ranks, n) for n in dn if n != :scenarios],
		[1:size(ranks, :sites)]
	)
	mn = ([size(ranks, k) for k in freq_dims]..., size(ranks, :sites))

	rank_frequencies = NamedDimsArray(zeros(mn...); zip(dn_subset, freq_elements)...)

	for rank in 1:n_ranks
		rank_frequencies[ranks = Int64(rank)] .= sum(ranks .== rank; dims = :scenarios)[scenarios = 1]
	end

	return rank_frequencies
end
function ranks_to_frequencies(
	ranks::NamedDimsArray{D, T, 3, A};
	n_ranks::Int64 = length(ranks.sites),
	agg_func = nothing,
)::NamedDimsArray where {D, T, A}
	if !isnothing(agg_func)
		return agg_func(ranks_to_frequencies(ranks, n_ranks))
	end

	return ranks_to_frequencies(ranks, n_ranks)
end
function ranks_to_frequencies(
	ranks::NamedDimsArray{D, T, 2, A};
	n_ranks::Int64 = length(ranks.sites),
	agg_func = nothing,
)::NamedDimsArray where {D, T, A}
	if !isnothing(agg_func)
		return agg_func(ranks_to_frequencies(ranks, n_ranks))
	end

	return ranks_to_frequencies(ranks, n_ranks)
end

"""
	location_selection_frequencies(ranks::NamedDimsArray; n_iv_locs::Int64=5)
	location_selection_frequencies(iv_log::NamedDimsArray{D,T,4,A}; dims::Union{Symbol,Vector{Symbol}}=:coral_id) where {D,T,A}

Determines the count of times each location was selected for a specific intervention over a
set of scenarios.

# Arguments
- `ranks` : Rankings of locations `rank_locations()`
- `n_iv_locs` : number of locations intervened at, for each decision point
- `iv_log` : Intervention logs
- `dims` : dimensions to sum selection frequencies over

# Returns
Number of times each location was selected for an intervention.
"""
function location_selection_frequencies(
	ranks::NamedDimsArray;
	n_iv_locs::Int64 = 5
)::NamedDimsArray
	ranks_frequencies = ranks_to_frequencies(ranks; n_ranks = n_iv_locs)
	loc_count = sum(ranks_frequencies[ranks = 1:n_iv_locs]; dims = :ranks)[ranks = 1]

	return loc_count
end
function location_selection_frequencies(
	iv_log::NamedDimsArray{D, T, 4, A};
	dims::Union{Symbol, Vector{Symbol}} = :coral_id
)::NamedDimsArray where {D, T, A}
	loc_count = dropdims(
		sum(dropdims(sum(iv_log; dims = dims); dims = dims) .> 0; dims = :scenarios);
		dims = :scenarios,
	)
	return loc_count
end

"""Drop single dimensions."""
_drop_single(x::AbstractMatrix) = dropdims(x; dims = (findall(size(x) .== 1)...,))

"""
	selection_score(ranks::NamedDimsArray{D,T,3,A}; dims::Vector{Symbol}=[:scenarios, :timesteps]) where {D,T,A}
	selection_score(ranks::NamedDimsArray{D,T,2,A}) where {D,T,A}
	selection_score(ranks::NamedDimsArray, dims::Vector{Symbol})

Calculates score ∈ [0, 1], where 1 is the highest score possible, indicative of the relative
desirability of each location.

The score reflects the location ranking and frequency of attaining a high rank.

# Arguments
- `ranks` : Rankings of locations from `rank_locations()`
- `dims` : Dimensions to sum over

# Returns
Selection score
"""
function selection_score(
	ranks::NamedDimsArray{D, T, 3, A};
	dims::Vector{Symbol} = [:scenarios, :timesteps]
)::NamedDimsArray where {D, T, A}
	return _drop_single(selection_score(ranks, dims))
end
function selection_score(
	ranks::NamedDimsArray{D, T, 2, A};
	dims::Vector{Symbol} = [:scenarios]
)::NamedDimsArray where {D, T, A}
	return selection_score(ranks, dims)
end
function selection_score(
	ranks::NamedDimsArray,
	dims::Vector{Symbol},
)::NamedDimsArray
    lowest_rank = maximum(ranks)  # 1 is best rank, n_sites + 1 is worst rank
    selection_score = dropdims(sum(lowest_rank .- ranks; dims=dims); dims=dims[1])
    return selection_score ./ ((lowest_rank - 1) * prod([size(ranks, d) for d in dims]))
end

"""
    decision_matrices(rs::ResultSet, criteria_row::DataFrameRow;
        loc_coral_cover = rs.site_max_coral_cover::Vector{Float64}, RCP::String = "45")

Calculates a decision matrix for a specified intervention, using a scenario specification 
and Domain alone. These can be visualised spatially using `viz.decision_matrices`.

# Arguments
- `rs` : ADRIA ResultSet
- `criteria_row` :  A row of a scenario dataframe, containing intervention criteria weights.
- `loc_coral_cover` : Relative coral cover to site k area (dims: nspecies*nsites), default 
is max cover over scenarios in rs.
-  `RCP` : RCP scenario to use (waves and dhws), default is "45".

# Returns
Selection score
"""
function decision_matrices(
    rs::ResultSet,
    criteria_row::DataFrameRow;
    mcda_func::Function;
    loc_coral_cover = rs.site_max_coral_cover::Vector{Float64},
    RCP::String = "45",
    env_aggregation::Function = summary_stat_env,
    site_ids = collect(1:length(rs.site_data.site_id))
    leftover_space = relative_leftover_space(loc_coral_cover) .* site_k_area(rs)

    wave_stress = env_aggregation(
        rs.wave_stats[RCP](; stat = "mean"), (:wave_scenario)
    )
    heat_stress = env_aggregation(
        rs.dhw_stats[RCP](; stat = "mean"), (:dhw_scenario)
    )

    conn_data = rs.connectivity_data[RCP]  # connectivity matrix
    connectivity_data = connectivity_strength(
        conn_data .* site_k_area(rs), loc_coral_cover, conn_data
    )
    # Strongest larval source location for each location
    strong_pred = connectivity_data.strongest_predecessor

    # Criteria for strongest larval sources to priority locations
    priority_source_criteria = ADRIA.decision.priority_predecessor_criteria(
        strong_pred, vec(rs.sim_constants["priority_sites"]), length(site_ids)
    )
    # Criteria for strongest larval sources to/members of priority zones
    zones_criteria = ADRIA.decision.zones_criteria(
        vec(rs.sim_constants["priority_zones"]), rs.site_data.zone_type, strong_pred,
        site_ids,
    )

    spec = model_spec(rs)
    crit = component_params(spec, CriteriaWeights)
    weights_seed_crit = criteria_params(crit, "(:seed, :weight)")
    weights_fog_crit = criteria_params(crit, "(:fog, :weight)")

    A, filtered_sites = ADRIA.decision.create_decision_matrix(
        site_ids,
        connectivity_data.in_conn,
        connectivity_data.out_conn,
        leftover_space[site_ids],
        wave_stress[site_ids],
        heat_stress[site_ids],
        rs.site_data.depth_med[site_ids],
        priority_source_criteria,
        zones_criteria,
        criteria_row.deployed_coral_risk_tol,
    )

    S, wse = ADRIA.decision.create_seed_matrix(
        A,
        criteria_row.seed_in_connectivity,
        criteria_row.seed_out_connectivity,
        criteria_row.seed_wave_stress,
        criteria_row.seed_heat_stress,
        criteria_row.seed_priority,
        criteria_row.seed_zone,
        criteria_row.seed_coral_cover_low,
        criteria_row.seed_depth;
        filter_space = -1.0,
    )
    SE = NamedDimsArray(
        S[:, 2:end];
        locations = rs.site_data.site_id[Int.(S[:, 1])],
        criteria = weights_seed_crit.fieldname,
    )

    S, wsh = ADRIA.decision.create_fog_matrix(
        A,
        site_k_area(rs)[site_ids][filtered_sites],
        criteria_row.fog_in_connectivity,
        criteria_row.fog_out_connectivity,
        criteria_row.fog_wave_stress,
        criteria_row.fog_heat_stress,
        criteria_row.fog_priority,
        criteria_row.fog_zone,
        criteria_row.fog_coral_cover_high,
    )
    SH = NamedDimsArray(
        S[:, 2:end];
        locations = rs.site_data.site_id[Int.(S[:, 1])],
        criteria = weights_fog_crit.fieldname,
    )

    # Create decision info struct
    decision_dict = Dict(
        :seed_matrix => SE,
        :seed_scores => vec(mcda_func(Matrix(SE))),
        :fog_matrix => SH,
        :fog_scores => vec(mcda_func(Matrix(SH))),
    )
    return decision_dict
end
