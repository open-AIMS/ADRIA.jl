using NamedDims, AxisKeys, YAXArrays

using ADRIA:
    Domain,
    sample,
    connectivity_strength,
    relative_leftover_space,
    site_k_area,
    n_locations,
    to_coral_spec,
    colony_mean_area,
    switch_RCPs!


"""
    rank_locations(domain::Domain, scenarios::DataFrame, area_to_seed::Float64; target_seed_sites=nothing, target_fog_sites=nothing)::NamedDimsArray
    rank_locations(domain::Domain,scenarios::DataFrame, area_to_seed::Float64, agg_func::Function,
        iv_type::Union{String,Int64}; target_seed_sites=nothing, target_fog_sites=nothing)::AbstractArray

Return location ranks for a given domain and scenarios.

# Arguments
- `domain` : The geospatial domain locations were selected from
- `scenarios` : Scenario specification
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
    dom::Domain,
    n_corals::Int64,
    scenarios::DataFrame;
    rcp=nothing,
    n_iv_locs=nothing,
    target_seed_sites=nothing,
    target_fog_sites=nothing,
)::NamedDimsArray
    n_locs = n_locations(dom)
    k_area_locs = site_k_area(dom)

    if isnothing(n_iv_locs)
        n_iv_locs = dom.sim_constants.n_site_int
    end

    if !isnothing(rcp)
        dom = switch_RCPs!(dom, rcp)
    end

    ranks_store = NamedDimsArray(
        zeros(n_locs, 2, nrow(scenarios)),
        sites=1:n_locs,
        intervention=["seed", "fog"],
        scenarios=1:nrow(scenarios)
    )

    target_site_ids = Int64[]
    if !isnothing(target_seed_sites)
        append!(target_site_ids, target_seed_sites)
    end

    if !isnothing(target_fog_sites)
        append!(target_site_ids, target_fog_sites)
    end

    if isnothing(target_seed_sites) && isnothing(target_fog_sites)
        target_site_ids = dom.site_ids
    end

    # Sum of coral cover (relative to k area) at each location and scenario
    sum_cover = repeat(sum(dom.init_coral_cover; dims=1), size(scenarios, 1))

    leftover_space_scens = relative_leftover_space(sum_cover.data) .* k_area_locs'

    area_weighted_conn = dom.conn .* site_k_area(dom)
    conn_cache = similar(area_weighted_conn)
    conn_names = ["seed_in_connectivity", "seed_out_connectivity", "seed_priority"]

    in_conn, out_conn, strong_pred = connectivity_strength(area_weighted_conn, collect(sum_cover[1, :]), conn_cache)

    axlist = (
        Dim{:scenarios}(1:nrow(scenarios)),
        Dim{:factors}(names(scenarios)),
    )
    scens = YAXArray(
        axlist,
        Matrix(scenarios)
    )

    seed_pref = SeedPreferences(dom, scens[1, :])
    fog_pref = FogPreferences(dom, scens[1, :])
    site_data = dom.site_data
    coral_habitable_locs = site_data.k .> 0.0
    for (scen_idx, scen) in enumerate(eachrow(scens))
        min_depth = scen[factors=At("depth_min")].data[1]
        depth_criteria::BitArray{1} = within_depth_bounds(
            site_data.depth_med,
            min_depth .+ scen[factors=At("depth_offset")].data[1],
            min_depth
        )

        valid_locs = coral_habitable_locs .& depth_criteria

        MCDA_approach = mcda_methods()[Int64(scen[factors=At("guided")][1])]
        leftover_space_m² = vec(leftover_space_scens[scen_idx, valid_locs])

        corals = to_coral_spec(scenarios[scen_idx, :])
        area_to_seed = mean(n_corals * colony_mean_area(corals.mean_colony_diameter_m[corals.class_id.==2]))

        seed_pref = SeedPreferences(dom, scen)
        fog_pref = FogPreferences(dom, scen)

        # Create shared decision matrix
        decision_mat = decision_matrix(dom.site_ids, seed_pref.names)

        # Set criteria values that do not change between time steps
        decision_mat[criteria=At("seed_depth")] = site_data.depth_med
        # Ensure what to do with this because it is usually empty
        # decision_mat[criteria=At("seed_zone")]

        # Remove locations that cannot support corals or are out of depth bounds
        # from consideration
        decision_mat = decision_mat[valid_locs, :]

        # Number of time steps in environmental layers to look ahead when making decisions
        Main.@infiltrate
        horizon::UnitRange{Int64} = 1:1+Int64(scen[factors=At("plan_horizon")][1])
        d_s::UnitRange{Int64} = 1:length(horizon)

        dhw_scen_idx = Int64(scen[factors=At("dhw_scenario")][1])
        wave_scen_idx = Int64(scen[factors=At("wave_scenario")][1])
        dhw_scens = dom.dhw_scens[:, :, dhw_scen_idx]
        wave_scens = dom.wave_scens[:, :, wave_scen_idx]

        @views env_horizon = decay[d_s] .* dhw_scen[horizon, considered_locs]
        decision_mat[criteria=At("seed_heat_stress")] = summary_stat_env(dhw_scens[env_horizon, valid_locs, dhw_scen_idx], :timesteps)

        @views env_horizon = decay[d_s] .* wave_scen[horizon, considered_locs]
        decision_mat[criteria=At("seed_wave_stress")] = summary_stat_env(wave_scens[env_horizon, valid_locs, wave_scen_idx], :timesteps)

        decision_mat[criteria=At("seed_coral_cover")] = collect(sum_cover[scen_idx, valid_locs])

        decision_mat[criteria=At(conn_names)] .= [
            in_conn[valid_locs] out_conn[valid_locs] strong_pred[valid_locs]
        ]

        # Identify valid, non-constant, columns for use in MCDA
        is_const = Bool[length(x) == 1 for x in unique.(eachcol(decision_mat.data))]
        valid_prefs = seed_pref.names[.!is_const]
        decision_mat = decision_mat[criteria=At(valid_prefs)]

        # Recreate preferences, removing criteria that are constant for this timestep
        sp = SeedPreferences(valid_prefs, seed_pref.weights[.!is_const], seed_pref.directions[.!is_const])
        fp = FogPreferences(valid_prefs, fog_pref.weights[.!is_const], fog_pref.directions[.!is_const])

        selected_seed_ranks = select_locations(
            sp,
            decision_mat,
            MCDA_approach,
            site_data.cluster_id[valid_locs],
            area_to_seed,
            leftover_space_m²,
            n_iv_locs,
            dom.sim_constants.max_members
        )
        if !isempty(selected_seed_ranks)
            ranks_store[selected_seed_ranks[:, 2], 1, scen_idx] .= 1:size(selected_seed_ranks, 1)
        end

        selected_fog_ranks = select_locations(
            fp, decision_mat, MCDA_approach, n_iv_locs
        )
        if !isempty(selected_fog_ranks)
            ranks_store[selected_fog_ranks[:, 2], 2, scen_idx] .= 1:n_iv_locs
        end
    end

    # Set filtered locations as n_locs+1 for consistency with time dependent ranks
    ranks_store[ranks_store.==0.0] .= length(dom.site_ids)+1

    return ranks_store
end
function rank_locations(
    domain::Domain,
    scenarios::DataFrame,
    sum_cover::NamedDimsArray,
    area_to_seed::Float64,
    agg_func::Function,
    iv_type::Union{Int64, Symbol, String};
    target_seed_sites=nothing,
    target_fog_sites=nothing,
)::AbstractArray
    ranks = rank_locations(
        domain,
        scenarios,
        sum_cover,
        area_to_seed;
        target_seed_sites=target_seed_sites,
        target_fog_sites=target_fog_sites,
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

    return agg_func(ranks(intervention=iv_id))
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
    n_ranks::Int64=length(ranks.sites),
    agg_func=nothing,
)::NamedDimsArray where {D,T,A}
    if !isnothing(agg_func)
        return agg_func(ranks_to_frequencies(ranks, n_ranks))
    end

    return ranks_to_frequencies(ranks, n_ranks)
end
function ranks_to_frequencies(
    ranks::NamedDimsArray{D,T,2,A};
    n_ranks::Int64=length(ranks.sites),
    agg_func=nothing,
)::NamedDimsArray where {D,T,A}
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
    n_iv_locs::Int64=5,
)::NamedDimsArray
    ranks_frequencies = ranks_to_frequencies(ranks; n_ranks=n_iv_locs)
    loc_count = sum(ranks_frequencies[ranks=1:n_iv_locs]; dims=:ranks)[ranks=1]

    return loc_count
end
function location_selection_frequencies(
    iv_log::NamedDimsArray{D,T,4,A};
    dims::Union{Symbol,Vector{Symbol}}=:coral_id,
)::NamedDimsArray where {D,T,A}
    loc_count = dropdims(
        sum(dropdims(sum(iv_log; dims=dims); dims=dims) .> 0; dims=:scenarios);
        dims=:scenarios,
    )
    return loc_count
end


"""Drop single dimensions."""
_drop_single(x::AbstractMatrix) = dropdims(x, dims=(findall(size(x) .== 1)...,))

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
    ranks::NamedDimsArray{D,T,3,A};
    dims::Vector{Symbol}=[:scenarios, :timesteps],
)::NamedDimsArray where {D,T,A}
    return _drop_single(selection_score(ranks, dims))
end
function selection_score(
    ranks::NamedDimsArray{D,T,2,A};
    dims::Vector{Symbol}=[:scenarios],
)::NamedDimsArray where {D,T,A}
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
