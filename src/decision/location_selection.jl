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
    switch_RCPs!,
    DataCube,
    copy_datacube


"""
    rank_locations(
        dom::Domain,
        n_corals::Int64,
        scenarios::DataFrame;
        rcp=nothing,
        min_iv_locs=nothing,
        max_members=nothing,
        target_seed_locs=nothing,
        target_fog_locs=nothing
    )::YAXArray

Determine location ranks for a given domain and intervention scenarios.
Locations are ranked by order of their determined deployment suitability according to
criteria weights. Values of 1 indicate highest rank.

Interventions are assessed separately such that seeding locations are not influenced by
fogging locations.

# Arguments
- `domain` : Domain dataset to assess
- `n_corals` : The total number of corals to deploy
- `scenarios` : Scenario specification
- `rcp` : RCP conditions to assess
- `min_iv_locs` : Minimum number of locations to intervene
- `max_members` : Maximum number of locations per cluster
- `target_seed_locs` : Locations to prioritize for seeding. Currently does nothing.
- `target_fog_locs` : Locations to prioritize for fogging. Currently does nothing.

# Returns
YAXArray[n_locations ⋅ [:seed, :fog] ⋅ n_scenarios], ranks from 1 to number of locations + 1
Values equal to number of locations + 1 indicate locations that unranked.
"""
function rank_locations(
    dom::Domain,
    n_corals::Int64,
    scenarios::DataFrame;
    rcp=nothing,
    min_iv_locs=nothing,
    max_members=nothing,
    target_seed_locs=nothing,
    target_fog_locs=nothing
)::YAXArray
    n_locs = n_locations(dom)
    k_area_locs = site_k_area(dom)

    if isnothing(min_iv_locs)
        min_iv_locs = scenarios.min_iv_locations
    else
        min_iv_locs = fill(min_iv_locs, nrow(scenarios))
    end

    if isnothing(max_members)
        max_members = scenarios.cluster_max_member
    else
        max_members = fill(max_members, nrow(scenarios))
    end

    if !isnothing(rcp)
        dom = switch_RCPs!(dom, rcp)
    end

    # Set filtered locations as n_locs+1 for consistency with time dependent ranks
    ranks_store = DataCube(
        fill(0.0, n_locs, 2, nrow(scenarios));
        locations=dom.site_ids,
        intervention=[:seed, :fog],
        scenarios=1:nrow(scenarios)
    )

    target_loc_ids = Int64[]
    if !isnothing(target_seed_locs)
        append!(target_loc_ids, target_seed_locs)
    end

    if !isnothing(target_fog_locs)
        append!(target_loc_ids, target_fog_locs)
    end

    if isnothing(target_seed_locs) && isnothing(target_fog_locs)
        target_loc_ids = dom.site_ids
    end

    # Sum of coral cover (relative to k area) at each location and scenario
    sum_cover = vec(sum(dom.init_coral_cover; dims=1).data)

    leftover_space_scens = relative_leftover_space(sum_cover) .* k_area_locs'

    area_weighted_conn = dom.conn.data .* site_k_area(dom)
    conn_cache = similar(area_weighted_conn)

    in_conn, out_conn, strong_pred = connectivity_strength(area_weighted_conn, sum_cover, conn_cache)

    scens = DataCube(
        Matrix(scenarios);
        scenarios=1:nrow(scenarios),
        factors=names(scenarios)
    )

    seed_pref = SeedPreferences(dom, scens[1, :])
    fog_pref = FogPreferences(dom, scens[1, :])

    α = 0.99

    site_data = dom.site_data
    coral_habitable_locs = site_data.k .> 0.0
    for (scen_idx, scen) in enumerate(eachrow(scens))
        # Decisions should place more weight on environmental conditions
        # closer to the decision point
        plan_horizon = Int64(scen[At("plan_horizon")])
        decay = α .^ (1:plan_horizon+1).^2

        min_depth = scen[factors=At("depth_min")].data[1]
        depth_offset = scen[factors=At("depth_offset")].data[1]

        depth_criteria = identify_within_depth_bounds(site_data.depth_med, min_depth, depth_offset)
        valid_locs = coral_habitable_locs .& depth_criteria
        considered_locs = findall(valid_locs)

        MCDA_approach = mcda_methods()[Int64(scen[factors=At("guided")][1])]
        leftover_space_m² = vec(leftover_space_scens[scen_idx, :])

        corals = to_coral_spec(scenarios[scen_idx, :])
        area_to_seed = mean(n_corals * colony_mean_area(corals.mean_colony_diameter_m[corals.class_id.==2]))

        seed_pref = SeedPreferences(dom, scen)
        fog_pref = FogPreferences(dom, scen)

        # Determine environmental projections
        dhw_scen_idx = Int64(scen[factors=At("dhw_scenario")][1])
        wave_scen_idx = Int64(scen[factors=At("wave_scenario")][1])
        dhw_scens = dom.dhw_scens[:, :, dhw_scen_idx]
        wave_scens = dom.wave_scens[:, :, wave_scen_idx]

        dhw_projection = weighted_projection(dhw_scens, 1, plan_horizon, decay, 75)
        wave_projection = weighted_projection(wave_scens, 1, plan_horizon, decay, 75)

        # Create shared decision matrix
        # Ignore locations that cannot support corals or are out of depth bounds
        # from consideration
        decision_mat = decision_matrix(
            dom.site_ids[valid_locs],
            seed_pref.names;
            depth=site_data.depth_med[valid_locs],
            in_connectivity=in_conn[valid_locs],
            out_connectivity=out_conn[valid_locs],
            heat_stress=dhw_projection[valid_locs],
            wave_stress=wave_projection[valid_locs],
            coral_cover=sum_cover[valid_locs]
        )

        # Ensure what to do with this because it is usually empty
        # seed_zone = strong_pred[valid_locs]

        min_locs = min_iv_locs[scen_idx]
        selected_seed_ranks = select_locations(
            seed_pref,
            decision_mat,
            MCDA_approach,
            site_data.cluster_id,
            area_to_seed,
            considered_locs,
            leftover_space_m²,
            min_locs,
            max_members[scen_idx]
        )

        if !isempty(selected_seed_ranks)
            ranks_store[locations=At(selected_seed_ranks), intervention=At(:seed), scenarios=scen_idx] .= 1:length(selected_seed_ranks)
        end

        selected_fog_ranks = select_locations(
            fog_pref, decision_mat, MCDA_approach, min_locs
        )
        if !isempty(selected_fog_ranks)
            ranks_store[locations=At(selected_fog_ranks), intervention=At(:fog), scenarios=scen_idx] .= 1:min_locs
        end
    end

    return ranks_store
end

"""
    aggregate_location_ranks(
        dom::Domain,
        n_corals::Int64,
        scenarios::DataFrame,
        agg_func::Function,
        iv_type::Symbol;
        rcp=nothing,
        min_iv_locs=nothing,
        max_members=nothing,
        target_seed_locs=nothing,
        target_fog_locs=nothing
    )::YAXArray

Return aggregate location ranks for a given domain and scenarios.

# Arguments
- `domain` : Domain dataset to assess
- `n_corals` : The total number of corals to deploy
- `scenarios` : Scenario specification
- `agg_func` : Aggregation function to apply, e.g `ranks_to_frequencies` or
    `ranks_to_location_order`
- `iv_type` : Intervention type to assess (`:seed` or `:fog`)
- `rcp` : RCP conditions to assess
- `min_iv_locs` : Minimum number of locations to intervene
- `max_members` : Maximum number of locations per cluster
- `target_seed_locs` : Locations to prioritize for seeding. Currently does nothing.
- `target_fog_locs` : Locations to prioritize for fogging. Currently does nothing.

# Returns
YAXArray[n_locations ⋅ [:seed, :fog]], with ranks aggregated according to `agg_func`.
"""
function aggregate_location_ranks(
    domain::Domain,
    n_corals::Int64,
    scenarios::DataFrame,
    agg_func::Function,
    iv_type::Symbol;
    rcp=nothing,
    target_seed_locs=nothing,
    target_fog_locs=nothing,
)::YAXArray
    throw("Unimplemented")
    ranks = rank_locations(
        domain,
        n_corals,
        scenarios;
        rcp=rcp,
        target_seed_locs=target_seed_locs,
        target_fog_locs=target_fog_locs,
    )
    return agg_func(ranks[intervention=At(iv_type)])
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

"""
    selection_score(ranks::YAXArray{T,3}; iv_type::Symbol)
    selection_score(ranks::YAXArray{T,4}; iv_type::Symbol; keep_time=false)

Calculates score ∈ [0, 1] indicative of the relative desirability of each location.
The score reflects the location ranking and frequency of attaining a high rank.

Scores of 1 indicates the location was always selected.

# Example
```julia
ranks = ADRIA.decision.rank_locations(dom, n_corals, scens)
# ranks will be a 3-dimensional array of [locations ⋅ interventions ⋅ scenarios]

# Identify locations that were most desirable for seeding over all scenarios
seed_scores = ADRIA.decision.selection_score(ranks, :seed)

# Identify locations that were most desirable for fogging over all scenarios
fog_scores = ADRIA.decision.selection_score(ranks, :fog)

# Selection scores can be assessed for a subset of scenarios, including a specific scenario
ADRIA.decision.selection_score(ranks[scenarios=1:4], :seed)

# Analysis includes the time dimension by default where scenario runs are being assessed
# such that the score for each location is returned
rs = ADRIA.scenario_runs(dom, scens, "45")
# rs.ranks will be a 4-dimensional array of
# [timesteps ⋅ locations ⋅ interventions ⋅ scenarios]

ADRIA.decision.selection_score(rs.ranks, 1; keep_time=false)
# 216-element YAXArray{Float32,1} with dimensions:
#   Dim{:sites} Sampled{Int64} 1:216 ForwardOrdered Regular Points
# Total size: 864.0 bytes

# Change to true to keep the time dimension (to get the score per time step)
ADRIA.decision.selection_score(rs.ranks, 1; keep_time=true)
# 75×216 YAXArray{Float32,2} with dimensions:
#   Dim{:timesteps} Sampled{Int64} 1:75 ForwardOrdered Regular Points,
#   Dim{:sites} Sampled{Int64} 1:216 ForwardOrdered Regular Points
# Total size: 63.28 KB
```

# Arguments
- `ranks` : Rankings of locations from `rank_locations()`
- `iv_type` : The intervention type to assess

# Returns
Selection score
"""
function selection_score(
    ranks::YAXArray{T, 3},
    iv_type::Union{Symbol,Int64},
)::YAXArray where {T<:Union{Int64, Float32, Float64}}
    # 1 is best rank, n_locs is worst rank, 0 are locations that were ignored
    # Determine the lowest rank for each scenario
    lowest_ranks = maximum([maximum(r)
                    for r in eachcol(ranks[intervention=At(iv_type)])])

    return _calc_selection_score(ranks, lowest_ranks, iv_type, (:scenarios, ))
end
function selection_score(
    ranks::YAXArray{T, 4},
    iv_type::Union{Symbol,Int64};
    keep_time=false
)::YAXArray where {T<:Union{Int64, Float32, Float64}}
    # [timesteps ⋅ locations ⋅ interventions ⋅ scenarios]
    # 1 is best rank, n_locs is worst rank, 0 values indicate locations that were ignored
    lowest_rank = maximum(ranks)

    # Determie dimensions to squash
    dims = keep_time ? (:scenarios, ) : (:scenarios, :timesteps)

    selection_score = _calc_selection_score(ranks, lowest_rank, iv_type, dims)

    return selection_score
end

"""
    _calc_selection_score(ranks::YAXArray, lowest_rank::Union{Int64,Float64,Float32}, iv_type::Union{Symbol,Int64}, dims::Tuple)::YAXArray

Determine the selection score.
Note: If `timesteps` are to be squashed the scores are normalized against the maximum score.

# Arguments
- `ranks` : The recorded/logged rankings
- `lowest_rank` : The identified lowest rank
- `iv_type` : The intervention type to assess (`:seed` or `:fog`)
- `dims` : Dimensions to squash

# Returns
YAXArray
"""
function _calc_selection_score(ranks::YAXArray, lowest_rank::Union{Int64,Float64,Float32}, iv_type::Union{Symbol,Int64}, dims::Tuple)::YAXArray
    # Subtract 1 from rank:
    # If the best rank is 1 and lowest is 3 (out of a single scenario):
    #     ([lowest rank] - [rank]) / [lowest rank]
    #     (3 - 1) / 3 = 0.6666
    # If we subtract 1 from the rank:
    #     (3 - 0) / 3 = 1.0
    #     (3 - 1) / 3 = 0.666
    #     (3 - 2) / 3 = 0.333
    # and unconsidered locations are given a score of 0.
    # Where the number of locations are > lowest considered rank,
    # the maximum of the rank and 0.0 is taken, so that unconsidered locations are still
    # given zero rank.
    # If, for example, there are 100 locations, and the lowest considered rank is 3:
    #     max((3 - 100) / 3, 0.0)

    tsliced = mapslices(x -> any(x .> 0) ? lowest_rank .- (x .- 1.0) : 0.0, ranks[intervention=At(iv_type)], dims="timesteps")
    selection_score = dropdims(
        sum(tsliced, dims=dims); dims=dims
    )

    times_ranked = size(ranks, :scenarios)

    # Determine number of decision years if timesteps are to be kept.
    # There's probably a better way of doing this but this works...
    if :timesteps in dims
        decisions = 0.0
        for s in axes(tsliced, :scenarios)
            for t in axes(tsliced, :timesteps)
                if any(tsliced[timesteps=At(t), scenarios=At(s)] .> 0.0)
                    decisions += 1.0
                end
            end
        end

        times_ranked = times_ranked * decisions
    end

    selection_score = max.(selection_score ./ (lowest_rank * times_ranked), 0.0)
    if :timesteps in dims
        selection_score ./= maximum(selection_score)
    end

    return selection_score
end

"""
    selection_frequency(ranks::YAXArray{T, 3}, iv_type::Union{Symbol,Int64})
    selection_frequency(ranks::YAXArray{T, 4}, iv_type::Union{Symbol,Int64})

Determine normalized frequency at which a location is selected for intervention, where
0 indicates the location is not selected at all, and 1 indicates the selection is
selected across all time and under all conditions modelled.

# Example
```julia
ranks = ADRIA.decision.rank_locations(dom, n_corals, scens)

# Identify locations that were selected for seeding over all scenarios
seed_freq = ADRIA.decision.selection_frequency(ranks, :seed)

# Selection scores can be assessed for a subset of scenarios, including a specific scenario
ADRIA.decision.selection_frequency(ranks[scenarios=1:4], :seed)

# Assess selection frequency for scenario runs
# Note: indices have to be used for now
rs = ADRIA.scenario_runs(dom, 16, "45")
scen_seed_freq = ADRIA.decision.selection_frequency(rs.ranks, 1)
```

# Arguments
- `ranks` : Ranks to assess
- `iv_type` : Intervention type to assess

# Return
Normalized selection frequency
"""
function selection_frequency(
    ranks::YAXArray{T, 3},
    iv_type::Union{Symbol,Int64}
)::YAXArray where {T<:Union{Int64, Float32, Float64}}
    # 0 is unconsidered locations, values > 0 indicate some rank was assigned.
    s = copy_datacube(ranks[intervention=At(iv_type)])
    s[s .> 0.0] .= 1.0

    selection_freq = dropdims(
        sum(s, dims=:scenarios); dims=:scenarios
    )

    norm_freq = selection_freq ./ maximum(selection_freq)
    norm_freq[norm_freq .> 0.0] .= 1.0 .- norm_freq[norm_freq .> 0.0]

    return norm_freq
end
function selection_frequency(
    ranks::YAXArray{T, 4},
    iv_type::Union{Symbol,Int64}
)::YAXArray where {T<:Union{Int64, Float32, Float64}}
    # 0 is unconsidered locations, values > 0 indicate some rank was assigned.
    s = copy_datacube(ranks[intervention=At(iv_type)])
    s[s .> 0.0] .= 1.0

    selection_freq = dropdims(
        sum(s, dims=(:scenarios, :timesteps)); dims=(:scenarios, :timesteps)
    )

    norm_freq = selection_freq ./ maximum(selection_freq)
    norm_freq[norm_freq .> 0.0] .= 1.0 .- norm_freq[norm_freq .> 0.0]

    return norm_freq
end