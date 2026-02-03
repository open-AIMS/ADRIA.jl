"""
    _get_ranks(rs::ResultSet, intervention::Int64; kwargs...)

Extracts results for a specific intervention (:seed or :fog)
"""
function _get_ranks(rs::ResultSet, intervention::Symbol; kwargs...)
    return slice_results(rs.ranks[intervention = At(intervention)]; kwargs...)
end

"""
    _collate_ranks(rs, selected)

Collates ranks into seed/fog ranking results into a common structure.
"""
function _collate_ranks(rs::ResultSet, selected; kwargs...)::YAXArray
    n_steps, n_locs = size(selected)

    ts = timesteps(rs)
    @assert length(ts) == n_steps

    r_ids = rs.loc_data.reef_siteid
    if haskey(kwargs, :sites)
        r_ids = r_ids[kwargs[:sites]]
    end

    if length(r_ids) != n_locs
        @warn "Length of reef ids do not match number of sites"
    end

    return DataCube(
        selected;
        timesteps=ts,
        locations=r_ids,
        scenarios=1:size(selected, 3)
    )
end

"""
    seed_ranks(rs::ResultSet; kwargs...)

# Arguments
- rs : ResultSet
- kwargs : named dimensions to slice across

# Returns
YAXArray[timesteps, sites, scenarios]

# Example
```julia
ADRIA.metrics.seed_ranks(rs; timesteps=1:10, scenarios=3:5)
```
"""
function seed_ranks(rs::ResultSet; kwargs...)
    selected = _get_ranks(rs, :seed; kwargs...)
    return _collate_ranks(rs, selected; kwargs...)
end

"""
    fog_ranks(rs::ResultSet; kwargs...)

# Arguments
- rs : ResultSet
- kwargs : named dimensions to slice across

# Returns
YAXArray[timesteps, sites, scenarios]

# Example
```julia
ADRIA.metrics.fog_ranks(rs; timesteps=1:10, scenarios=3:5)
```
"""
function fog_ranks(rs::ResultSet; kwargs...)
    selected = _get_ranks(rs, :fog; kwargs...)
    return _collate_ranks(rs, selected; kwargs...)
end

"""
    _collate_ranked_locs(data::YAXArray)::Matrix{Int64}

Collates number of ranked locations.
"""
function _collate_ranked_locs(data::YAXArray)::Matrix{Int64}
    locs = zeros(Int64, size.([data], (1, 3))...)

    Threads.@threads for scen in axes(data, :scenarios)
        scen_ranks = data[:, :, scen]

        if all(scen_ranks .== 0.0)
            continue
        end

        for t in axes(scen_ranks, :timesteps)
            locs[t, scen] = count(scen_ranks[t, :] .> 0.0)
        end
    end

    return locs
end

"""
    n_seed_locations(rs::ResultSet; kwargs...)::Matrix{Int64}

Determine the number of locations seeded at each time step, for each scenario.

# Returns
YAXArray[timesteps ⋅ scenarios] indicating the number of locations seeded at each time step.
"""
function n_seed_locations(rs::ResultSet; kwargs...)::YAXArray{Int64}
    ranked_locs = seed_ranks(rs; kwargs...)

    return DataCube(
        _collate_ranked_locs(ranked_locs);
        timesteps=collect(ranked_locs.timesteps),
        scenarios=collect(ranked_locs.scenarios)
    )
end

"""
    n_fog_locations(rs::ResultSet; kwargs...)::Matrix{Int64}

Determine the number of locations fogged at each time step, for each scenario.

# Returns
YAXArray[timesteps ⋅ scenarios] indicating the number of locations fogged at each time step.
"""
function n_fog_locations(rs::ResultSet; kwargs...)::YAXArray{Int64}
    ranked_locs = fog_ranks(rs; kwargs...)

    return DataCube(
        _collate_ranked_locs(ranked_locs);
        timesteps=collect(ranked_locs.timesteps),
        scenarios=collect(ranked_locs.scenarios)
    )
end

"""
    top_n_seeded_sites(rs::ResultSet, n::Int64; kwargs...)

Get the top n seeded sites over time by their unique location id.
Lower rank values are better (e.g., 1 = first choice)

# Arguments
- rs : ResultSet
- n : `n` locations to retrieve
- kwargs : dimensions to slice across

# Returns
YAXArray[locations, [loc_id, loc_name, rank], scenarios]
"""
function top_n_seeded_sites(rs::ResultSet, n::Int64; kwargs...)::YAXArray
    ranked_locs = seed_ranks(rs; kwargs...)

    r_ids = rs.loc_data.reef_siteid
    min_rank = length(r_ids) + 1

    c_ranks = collect(dropdims(mean(ranked_locs; dims=1); dims=1))

    n_scenarios = size(ranked_locs, 3)
    top_sites = Array{Union{String,Int32,Float32,Missing}}(undef, n, 3, n_scenarios)
    for scen in axes(ranked_locs, 3)
        flat = vec(c_ranks[:, scen])
        flat[flat .== 0.0] .= min_rank

        idx = collect(partialsortperm(flat, 1:n))

        rank_score = flat[idx]
        if all(rank_score .== min_rank)
            top_sites[:, :, scen] .= missing
            continue
        end

        top_sites[:, 1, scen] .= Int32.(idx)
        top_sites[:, 2, scen] .= r_ids[idx]
        top_sites[:, 3, scen] .= rank_score
    end

    return DataCube(
        top_sites;
        ranks=1:n,
        locations=[:loc_id, :loc_name, :rank],
        scenarios=1:n_scenarios
    )
end

"""
    top_N_sites(rs::ResultSet; N::Int64; metric::relative_cover)
    top_N_sites(data::AbstractArray{Real}, N::Int64; stat=mean)

Return the top `N` sites according to the provided metric (defaulting to `mean` of `relative_cover`).

# Arguments
- rs : ResultSet
- N : Number of best performing sites to be selected
- metric : Metric to use to order sites from best to worst,
           must take ResultSet as input
- stat : Summary statistic to use for comparison (default: mean)

# Returns
YAXArray[:scenarios, :locations], where `locations` indicates order of location ranking.

# Example
```julia
ADRIA.metrics.top_N_sites(rs, 5)
ADRIA.metrics.top_N_sites(rs, 5; metric=ADRIA.metric.relative_cover)
ADRIA.metrics.top_N_sites(rs, 5; metric=ADRIA.metric.relative_cover, stat=median)
```
"""
function top_N_sites(rs::ResultSet, N::Int64; metric=relative_cover, stat=mean)::YAXArray
    return top_N_sites(metric(rs), N; stat=stat)
end
function top_N_sites(data::AbstractArray{<:Real}, N::Int64; stat=mean)
    stat_m = dropdims(stat(data; dims=:timesteps); dims=:timesteps)

    top_N_sites = zeros(Int64, size(stat_m, :scenarios), N)
    for scen in axes(stat_m, :scenarios)
        # sort each scenario according to metric and get indexes
        inds = sortperm(stat_m[:, scen]; rev=true)
        top_N_sites[scen, :] = inds[1:N]
    end

    return DataCube(
        top_N_sites;
        scenarios=1:size(stat_m, :scenarios),
        locations=data.locations  # Note: assumes data holds location dimension
    )
end
