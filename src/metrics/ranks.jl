"""
    _get_ranks(rs::ResultSet, intervention::Int64; kwargs...)

Extracts results for a specific intervention (seeding [1] or shading [2])
"""
function _get_ranks(rs::ResultSet, intervention::Int64; kwargs...)
    return slice_results(rs.ranks[intervention=intervention]; kwargs...)
end

"""
    _collate_ranks(rs, selected)

Collates ranks into seed/shade ranking results into a common structure.
"""
function _collate_ranks(rs::ResultSet, selected; kwargs...)::NamedDimsArray
    n_steps, n_locations = size(selected)

    ts = timesteps(rs)
    @assert length(ts) == n_steps

    r_ids = rs.location_data.reef_siteid
    if haskey(kwargs, :locations)
        r_ids = r_ids[kwargs[:locations]]
    end

    if length(r_ids) != n_locations
        @warn "Length of reef ids do not match number of locations"
    end

    return NamedDimsArray(selected; zip((:timesteps, :locations, :scenarios), (ts, r_ids, 1:size(selected, 3)))...)
end


"""
    seed_ranks(rs::ResultSet; kwargs...)

# Arguments
- rs : ResultSet
- kwargs : named dimensions to slice across

# Returns
NamedDimsArray[timesteps, locations, scenarios]

# Example
```julia
ADRIA.metrics.seed_ranks(rs; timesteps=1:10, scenarios=3:5)
```
"""
function seed_ranks(rs::ResultSet; kwargs...)
    selected = _get_ranks(rs, 1; kwargs...)
    return _collate_ranks(rs, selected; kwargs...)
end

"""
    top_n_seeded_locations(rs::ResultSet, n::Int64; kwargs...)

Get the top n seeded locations over time by their unique location id.
Lower rank values are better (e.g., 1 = first choice)

# Arguments
- rs : ResultSet
- n : `n` locations to retrieve
- kwargs : dimensions to slice across

# Returns
NamedDimsArray[locations, [Site Index, Unique ID, Rank], scenarios]
"""
function top_n_seeded_locations(rs::ResultSet, n::Int64; kwargs...)
    ranked_locations = seed_ranks(rs; kwargs...)

    r_ids = rs.location_data.reef_siteid
    min_rank = length(r_ids) + 1

    c_ranks = mean(ranked_locations, dims=1)
    top_locations = Array{Union{String,Int32,Float32,Missing}}(undef, n, 3, size(ranked_locations, 3))
    for scen in axes(ranked_locations, 3)
        flat = vec(c_ranks[1, :, scen])

        idx = partialsortperm(flat, 1:n)

        rank_score = flat[idx]
        if all(rank_score .== min_rank)
            top_locations[:, :, scen] .= missing
            continue
        end

        top_locations[:, 1, scen] .= Int32.(idx)
        top_locations[:, 2, scen] .= r_ids[idx]
        top_locations[:, 3, scen] .= rank_score
    end

    return NamedDimsArray(top_locations, (:locations, :location_ranks, :scenarios))
end

"""
    shade_ranks(rs::ResultSet; kwargs...)

# Arguments
- rs : ResultSet
- kwargs : named dimensions to slice across

# Returns
NamedDimsArray[timesteps, locations, scenarios]

# Example
```julia
ADRIA.metrics.shade_ranks(rs; timesteps=1:10, scenarios=3:5)
```
"""
function shade_ranks(rs::ResultSet; kwargs...)
    selected = _get_ranks(rs, 2; kwargs...)
    return _collate_ranks(rs, selected; kwargs...)
end

"""
    top_N_locations(rs::ResultSet; N::Int64; metric::relative_cover)
    top_N_locations(data::AbstractArray{Real}, N::Int64; stat=mean)

Return the top `N` locations according to the provided metric (defaulting to `mean` of `relative_cover`).


# Arguments
- rs : ResultSet
- N : No. of best performing locations to be selected
- metric : metric to use to order locations from best to worst,
           must take ResultSet as input
- stat : summary statistic to use for comparison (default: mean)

# Returns
NamedDimsArray[:scenarios,:location_ids], where `location_ids` indicates order of location ranking as well.

# Example
```julia
ADRIA.metrics.top_N_locations(rs, 5)
ADRIA.metrics.top_N_locations(rs, 5; metric=ADRIA.metric.relative_cover)
ADRIA.metrics.top_N_locations(rs, 5; metric=ADRIA.metric.relative_cover, stat=median)
```
"""
function top_N_locations(rs::ResultSet, N::Int64; metric=relative_cover, stat=mean)::NamedDimsArray
    return top_N_locations(metric(rs), N; stat=stat)
end
function top_N_locations(data::AbstractArray{<:Real}, N::Int64; stat=mean)
    stat_m = dropdims(stat(data, dims=:timesteps), dims=:timesteps)

    top_N_locations = zeros(Int64, size(stat_m, :scenarios), N)
    for scen in axes(stat_m, :scenarios)
        # sort each scenario according to metric and get indexes
        inds = sortperm(stat_m[:, scen], rev=true)
        top_N_locations[scen, :] = inds[1:N]
    end

    return NamedDimsArray(top_N_locations, (:scenarios, :location_ids))
end
