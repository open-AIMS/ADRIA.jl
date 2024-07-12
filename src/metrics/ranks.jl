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
function _collate_ranks(rs::ResultSet, selected; kwargs...)::YAXArray
    n_steps, n_sites = size(selected)

    ts = timesteps(rs)
    if :timesteps in keys(kwargs)
        try
            ts = ts[kwargs[:timesteps]]
        catch err
            if err isa BoundsError
                error("Requested timesteps not in bounds")
            end

            rethrow(err)
        else
            @assert length(ts) == n_steps "Mismatch between dataset and requested time range"
        end
    end

    r_ids = rs.site_data.reef_siteid
    if haskey(kwargs, :sites)
        r_ids = r_ids[kwargs[:sites]]
    end

    if length(r_ids) != n_sites
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

Obtain the logged ranks for each time step, location, and scenario.

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
    selected = _get_ranks(rs, 1; kwargs...)
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
    selected = _get_ranks(rs, 2; kwargs...)
    return _collate_ranks(rs, selected; kwargs...)
end

"""
    _collate_ranked_locs(data::YAXArray)::Matrix{Int64}

Collates number of ranked locations.
"""
function _collate_ranked_locs(data::YAXArray)::Matrix{Int64}
    locs = zeros(Int64, size.([data], (1,3))...)

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

    r_ids = rs.site_data.reef_siteid
    min_rank = length(r_ids) + 1

    c_ranks = collect(dropdims(mean(ranked_locs, dims=1), dims=1))

    top_sites = Array{Union{String,Int32,Float32,Missing}}(undef, n, 3, size(ranked_locs, 3))
    fill!(top_sites, missing)

    for scen in axes(ranked_locs, 3)
        flat = vec(c_ranks[:, scen])
        flat[flat .== 0.0] .= min_rank

        idx = collect(partialsortperm(flat, 1:n))

        rank_score = flat[idx]

        top_sites[:, 1, scen] .= Int32.(idx)
        top_sites[:, 2, scen] .= r_ids[idx]
        top_sites[:, 3, scen] .= rank_score
    end

    return DataCube(
        top_sites;
        ranks=collect(1:n),
        locations=[:loc_id, :loc_name, :rank],
        scenarios=1:size(ranked_locs, 3)
    )
end

"""
    shade_ranks(rs::ResultSet; kwargs...)

# Arguments
- rs : ResultSet
- kwargs : named dimensions to slice across

# Returns
YAXArray[timesteps, sites, scenarios]

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
    top_N_locs(rs::ResultSet; N::Int64; metric::relative_cover)
    top_N_locs(data::AbstractArray{Real}, N::Int64; stat=mean)

Determine the top `N` locations according to the provided metric (defaulting to `mean` of `relative_cover`).
Assumes higher metric scores should be ranked higher.

# Arguments
- rs : ResultSet
- N : Number of best performing sites to be selected
- metric : Metric to use to order sites from best to worst,
           must take ResultSet as input
- stat : Summary statistic to use for comparison (default: mean)

# Returns
YAXArray[:scenarios, :rank], set of location indices where rows relate to the scenario and
columns indicate the location rank from 1 to N.

# Example
```julia
...
scens = ADRIA.sample(dom, 64)
rs = ADRIA.run_scenarios(dom, scens, "45")

ADRIA.metrics.top_N_locs(rs, 5).data
# 64×5 Matrix{Int64}:
#  291  228  267  288  270
#  211   60   24   54   43
#    ⋮
#  106  118   67  144  102

ADRIA.metrics.top_N_locs(rs, 5; metric=ADRIA.metric.relative_cover)
ADRIA.metrics.top_N_locs(rs, 5; metric=ADRIA.metric.relative_cover, stat=median)
```
"""
function top_N_locs(rs::ResultSet, N::Int64; metric=relative_cover, stat=mean)::YAXArray
    return top_N_locs(metric(rs), N; stat=stat)
end
function top_N_locs(data::AbstractArray{<:Real}, N::Int64; stat=mean)
    stat_m = dropdims(stat(data, dims=:timesteps), dims=:timesteps)

    top_locs = zeros(Int64, size(stat_m, :scenarios), N)
    for scen in axes(stat_m, :scenarios)
        # Sort each scenario according to metric and get indexes
        inds = sortperm(stat_m[:, scen], rev=true)
        top_locs[scen, :] = inds[1:N]
    end

    return DataCube(
        top_locs;
        scenarios=1:size(stat_m, :scenarios),
        rank=1:N
    )
end
