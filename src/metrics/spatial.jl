using Bootstrap
using Random

"""Functions and methods to produce location-level summaries."""

"""
    per_loc(metric, data::YAXArray{D,T,N,A})::YAXArray where {D,T,N,A}

Alias for summarize(data, [:scenarios, :timesteps], metric). Get metric results applied to
the location-level at indicated time (or across timesteps).

# Arguments
- metric : Any function (nominally from the Statistics package) to be applied to `data`
- data : Data set to apply metric to
- timesteps : timesteps to apply `metric` across

# Returns
Named Vector of \$N\$ elements, where \$N\$ is the number of sites.
"""
function per_loc(metric, data::YAXArray{D,T,N,A})::YAXArray where {D,T,N,A}
    return summarize(data, [:scenarios, :timesteps], metric)
end
function per_loc(
    metric, data::YAXArray{D,T,N,A}, timesteps::Union{UnitRange,Int64}
)::YAXArray where {D,T,N,A}
    return summarize(data[timesteps=timesteps], [:scenarios, :timesteps], metric)
end

"""
    loc_trajectory(metric, data::YAXArray{D,T,N,A})::YAXArray where {D,T,N,A}

Alias for summarize(data, [:scenarios], metric). Collate trajectory for each location.

# Examples
```julia
using Statistics

rs = ADRIA.load_results("some results")
tac = ADRIA.metrics.total_absolute_cover(rs)

# Get median trajectory for each site
ADRIA.metrics.loc_trajectory(median, tac)
#75×216 YAXArray{Float64,2} with dimensions:
#  Dim{:timesteps} Categorical{Any} Any[1, 2, …, 74, 75] Unordered,
#  Dim{:sites} Categorical{Any} Any[1, 2, …, 215, 216] Unordered
#Total size: 126.56 KB

# Get upper 95% CI for each site
ADRIA.metrics.loc_trajectory(x -> quantile(x, 0.975), tac)
#75×216 YAXArray{Float64,2} with dimensions:
#  Dim{:timesteps} Categorical{Any} Any[1, 2, …, 74, 75] Unordered,
#  Dim{:sites} Categorical{Any} Any[1, 2, …, 215, 216] Unordered
#Total size: 126.56 KB

# Arguments
- metric : Any function (nominally from the Statistics package) to be applied to `data`
- data : Data set to apply metric to

# Returns
2D array of \$T ⋅ S\$, where \$T\$ is total number of time steps and \$S\$ is number of sites
"""
function loc_trajectory(
    metric, data::YAXArray{D,T,N,A}
)::YAXArray where {D,T,N,A}
    return summarize(data, [:scenarios], metric)
end
function loc_trajectory(
    metric, data::YAXArray{D,T,N,A}, timesteps::Union{UnitRange,Int64}
)::YAXArray where {D,T,N,A}
    return summarize(data[timesteps=timesteps], [:scenarios], metric)
end

"""
    summarize(data::YAXArray{<:Real}, alongs_axis::Vector{Symbol}, metric::Function)::YAXArray{<:Real}
    summarize(data::YAXArray{<:Real}, alongs_axis::Vector{Symbol}, metric::Function, timesteps::Union{UnitRange,Vector{Int64},BitVector})::YAXArray{<:Real}

Apply summary metric along some axis of a data set across some or all timesteps.

# Arguments
- `data` : Data set to apply metric to.
- `alongs_axis` : which axis will be replaced with (:) when slicing.
- `metric` : Any function (nominally from the Statistics package) to be applied to `data`.
- `timesteps` : timesteps to apply `metric` across.

# Returns
YAXArray with summary metric for the remaining axis.
"""
function summarize(
    data::YAXArray{D,T,N,A}, alongs_axis::Vector{Symbol}, metric::Function
)::YAXArray where {D,T,N,A}
    alongs = sort([axis_index(data, axis) for axis in alongs_axis])

    # Use of JuliennedArrays is in an attempt to speed up calculation of summary statistics.
    #   We see a small but still worthwhile improvement in practice.
    #   see: https://stackoverflow.com/a/62040897
    data_slices = JuliennedArrays.Slices(data.data, alongs...)
    summarized_data = map(metric, data_slices)

    new_dims = setdiff(axes_names(data), alongs_axis)
    new_axis = [axis_labels(data, ax) for ax in new_dims]

    return DataCube(summarized_data; NamedTuple{Tuple(new_dims)}(new_axis)...)
end
function summarize(
    data::YAXArray{D,T,N,A},
    alongs_axis::Vector{Symbol},
    metric::Function,
    timesteps::Union{UnitRange,Vector{Int64},BitVector},
)::YAXArray where {D,T,N,A}
    return summarize(data[timesteps=timesteps], alongs_axis, metric)
end

"""
    cf_difference_map(rs, scens, metric)

Give the mean bootstraped difference from the counterfactual for a given metric.
"""
function cf_difference_map(rs, scens, metric)
    # Extract some metric
    outcomes = metric(rs)
    locations = rs.site_data

    # Mean over all timesteps
    outcomes_mean = dropdims(mean(outcomes, dims=:timesteps), dims=:timesteps)

    cf_outcomes = outcomes_mean[scenarios=scens.guided .== -1]
    ug_outcomes = outcomes_mean[scenarios=scens.guided .== 0]
    gd_outcomes = outcomes_mean[scenarios=scens.guided .> 0]

    n_locs = length(outcomes.sites)

    ug_result = ZeroDataCube(; T=Float64, sites=locations.reef_siteid)
    gd_result = ZeroDataCube(; T=Float64, sites=locations.reef_siteid)
    n_scens = size(cf_outcomes, :scenarios)
    for loc in 1:n_locs
        cf_shuf_set1 = shuffle(1:n_scens)
        ug_shuf_set = shuffle(1:n_scens)
        ug_diff = collect(ug_outcomes[loc, ug_shuf_set] .- cf_outcomes[loc, cf_shuf_set1])
        ug_result[loc] = bootstrap(median, ug_diff, BalancedSampling(100)).t0[1]

        cf_shuf_set2 = shuffle(1:n_scens)
        gd_shuf_set = shuffle(1:n_scens)
        gd_diff = collect(gd_outcomes[loc, gd_shuf_set] .- cf_outcomes[loc, cf_shuf_set2])
        gd_result[loc] = bootstrap(median, gd_diff, BalancedSampling(100)).t0[1]
    end

    return ug_result, gd_result
end
