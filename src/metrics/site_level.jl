"""Functions and methods to produce location-level summaries."""

"""
    per_loc(metric, data::NamedDimsArray{D,T,N,A})::NamedDimsArray where {D,T,N,A}

Alias for summarize(data, [:scenarios, :timesteps], metric). Get metric results applied to
the location-level at indicated time (or across timesteps).

# Arguments
- metric : Any function (nominally from the Statistics package) to be applied to `data`
- data : Data set to apply metric to
- timesteps : timesteps to apply `metric` across

# Returns
Named Vector of \$N\$ elements, where \$N\$ is the number of sites.
"""
function per_loc(metric, data::NamedDimsArray{D,T,N,A})::NamedDimsArray where {D,T,N,A}
    return summarize(data, [:scenarios, :timesteps], metric)
end
function per_loc(
    metric, data::NamedDimsArray{D,T,N,A}, timesteps::Union{UnitRange,Int64}
)::NamedDimsArray where {D,T,N,A}
    return summarize(data[timesteps=timesteps], [:scenarios, :timesteps], metric)
end

"""
    loc_trajectory(metric, data::NamedDimsArray{D,T,N,A})::NamedDimsArray where {D,T,N,A}

Alias for summarize(data, [:scenarios], metric). Collate trajectory for each location.

# Examples
```julia
using Statistics

rs = ADRIA.load_results("some results")
tac = ADRIA.metrics.total_absolute_cover(rs)

# Get median trajectory for each site
ADRIA.metrics.loc_trajectory(median, tac)
# 2-dimensional NamedDimsArray(KeyedArray(...)) with keys:
# ↓   timesteps ∈ 75-element Vector{Any}
# →   sites ∈ 216-element Vector{Any}
# And data, 75×216 Matrix{Float32}:

# Get upper 95% CI for each site
ADRIA.metrics.loc_trajectory(x -> quantile(x, 0.975), tac)
# 2-dimensional NamedDimsArray(KeyedArray(...)) with keys:
# ↓   timesteps ∈ 75-element Vector{Any}
# →   sites ∈ 216-element Vector{Any}
# And data, 75×216 Matrix{Float32}:
```

# Arguments
- metric : Any function (nominally from the Statistics package) to be applied to `data`
- data : Data set to apply metric to

# Returns
2D array of \$T ⋅ S\$, where \$T\$ is total number of time steps and \$S\$ is number of sites
"""
function loc_trajectory(
    metric, data::NamedDimsArray{D,T,N,A}
)::NamedDimsArray where {D,T,N,A}
    return summarize(data, [:scenarios], metric)
end
function loc_trajectory(
    metric, data::NamedDimsArray{D,T,N,A}, timesteps::Union{UnitRange,Int64}
)::NamedDimsArray where {D,T,N,A}
    return summarize(data[timesteps=timesteps], [:scenarios], metric)
end

"""
    summarize(data::NamedDimsArray{<:Real}, alongs_axis::Vector{Symbol}, metric::Function)::NamedDimsArray{<:Real}
    summarize(data::NamedDimsArray{<:Real}, alongs_axis::Vector{Symbol}, metric::Function, timesteps::Union{UnitRange,Vector{Int64},BitVector})::NamedDimsArray{<:Real}

Apply summary metric along some axis of a data set across some or all timesteps.

# Arguments
- `data` : Data set to apply metric to.
- `alongs_axis` : which axis will be replaced with (:) when slicing.
- `metric` : Any function (nominally from the Statistics package) to be applied to `data`.
- `timesteps` : timesteps to apply `metric` across.

# Returns
NamedDimsArray with summary metric for the remaining axis.
"""
function summarize(
    data::NamedDimsArray{D,T,N,A}, alongs_axis::Vector{Symbol}, metric::Function
)::NamedDimsArray where {D,T,N,A}
    # Get length of timestep dimension directly
    #   `map` erroneously extracts every single element from the NamedDimsArray
    #   so we use `timesteps` to subset the dataset.
    timesteps = axes(data, :timesteps)
    _data = data[timesteps=timesteps]  # Note: use axes instead of axiskeys for speed
    alongs = sort([NamedDims.dim(data, axis) for axis in alongs_axis])

    keep_axis = 1:ndims(_data) .∉ [alongs]
    dims = dimnames(_data)[keep_axis]
    axis = axiskeys(_data)[keep_axis]

    # Use of JuliennedArrays is in an attempt to speed up calculation of summary statistics.
    #   We see a small but still worthwhile improvement in practice.
    #   see: https://stackoverflow.com/a/62040897
    data_slices = JuliennedArrays.Slices(_data, alongs...)
    summarized_data = map(metric, data_slices)
    summarized_keys = NamedTuple{Tuple(dims)}(axis)

    return NamedDimsArray(summarized_data; summarized_keys...)
end
function summarize(
    data::NamedDimsArray{D,T,N,A},
    alongs_axis::Vector{Symbol},
    metric::Function,
    timesteps::Union{UnitRange,Vector{Int64},BitVector},
)::NamedDimsArray where {D,T,N,A}
    return summarize(data[timesteps=timesteps], alongs_axis, metric)
end
