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
Named Vector of \$N\$ elements, where \$N\$ is the number of locations.
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
#  Dim{:locations} Categorical{Any} Any[1, 2, …, 215, 216] Unordered
#Total size: 126.56 KB

# Get upper 95% CI for each site
ADRIA.metrics.loc_trajectory(x -> quantile(x, 0.975), tac)
#75×216 YAXArray{Float64,2} with dimensions:
#  Dim{:timesteps} Categorical{Any} Any[1, 2, …, 74, 75] Unordered,
#  Dim{:locations} Categorical{Any} Any[1, 2, …, 215, 216] Unordered
#Total size: 126.56 KB
```

# Arguments
- metric : Any function (nominally from the Statistics package) to be applied to `data`
- data : Data set to apply metric to

# Returns
2D array of \$T ⋅ S\$, where \$T\$ is total number of time steps and \$S\$ is number of locations
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
    # The approach using JuliennedArrays is faster then directly passing the YAXArray to
    # mapslices only when we `read(data)`. Since `read(data)` loads all the data from disk
    # into memory, we only want to use that approach when there is a reasonable amount of
    # available space in RAM
    data_size = Base.summarysize(data)
    free_ram = Sys.free_memory()
    proportional_usage = data_size / free_ram

    # Only use this approach when data occupies more than 70% of available space in RAM
    if proportional_usage > 0.7 * free_ram
        # `D.` is ensuring the returned YAXArray has the same type as the input `data`
        return D.(mapslices(metric, data; dims=alongs_axis))
    end

    alongs = sort([axis_index(data, axis) for axis in alongs_axis])

    # Use of JuliennedArrays is in an attempt to speed up calculation of summary statistics.
    # We see a small but still worthwhile improvement in practice.
    # see: https://stackoverflow.com/a/62040897
    data_slices = JuliennedArrays.Slices(read(data), alongs...)
    summarized_data = map(metric, data_slices)

    new_dims = setdiff(axes_names(data), alongs_axis)
    new_axis = [axis_labels(data, ax) for ax in new_dims]

    return DataCube(summarized_data; NamedTuple{Tuple(new_dims)}(new_axis)...)
end
function summarize(
    data::YAXArray{D,T,N,A},
    alongs_axis::Vector{Symbol},
    metric::Function,
    timesteps::Union{UnitRange,Vector{Int64},BitVector}
)::YAXArray where {D,T,N,A}
    return summarize(data[timesteps=timesteps], alongs_axis, metric)
end
