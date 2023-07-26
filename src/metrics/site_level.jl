"""Functions and methods to produce location-level summaries."""


"""
    per_loc(metric, data::NamedDimsArray{D,T,N,A})::NamedDimsArray where {D,T,N,A}

Get metric results applied to the location-level at indicated time (or across timesteps).

# Arguments
- metric : Any function (nominally from the Statistics package) to be applied to `data`
- data : Data set to apply metric to
- timesteps : timesteps to apply `metric` across

# Returns
Named Vector of \$N\$ elements, where \$N\$ is the number of sites.
"""
function per_loc(metric, data::NamedDimsArray{D,T,N,A})::NamedDimsArray where {D,T,N,A}
    # Get length of timestep dimension directly
    #   `map` erroneously extracts every single element from the NamedDimsArray
    #   so we use `tf` to subset the dataset.
    # Use of JuliennedArrays is in an attempt to speed up calculation of summary statistics.
    #   We see a small but still worthwhile improvement in practice.
    #   see: https://stackoverflow.com/a/62040897
    tf = axes(data, :timesteps)
    s::Vector{eltype(data)} = map(metric,
        JuliennedArrays.Slices(data[timesteps=tf], dim(data, :timesteps), dim(data, :scenarios))
    )

    return NamedDimsArray(s, sites=axiskeys(data, :sites))
end
function per_loc(metric, data::NamedDimsArray{D,T,N,A}, timesteps::Union{UnitRange,Int64})::NamedDimsArray where {D,T,N,A}
    return per_loc(metric, data[timesteps=timesteps])
end

"""
    loc_trajectory(metric, data::NamedDimsArray{D,T,N,A})::NamedDimsArray where {D,T,N,A}

Collate trajectory for each location.

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
2D array of of \$T ⋅ S\$, where \$T\$ is total number of time steps and \$S\$ is number of sites
"""
function loc_trajectory(metric, data::NamedDimsArray{D,T,N,A})::NamedDimsArray where {D,T,N,A}
    tf = axes(data, :timesteps)  # Note: use axes instead of axiskeys for speed
    s::Matrix{eltype(data)} = map(metric,
        JuliennedArrays.Slices(data[timesteps=tf], dim(data, :scenarios))
    )

    return NamedDimsArray(s, timesteps=axiskeys(data, :timesteps), sites=axiskeys(data, :sites))
end
function loc_trajectory(metric, data::NamedDimsArray{D,T,N,A}, timesteps::Union{UnitRange,Int64})::NamedDimsArray where {D,T,N,A}
    return loc_trajectory(metric, data[timesteps=timesteps])
end
