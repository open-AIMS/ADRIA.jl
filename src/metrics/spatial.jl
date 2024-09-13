"""Functions and methods to produce location-level summaries."""

using Bootstrap
using Random

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
    available_ram = Sys.free_memory() * 0.7
    data_size = sizeof(eltype(data)) * length(data)

    # Only use this approach when data occupies more than 70% of available RAM
    if data_size > available_ram
        # `D.` is ensuring the returned YAXArray has the same type as the input `data`
        return fill_metadata!(D.(mapslices(metric, data; dims=alongs_axis)), metadata(data))
    end

    alongs = sort([axis_index(data, axis) for axis in alongs_axis])

    # Use of JuliennedArrays is in an attempt to speed up calculation of summary statistics.
    # We see a small but still worthwhile improvement in practice.
    # see: https://stackoverflow.com/a/62040897
    data_slices = JuliennedArrays.Slices(read(data), alongs...)
    summarized_data = map(metric, data_slices)

    new_dims = setdiff(axes_names(data), alongs_axis)
    new_axis = [axis_labels(data, ax) for ax in new_dims]

    return fill_metadata!(
        DataCube(summarized_data; NamedTuple{Tuple(new_dims)}(new_axis)...), metadata(data)
    )
end
function summarize(
    data::YAXArray{D,T,N,A},
    alongs_axis::Vector{Symbol},
    metric::Function,
    timesteps::Union{UnitRange,Vector{Int64},BitVector}
)::YAXArray where {D,T,N,A}
    return summarize(data[timesteps=timesteps], alongs_axis, metric)
end

"""
    cf_difference_loc(rs, scens, metric)

Mean bootstrapped differences (counterfactual - guided) and (counterfactual - unguided) for
each location.

# Arguments
- `outcome` : Metric outcome with dimensions (:timesteps, :locations, :scenarios)
- `scens` : Scenarios DataFrame

# Returns
Two elements tuple with mean bootstrapped difference (counterfactual - guided) and
(counterfactual - unguided) for each location.
"""
function cf_difference_loc(
    outcome::YAXArray{T,3}, scens::DataFrame; conf::Float64=0.95
)::Tuple where {T}
    # Mean over all timesteps
    outcomes_agg = dropdims(mean(outcome; dims=:timesteps); dims=:timesteps)

    # Counterfactual, guided and unguided outcomes
    cf_outcomes = outcomes_agg[scenarios=scens.guided .== -1]
    gd_outcomes = outcomes_agg[scenarios=scens.guided .> 0]
    ug_outcomes = outcomes_agg[scenarios=scens.guided .== 0]

    # Find the smallest set of outcomes because they might not have the same size
    n_cf_outcomes = size(cf_outcomes, :scenarios)
    n_gd_outcomes = size(gd_outcomes, :scenarios)
    n_ug_outcomes = size(ug_outcomes, :scenarios)
    min_n_outcomes = min(n_cf_outcomes, n_gd_outcomes, n_ug_outcomes)

    # Build (counterfactual-guided) and (counterfactual-unguided) result DataCubes
    _locations = axis_labels(outcome, :locations)
    gd_result = ZeroDataCube(;
        T=T, value=[:lower_bound, :value, :upper_bound], locations=_locations
    )
    ug_result = ZeroDataCube(;
        T=T, value=[:lower_bound, :value, :upper_bound], locations=_locations
    )

    n_locs = length(_locations)
    for loc in 1:n_locs
        cf_shuf_set::Vector{Int64} = shuffle(1:n_cf_outcomes)[1:min_n_outcomes]
        gd_shuf_set::Vector{Int64} = shuffle(1:n_gd_outcomes)[1:min_n_outcomes]
        ug_shuf_set::Vector{Int64} = shuffle(1:n_ug_outcomes)[1:min_n_outcomes]

        # Counterfactual - Guided
        gd_diff = collect(gd_outcomes[loc, gd_shuf_set] .- cf_outcomes[loc, cf_shuf_set])
        cf_gd_bootstrap = bootstrap(median, gd_diff, BalancedSampling(100))
        gd_result[[2, 1, 3], loc] .= confint(cf_gd_bootstrap, PercentileConfInt(conf))[1]

        # Counterfactual - Unguided
        ug_diff = collect(ug_outcomes[loc, ug_shuf_set] .- cf_outcomes[loc, cf_shuf_set])
        cf_ug_bootstrap = bootstrap(median, ug_diff, BalancedSampling(100))
        ug_result[[2, 1, 3], loc] .= confint(cf_ug_bootstrap, PercentileConfInt(conf))[1]
    end

    return gd_result, ug_result
end
