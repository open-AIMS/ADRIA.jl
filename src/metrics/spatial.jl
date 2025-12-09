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
    ensemble_loc_difference(outcome::YAXArray{T,3}, scens::DataFrame; agg_metric::Union{Function,AbstractFloat}=median, diff_target=:guided, conf::Float64=0.95, rng::AbstractRNG=Random.GLOBAL_RNG)::YAXArray where {T}

Mean bootstrapped difference (counterfactual - target) between some outcome aggregated for
each location.

# Arguments
- `outcome` : Metric outcome with dimensions (:timesteps, :locations, :scenarios).
- `scens` : Scenarios DataFrame.
- `agg_metric` : Metric used to aggregate scenarios when comparing between counterfactual and
target. If it is an `AbstractFloat` between 0 and 1, it uses the `bs_metric`-th quantile.
Defaults to `median`.
- `diff_target` : Target group of scenarios to compare with. Valid options are `:guided` and
`:unguided`. Defaults to `:guided`
- `conf` : Percentile used for the confidence interval. Defaults to 0.95.
- `rng` : Pseudorandom number generator.

# Example
```
# Load domain
dom = ADRIA.load_domain(path_to_domain, "<RCP>")

# Create scenarios
num_scens = 2^6
scens = ADRIA.sample(dom, num_scens)

# Run model
rs = ADRIA.run_scenarios(dom, scens, "45")

# Calculate difference to the counterfactual for given metric
_relative_cover = metrics.relative_cover(rs)

# Compute difference between guided and counterfactual using the 0.6-th quantile
gd_res = metrics.ensemble_loc_difference(r_cover, scens; agg_metric=0.6)

# Compute difference between unguided and counterfactual using the median
ug_res = metrics.ensemble_loc_difference(r_cover, scens; diff_target=:unguided)

# Plot maps of difference to the counterfactual
ADRIA.viz.map(rs, gd_res[summary=At(:agg_value)]; diverging=true)
ADRIA.viz.map(rs, ug_res[summary=At(:agg_value)]; diverging=true)
```

# Returns
Vector with bootstrapped difference (counterfactual - guided) for each location.
"""
function ensemble_loc_difference(
    outcome::YAXArray{T,3},
    scens::DataFrame;
    agg_metric::Union{Function,AbstractFloat}=median,
    diff_target=:guided,
    conf::Float64=0.95,
    rng::AbstractRNG=Random.GLOBAL_RNG
)::YAXArray where {T}
    is_quantile_metric = isa(agg_metric, AbstractFloat)
    if is_quantile_metric && !(0 <= agg_metric <= 1)
        error("When metric is a number, it must be within the interval [0,1]")
    end

    # Mean over all timesteps
    outcomes_agg = dropdims(mean(outcome; dims=:timesteps); dims=:timesteps)

    # Counterfactual, target outcomes
    cf_outcomes = outcomes_agg[scenarios=scens.guided .== -1]
    target_outcomes = if diff_target == :guided
        outcomes_agg[scenarios=scens.guided .> 0]
    elseif diff_target == :unguided
        outcomes_agg[scenarios=scens.guided .== 0]
    else
        error("Invalid diff_target value. Valid values are :guided and :unguided.")
    end

    # Find the smallest set of outcomes because they might not have the same size
    n_cf_outcomes = size(cf_outcomes, :scenarios)
    n_target_outcomes = size(target_outcomes, :scenarios)
    min_n_outcomes = min(n_cf_outcomes, n_target_outcomes)

    # Build (counterfactual-guided) and (counterfactual-unguided) result DataCubes
    _locations = axis_labels(outcome, :locations)
    results = ZeroDataCube(;
        T=T, summary=[:lower_bound, :agg_value, :upper_bound], locations=_locations
    )

    n_locs = length(_locations)
    for loc in 1:n_locs
        cf_shuf_set::Vector{Int64} = shuffle(rng, 1:n_cf_outcomes)[1:min_n_outcomes]
        target_shuf_set::Vector{Int64} = shuffle(rng, 1:n_target_outcomes)[1:min_n_outcomes]

        @views target_diff =
            collect(target_outcomes[loc, target_shuf_set]) .-
            collect(cf_outcomes[loc, cf_shuf_set])

        bootstrap_func(x) = is_quantile_metric ? quantile(x, agg_metric) : agg_metric(x)
        cf_target_bootstrap = bootstrap(bootstrap_func, target_diff, BalancedSampling(100))
        results[[2, 1, 3], loc] .= confint(cf_target_bootstrap, PercentileConfInt(conf))[1]
    end

    return results
end
