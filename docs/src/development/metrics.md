# Contributing a metric

At a minimum, all metrics must define:

1. a "private" function (i.e., one that starts with an underscore: `_`) that performs the calculation with expected arguments
2. a "private" function that accepts a `ResultSet` as its first argument.
3. a "public" `Metric` function (i.e., no underscore) with some additional metadata

The `Metric` type allows metadata regarding the expected dimension names and unit of measure.
Note that the unit of measure is optional and can be left out.

The core metric calculation logic is stored in the `ADRIAIndicators.jl` package.

Below is the implementation of the `total_absolute_cover` metric.

```julia
function _total_absolute_cover(relative_cover::YAXArray{<:Real,3}, k_area::Vector{<:Real})::YAXArray
    return DataCube(
        relative_cover.data .* k_area', parentmodule(metrics).axes_names(relative_cover)
    )
end
function _total_absolute_cover(rs::ResultSet)::AbstractArray{<:Real}
    return _total_absolute_cover(relative_cover(rs), loc_k_area(rs))
end
total_absolute_cover = Metric(
    _total_absolute_cover,
    (:timesteps, :locations, :scenarios),
    (:timesteps, :locations, :scenarios),
    "Cover",
    IS_NOT_RELATIVE,
    UNIT_AREA
)

# Unit of measure is optional, in cases where the values are non-dimensional
# some_example_metric = Metric(_some_example_metric, (:timesteps, :scenarios))
```
