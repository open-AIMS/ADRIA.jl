# Contributing a metric

At a minimum, all metrics must define:

1. a "private" function (i.e., one that starts with an underscore: `_`) that performs the calculation with expected arguments
2. a "private" function that accepts a `ResultSet` as its first argument.
3. a "public" `Metric` function (i.e., no underscore) with some additional metadata

The `Metric` type allows metadata regarding the expected dimension names and unit of measure.
Note that the unit of measure is optional and can be left out.

Below is the implementation of the `total_absolute_cover` metric.

```julia
function _total_absolute_cover(X::AbstractArray{<:Real}, loc_area::Vector{<:Real})::AbstractArray{<:Real}
    return dropdims(sum(X, dims=(:groups, :sizes)), dims=(:groups, :sizes)) .* loc_area'
end
function _total_absolute_cover(rs::ResultSet)::AbstractArray{<:Real}
    return rs.outcomes[:total_absolute_cover]
end
total_absolute_cover = Metric(_total_absolute_cover, (:timesteps, :sites, :scenarios), "m²")

# Unit of measure is optional, in cases where the values are non-dimensional
# some_example_metric = Metric(_some_example_metric, (:timesteps, :scenarios))
```
