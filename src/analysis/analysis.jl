module analysis

using Statistics, DataFrames
using ADRIA
import ADRIA: ResultSet


"""
    normalize(data::Matrix)::Matrix

Normalize a matrix so that the data is ∈ [0, 1] relative to values in each column.
"""
function normalize(data::Matrix)::Matrix
    limits = extrema.(eachcol(data))

    scaled = hcat([(d .- mi) ./ (ma - mi) for (d, (mi, ma)) in zip(eachcol(data), limits)]...)
    replace!(scaled, NaN => 0.0)
    return scaled
end


include("pareto.jl")


end