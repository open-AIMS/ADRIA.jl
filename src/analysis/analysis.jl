module analysis

using Statistics, DataFrames
using ADRIA
import ADRIA: ResultSet


"""
    normalize(data::AbstractArray)::AbstractArray
    normalize(data::Vector)::Vector

Normalize a matrix (∈ [0, 1]) on a per-column basis.
"""
function normalize(data::AbstractMatrix)::AbstractMatrix
    return hcat(normalize.(eachcol(data))...)
end
function normalize(data::AbstractVector)::AbstractVector
    limits = extrema(data)
    scaled = let (mi, ma) = limits
        (data .- mi) ./ (ma - mi)
    end

    replace!(scaled, NaN => 0.0)
    return scaled
end


include("pareto.jl")
include("sensitivity.jl")


end