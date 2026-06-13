using StaticArrays

"""
    dominates(x::Vector{<:Real}, y::Vector{<:Real})::Vector

Adapted from:
https://discourse.julialang.org/t/fast-optimized-non-dominated-sorting-algorithms/86793/7

Original function name is `dominates2()`
"""
function dominates(x::AbstractVector{<:Real}, y::AbstractVector{<:Real})::Bool
    strict_inequality_found = false
    for i in eachindex(x)
        y[i] > x[i] && return false
        strict_inequality_found |= x[i] > y[i]
    end

    return strict_inequality_found
end

"""
    nds(X::AbstractArray{<:Real}, dist::Int64=0)::Vector{Vector{<:Int}}

Naive n-dimensional non-dominated sorting.

Adapted from:
https://discourse.julialang.org/t/fast-optimized-non-dominated-sorting-algorithms/86793/7

Original function name is `nds4()`

# Arguments
X : outcomes, where rows are scenarios and columns are metric results.
dist : distance from front, where 0 is on the frontier.

# Returns
Vector of Vectors with row indices for each `dist` from frontier, where 0 is on the frontier.
"""
function nds(X::AbstractArray{<:Real}, dist::Int64=0)::Vector{Vector{<:Int}}
    fronts = Vector{Int64}[]
    ind = collect(axes(X, 1))
    a = SVector{size(X, 2)}.(eachrow(X))
    count::Int64 = 0
    while !isempty(a)
        red::BitVector = [all(x -> !dominates(x, y), a) for y in a]
        push!(fronts, ind[red])

        if count == dist
            break
        end

        count += 1

        deleteat!(ind, red)
        deleteat!(a, red)
    end

    return fronts
end
