"""
    CoralGrowth(n_locs::Integer, n_groups::Integer, n_sizes::Integer, ode_p::NamedTuple)

Coral growth specification for growth ODE model.
"""
struct CoralGrowth{A<:Integer,T<:NamedTuple}
    n_locs::A
    n_groups::A
    n_sizes::A
    n_group_and_size::A

    ode_p::T
end

"""
    CoralGrowth(n_locs)

Implements temporary hardcoded caches for a scenario with 35 'species' (split into 5 groups).
"""
function CoralGrowth(n_locs::Int64)
    n_groups::Int64, n_sizes::Int64 = size(linear_extensions())
    return CoralGrowth(n_locs, n_groups, n_sizes)
end
function CoralGrowth(n_locs::Int64, n_groups::Int64, n_sizes::Int64)::CoralGrowth
    n_group_and_size::Int64 = n_groups * n_sizes

    # Store specific indices for use in growth ODE function
    # The following splits into "small, mid and large" indices based on number of size classes and functional groups
    indices_array = zeros(Int64, (n_groups, n_sizes))
    for n_g in 0:(n_groups - 1)
        indices_array[n_g + 1, :] .= collect((n_g * n_sizes + 1):(n_g * n_sizes + n_sizes))
    end

    mid_ind = collect(indices_array[:, 2:(end - 1)])[:]
    size_mid = prod(size(indices_array[:, 2:(end - 1)]))

    small::SVector{n_groups} = indices_array[:, 1]
    mid::SVector{size_mid} = mid_ind
    large::SVector{n_groups} = indices_array[:, end]
    p = @NamedTuple{
        small::StaticArrays.SVector{n_groups,Int64},           # indices for small size classes
        mid::StaticArrays.SVector{size_mid,Int64},            # indices for mid-size corals
        large::StaticArrays.SVector{n_groups,Int64},           # indices for large corals
        rec::Matrix{Float64},                           # recruitment values, where `s` relates to available space (not max carrying capacity)
        sXr::Array{Float64},                           # s * X * r
        X_mb::Array{Float64},                          # X * mb
        r::Matrix{Float64},                             # growth rate
        mb::Matrix{Float64}                             # background mortality
    }((                        # cache matrix to hold X (current coral cover)
        # cached indices
        small, mid, large,

        # cache matrices
        zeros(n_groups, n_locs),  # rec
        zeros(n_groups, n_sizes, n_locs), # sXr
        zeros(n_groups, n_sizes, n_locs),  # X_mb
        zeros(n_groups, n_sizes),  # r
        zeros(n_groups, n_sizes)   # mb
    ))

    return CoralGrowth(n_locs, n_groups, n_sizes, n_group_and_size, p)
end
