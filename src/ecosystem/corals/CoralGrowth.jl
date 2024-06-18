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
    CoralGrowth(n_sites)

Implements temporary hardcoded caches for a scenario with 35 'species' (split into 5 groups).
"""
function CoralGrowth(n_locs::Int64)::CoralGrowth
    # TODO Check this
    n_groups, n_sizes = size(ADRIA.bin_edges(), 1), size(ADRIA.bin_edges(), 2) - 1
    n_group_and_size = n_sizes * n_groups
    # n_group_and_size, n_groups = 35, 5
    # n_sizes = Int64(n_group_and_size / n_groups)

    # Store specific indices for use in growth ODE function
    small::SVector = SVector{n_groups}(collect(1:n_sizes:n_group_and_size))

    mid_size_starts = collect(2:n_sizes:n_group_and_size)
    mid_size_ends = collect(n_sizes-1:n_sizes:n_group_and_size)
    n_mid = (n_sizes - 2) * n_groups
    mid::SVector = SVector{n_mid}(vcat([s:e for (s, e) in zip(mid_size_starts, mid_size_ends)]...))

    large::SVector = SVector{n_groups}(collect(n_sizes:n_sizes:n_group_and_size))

    # Values for ode_p
    # TODO: The named tuple was for use with ODE solvers, which is no unneeded.
    # These caches should all be moved into the `CoralGrowth` struct.
    p = @NamedTuple{
        small::StaticArrays.SVector{n_groups,Int64},  # indices for small size classes
        mid::StaticArrays.SVector{n_mid,Int64},   # indices for mid-size corals
        large::StaticArrays.SVector{n_groups,Int64},  # indices for large corals
        rec::Matrix{Float64},                  # recruitment values, where `s` relates to available space (not max carrying capacity)
        sXr::Array{Float64,3},                  # s * X * r
        X_mb::Array{Float64,3},                 # X * mb
        r::Matrix{Float64},                    # growth rate
        mb::Matrix{Float64}                    # background mortality
    }((                        # cache matrix to hold X (current coral cover)
        # Cached indices
        small, mid, large,

        # Cache matrices
        zeros(n_groups, n_locs),           # rec
        zeros(n_groups, n_sizes, n_locs),  # sXr
        zeros(n_groups, n_sizes, n_locs),  # X_mb
        zeros(n_groups, n_sizes),          # r
        zeros(n_groups, n_sizes)           # mb
    ))

    return CoralGrowth(n_locs, n_groups, n_sizes, n_group_and_size, p)
end
