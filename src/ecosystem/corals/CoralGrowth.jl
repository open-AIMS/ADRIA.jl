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
    n_group_and_size = n_classes * n_groups
    # n_group_and_size, n_groups = 35, 5
    # n_sizes = Int64(n_group_and_size / n_groups)

    # Store specific indices for use in growth ODE function
    # These are specific to the 35 "species"/ 5 group formulation
    small::SVector = @SVector [1, 8, 15, 22, 29]
    mid::SVector = SVector{25}(collect([2:6; 9:13; 16:20; 23:27; 30:34]))
    large::SVector = @SVector [7, 14, 21, 28, 35]

    # Values for ode_p
    # TODO: The named tuple was for use with ODE solvers, which is no unneeded.
    # These caches should all be moved into the `CoralGrowth` struct.
    p = @NamedTuple{
            small::StaticArrays.SVector{5,Int64},  # indices for small size classes
            mid::StaticArrays.SVector{25,Int64},   # indices for mid-size corals
            large::StaticArrays.SVector{5,Int64},  # indices for large corals
            rec::Matrix{Float64},                  # recruitment values, where `s` relates to available space (not max carrying capacity)
            sXr::Array{Float64, 3},                  # s * X * r
            X_mb::Array{Float64, 3},                 # X * mb
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
