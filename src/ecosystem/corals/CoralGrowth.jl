"""
    CoralGrowth(n_sites::Integer, n_species::Integer, n_groups::Integer, ode_p::NamedTuple)

Coral growth specification for growth ODE model.
"""
struct CoralGrowth{A<:Integer,T<:NamedTuple}
    n_sites::A
    n_species::A
    n_groups::A

    ode_p::T
end


"""
    CoralGrowth(n_sites)

Implements temporary hardcoded caches for a scenario with 35 'species' (split into 5 groups).
"""
function CoralGrowth(n_sites::Int64)::CoralGrowth
    n_species, n_groups = 35, 5

    # Store specific indices for use in growth ODE function
    # These are specific to the 36 "species"/ 6 group formulation
    small::SVector = @SVector [1, 8, 15, 22, 29]
    mid::SVector = SVector{25}(collect([2:6; 9:13; 16:20; 23:27; 30:34]))
    large::SVector = @SVector [7, 14, 21, 28, 35]

    p = @NamedTuple{
            small::StaticArrays.SVector{5,Int64},           # indices for small size classes
            mid::StaticArrays.SVector{25,Int64},            # indices for mid-size corals
            large::StaticArrays.SVector{5,Int64},           # indices for large corals
            rec::Matrix{Float64},                           # recruitment values, where `s` relates to available space (not max carrying capacity)
            sXr::Matrix{Float64},                           # s * X * r
            X_mb::Matrix{Float64},                          # X * mb
            r::Vector{Float64},                             # growth rate
            mb::Vector{Float64}                             # background mortality
        }((                        # cache matrix to hold X (current coral cover)
        # cached indices
        small, mid, large,

        # cache matrices
        zeros(n_groups, n_sites),  # rec
        zeros(n_species, n_sites), # sXr
        zeros(n_species, n_sites),  # X_mb
        zeros(n_species),  # r
        zeros(n_species)   # mb
    ))

    return CoralGrowth(n_sites, n_species, n_groups, p)
end
