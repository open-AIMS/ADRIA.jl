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

Implements temporary hardcoded caches for a scenario with 36 'species' (split into 6 groups).
"""
function CoralGrowth(n_sites::Int64)::CoralGrowth
    n_species, n_groups = 36, 6

    # Store specific indices for use in growth ODE function
    # These are specific to the 36 "species"/ 6 group formulation
    small::SVector = @SVector [1, 7, 13, 19, 25, 31]
    mid::SVector = SVector{24}(collect([2:5; 8:11; 14:17; 20:23; 26:29; 32:35]))
    large::SVector = @SVector [6, 12, 18, 24, 30, 36]

    p = @NamedTuple{
        r::Matrix{Float64},   # growth rate
        k::Vector{Float64},   # max carrying capacity
        mb::Matrix{Float64},  # background mortality
        small::StaticArrays.SVector{6,Int64},           # indices for small size classes
        mid::StaticArrays.SVector{24,Int64},            # indices for mid-size corals
        large::StaticArrays.SVector{6,Int64},           # indices for large corals
        rec::Matrix{Float64},                            # recruitment values, where `s` relates to available space (not max carrying capacity)
        sXr::Matrix{Float64},                            # s * X * r
        X_mb::Matrix{Float64},                           # X * mb
        cover::Vector{Float64}}((                        # cache matrix to hold X (current coral cover)
        # r, s, mb,
        zeros(n_species, 1), zeros(n_sites), zeros(n_species, 1),


        # cached indices
        small, mid, large,

        # cache matrices
        # rec, sXr, 
        # X_mb, cover
        zeros(n_groups, n_sites),zeros(n_species, n_sites), 
        zeros(n_species, n_sites), zeros(n_sites)
    ))

    return CoralGrowth(n_sites, n_species, n_groups, p)
end
