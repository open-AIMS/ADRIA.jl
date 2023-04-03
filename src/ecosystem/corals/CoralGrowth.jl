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
    small_massives::SVector = SVector{3}(collect(26:28))
    small::SVector = @SVector [1, 7, 13, 19, 25, 31]
    mid::SVector = SVector{19}(collect([2:4; 8:10; 14:17; 20:23; 29; 32:35]))
    large::SVector = @SVector [18, 24, 30, 36]

    acr_5_11::SVector = @SVector [5, 11]
    acr_6_12::SVector = @SVector [6, 12]

    p = @NamedTuple{
        r::Matrix{Float64},   # growth rate
        k::Vector{Float64},   # max carrying capacity
        mb::Matrix{Float64},  # background mortality
        comp::Float64,        # competition between small and large 
        sm_comp::Matrix{Float64},  # tmp store for competition between tab and small massives
        small_massives::StaticArrays.SVector{3,Int64},  # index locations for small massives
        small::StaticArrays.SVector{6,Int64},           # indices for small size classes
        mid::StaticArrays.SVector{19,Int64},            # indices for mid-size corals
        large::StaticArrays.SVector{4,Int64},           # indices for large corals
        acr_5_11::StaticArrays.SVector{2,Int64},          # size 5 Tabular Acropora (enhanced and unenhanced)
        acr_6_12::StaticArrays.SVector{2,Int64},        # size 6 Tabular Acropora (enhanced and unenhanced)
        rec::Matrix{Float64},                            # recruitment values, where `s` relates to available space (not max carrying capacity)
        M_sm::Matrix{Float64},                           # mortality for small massive corals due to competition and background mortality
        sXr::Matrix{Float64},                            # s * X * r
        X_mb::Matrix{Float64},                           # X * mb
        cover::Vector{Float64}}((                        # cache matrix to hold X (current coral cover)
        # r, s, mb, comp, sm_comp
        zeros(n_species, 1), zeros(n_sites), zeros(n_species, 1), 0.3, zeros(2, n_sites),


        # cached indices
        small_massives, small, mid, large,
        acr_5_11, acr_6_12,

        # cache matrices
        # rec, M_sm, 
        # sXr, X_mb, cover
        zeros(n_groups, n_sites), zeros(3, n_sites),
        zeros(n_species, n_sites), zeros(n_species, n_sites), zeros(n_sites)
    ))

    return CoralGrowth(n_sites, n_species, n_groups, p)
end
