"""
    CoralGrowth(n_sites::Integer, n_species::Integer, n_groups::Integer, ode_p::NamedTuple)

Coral growth specification for growth ODE model.
"""
struct CoralGrowth{A<:Integer, T<:NamedTuple}
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
    small_r::SVector = @SVector [1, 2, 3, 4, 5, 6]
    small::SVector = @SVector [1, 7, 13, 19, 25, 31]
    mid::SVector = SVector{19}(collect([2:4; 8:10; 14:17; 20:23; 29; 32:35]))
    large::SVector = @SVector [18, 24, 30, 36]

    acr_5::SVector = @SVector [5, 11]
    acr_6::SVector = @SVector [6, 12]

    p = @NamedTuple{
            r::Matrix{Float64},   # growth rate
            k::Vector{Float64},   # max carrying capacity
            mb::Matrix{Float64},  # background mortality
            comp::Float64,        # competition between small and large 
            r_comp::Matrix{Float64},  # tmp store for competition between tab and small massives
            small_massives::StaticArrays.SVector{3, Int64},  # index locations for small massives
            small::StaticArrays.SVector{4, Int64},           # indices for small size classes
            mid::StaticArrays.SVector{19, Int64},            # indices for mid-size corals
            large::StaticArrays.SVector{4, Int64},           # indices for large corals
            acr_5::StaticArrays.SVector{2, Int64},          # size 5 Tabular Acropora (enhanced and unenhanced)
            acr_6::StaticArrays.SVector{2, Int64},        # size 6 Tabular Acropora (enhanced and unenhanced)
            small_r::StaticArrays.SVector{4, Int64},         # growth rate for small corals
            rec::Matrix{Float64},                            # recruitment values
            sigma::Matrix{Float64},                          # available space, i.e., [max carrying cap] - [current coral cover]
            sX_acr_5::Matrix{Float64},                      # cache store for k * X_{sel_en}, where `k` relates to available space (not max carrying capacity)
            M_sm::Matrix{Float64},                           # Coral cover of tabular corals
            sXr::Matrix{Float64},                            # s * X * r
            X_mb::Matrix{Float64},                           # X * mb
            cover::Vector{Float64}}((                        # cache matrix to hold X (current coral cover)
        # r, k, mb, comp, r_comp
        zeros(n_species, 1), zeros(n_sites), zeros(n_species, 1), 0.3, zeros(2, n_sites),


        # cached indices
        small_massives, small, mid, large,
        acr_5, acr_6, small_r,

        # cache matrices
        # rec, sigma, sX_sel_en, M_sm, 
        # sXr, X_mb, cover
        zeros(n_groups, n_sites), zeros(1, n_sites), zeros(2, n_sites), zeros(3, n_sites),
        zeros(n_species, n_sites), zeros(n_species, n_sites), zeros(n_sites)
    ))

    # p = (
    #     # static inputs
    #     r=rand(n_species, 1), P=0.3, mb=rand(n_species, 1), comp=0.3,

    #     # cached indices
    #     small_massives=small_massives, small=small, mid=mid, large=large,
    #     sel_en=sel_en, sel_unen=sel_unen, encrusting=encrusting, small_r=small_r, enc=enc,

    #     # cache matrices
    #     rec=zeros(n_groups, n_sites), k=zeros(1, n_sites), kX_sel_en=zeros(2, n_sites), X_tab=zeros(1, n_sites),
    #     kXr=zeros(n_species, n_sites), k_rec=zeros(n_groups, n_sites), X_mb=zeros(n_species, n_sites)
    # )

    return CoralGrowth(n_sites, n_species, n_groups, p)
end
