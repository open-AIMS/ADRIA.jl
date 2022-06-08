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
    small_r::SVector = @SVector [1, 2, 4, 6]
    small::SVector = @SVector [1, 7, 19, 31]
    mid::SVector = SVector{19}(collect([2:4; 8:10; 14:17; 20:23; 29; 32:35]))
    large::SVector = @SVector [18, 24, 30, 36]

    enc::SVector = @SVector [3, 5]
    encrusting::SVector = @SVector [13, 25]
    sel_en::SVector = @SVector [5, 11]
    sel_unen::SVector = @SVector [6, 12]

    p = @NamedTuple{r::Matrix{Float64}, P::Float64, mb::Matrix{Float64}, comp::Float64, small_massives::StaticArrays.SVector{3, Int64},
                     small::StaticArrays.SVector{4, Int64}, mid::StaticArrays.SVector{19, Int64}, large::StaticArrays.SVector{4, Int64},
                     sel_en::StaticArrays.SVector{2, Int64}, sel_unen::StaticArrays.SVector{2, Int64}, encrusting::StaticArrays.SVector{2, Int64},
                     small_r::StaticArrays.SVector{4, Int64}, enc::StaticArrays.SVector{2, Int64}, rec::Matrix{Float64}, k::Matrix{Float64},
                     kX_sel_en::Matrix{Float64}, X_tab::Matrix{Float64}, kXr::Matrix{Float64}, k_rec::Matrix{Float64},
                     X_mb::Matrix{Float64}}((
        zeros(n_species, 1), 0.3, zeros(n_species, 1), 0.3,

        # cached indices
        small_massives, small, mid, large,
        sel_en, sel_unen, encrusting, small_r, enc,

        # cache matrices
        zeros(n_groups, n_sites), zeros(1, n_sites), zeros(2, n_sites), zeros(1, n_sites),
        zeros(n_species, n_sites), zeros(n_groups, n_sites), zeros(n_species, n_sites)
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
