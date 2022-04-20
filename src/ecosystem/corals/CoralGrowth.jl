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

Implements temporary caches for a scenario with 36 'species' (split into 6 groups)
"""
function CoralGrowth(n_sites)
    n_species, n_groups = 36, 6

    # Store specific indices for use in growth ODE function
    # These are specific to the 36 "species"/ 6 group formulation
    small_massives = SVector{3}(collect(26:28))
    small_r = @SVector [1, 2, 4, 6]
    small = @SVector [1, 7, 19, 31]
    mid = SVector{19}(collect([2:4; 8:10; 14:17; 20:23; 29; 32:35]))
    large = @SVector [18, 24, 30, 36]

    enc = @SVector [3, 5]
    encrusting = @SVector [13, 25]
    sel_en = @SVector [5, 11]
    sel_unen = @SVector [6, 12]

    p = (
        # static inputs
        r=rand(n_species, 1), P=0.3, mb=rand(n_species, 1), comp=0.3,

        # cached indices
        small_massives=small_massives, small=small, mid=mid, large=large,
        sel_en=sel_en, sel_unen=sel_unen, encrusting=encrusting, small_r=small_r, enc=enc,

        # cache matrices
        rec=zeros(n_groups, n_sites), k=zeros(1, n_sites), kX_sel_en=zeros(2, n_sites), X_tab=zeros(1, n_sites),
        kXr=zeros(n_species, n_sites), k_rec=zeros(n_groups, n_sites), X_mb=zeros(n_species, n_sites)
    )

    return CoralGrowth(n_sites, n_species, n_groups, p)
end
