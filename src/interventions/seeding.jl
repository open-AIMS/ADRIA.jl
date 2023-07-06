"""
    distribute_seeded_corals(total_site_area::Vector{Float64},
        prefseedsites::Vector{Int64}, available_space::Vector{Float64},
        seeded_area::NamedDimsArray)::NamedDimsArray

Calculate proportion of Tabular Acropora (TA) and Corymbose Acropora (CA) to
be seeded at each of the n_site_int seeding sites selected. Distributes seeded
corals according to current available space at each selected site.

# Arguments
- total_site_area : nsites*1, total area at each site in m²  .
- prefseedsites : n_site_int*1, indices for the selected seeding sites.
- available_space : nsites*1, current available space at each site in m².
- seeded_area : area (in m²) of each coral type to be seeded with dim taxa.

# Returns
scaled_seed : NamedDimsArray [taxa to seed ⋅ locations]
"""
function distribute_seeded_corals(total_site_area::Vector{Float64},
    prefseedsites::Vector{Int64}, available_space::Vector{Float64},
    seeded_area::NamedDimsArray)::NamedDimsArray

    # Extract site area for sites selected
    site_area_seed = total_site_area[prefseedsites]

    # Proportion of available space on each site relative to total space available on these
    # sites
    prop_area_avail = available_space[prefseedsites] ./ sum(available_space[prefseedsites])

    # Distribute seeded corals (as area) across sites according to available space 
    # proportions:
    #     proportion * (area of 1 coral * num seeded corals)
    # Convert to relative cover proportion by dividing by site area
    scaled_seed = ((prop_area_avail .* seeded_area') ./ site_area_seed)'

    return scaled_seed
end


function seed_corals!(Y_pstep, a_adapt, total_location_area, prefseedsites, leftover_space_m²,
    seeded_area, Yseed, seed_sc_TA, seed_sc_CA, c_dist_t)
    # Calculate proportion to seed based on current available space
    scaled_seed = distribute_seeded_corals(total_location_area, prefseedsites, leftover_space_m², seeded_area)

    # Seed each site with TA or CA
    @views Y_pstep[seed_sc_TA, prefseedsites] .+= scaled_seed.TA
    @views Y_pstep[seed_sc_CA, prefseedsites] .+= scaled_seed.CA

    # Log seed values/sites (these values are relative to site area)
    Yseed[1, prefseedsites] .= scaled_seed.TA
    Yseed[2, prefseedsites] .= scaled_seed.CA

    w_TA::Vector{Float64} = scaled_seed.TA ./ Y_pstep[seed_sc_TA, prefseedsites]
    w_CA::Vector{Float64} = scaled_seed.CA ./ Y_pstep[seed_sc_CA, prefseedsites]

    # Update critical DHW distribution for deployed size classes
    @floop for (i, loc) in enumerate(prefseedsites)
        # Previous distributions
        TA_c_dist_t = c_dist_t[seed_sc_TA, loc]
        CA_c_dist_t = c_dist_t[seed_sc_CA, loc]

        # Priors
        wta = Float64[w_TA[i], 1.0-w_TA[i]]
        wca = Float64[w_CA[i], 1.0-w_CA[i]]

        # TN distributions for deployed corals
        TN_TA = TruncatedNormal(a_adapt[seed_sc_TA], std.(TA_c_dist_t), 0.0, maximum(TA_c_dist_t))
        TN_CA = TruncatedNormal(a_adapt[seed_sc_CA], std.(CA_c_dist_t), 0.0, maximum(CA_c_dist_t))

        c_dist_t[seed_sc_TA, loc] = MixtureModel([TA_c_dist_t, TN_TA], wta)
        c_dist_t[seed_sc_CA, loc] = MixtureModel([CA_c_dist_t, TN_CA], wca)
    end
end
