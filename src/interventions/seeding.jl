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


function seed_corals!(Y_pstep, a_adapt, total_location_area, prefseedsites,
    leftover_space_m², seeded_area, Yseed, seed_sc, c_dist_t)

    # Calculate proportion to seed based on current available space
    scaled_seed = distribute_seeded_corals(total_location_area, prefseedsites, vec(leftover_space_m²), seeded_area)

    # Seed each site and log
    @views Y_pstep[seed_sc, prefseedsites] .+= scaled_seed
    Yseed[:, prefseedsites] .= scaled_seed

    w_taxa::Matrix{Float64} = scaled_seed ./ Y_pstep[seed_sc, prefseedsites]

    # Update critical DHW distribution for deployed size classes
    @floop for (i, loc) in enumerate(prefseedsites)
        # Previous distributions
        c_dist = c_dist_t[seed_sc, loc]

        # Priors (weights based on cover for each species)
        wta = collect(zip(w_taxa[:, i], 1.0 .- w_taxa[:, i]))

        # Truncated normal distributions for deployed corals
        # assume same stdev and bounds as original
        tn = TruncatedNormal.(a_adapt[seed_sc], std.(c_dist), 0.0, maximum.(c_dist))

        # Create new distributions by mixing previous and current distributions
        c_dist_t[seed_sc, loc] = map((t, t1, w) -> MixtureModel([t, t1], [w...]), c_dist, tn, wta)
    end
end
