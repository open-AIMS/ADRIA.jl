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
