"""
    distribute_seeded_corals(total_site_area::Vector{Float64},
        prefseedsites::Vector{Int64},available_space::Vector{Float64},
        n_to_seed::Vector{Int64},col_area_seed::Vector{Float64})

Calculate proportion of TA and CA corals to be seeded at each of the
nsiteint seeding sites selected. Distributes seeded corals according to 
current available space at each selected site.

# Arguments
- total_site_area : nsites*1, total area at each site in m^2.
- prefseedsites : nsiteint*1, indices for the selected seeding sites.
- available_space : nsites*1, current available space at each site in m^2.
- n_to_seed : 1*2, number of each coral type to be seeded with:
    - n_to_seed.nTA = no. of TA corals to be seeded.
    - n_to_seed.nCA = no. of CA corals to be seeded.
- col_area_seed : 1*2, area of each coral type to be seeded with:
    - col_area_seed.areaTA = colony area of a TA coral.
    - col_area_seed.areaCA = colony area of a CA coral.

# Returns
- Named tuple (seedTAprop = scaled_seed_TA, seedCAprop = scaled_seed_CA), where:
    - seedTAprop = nsiteint*1 vector of proportions of TA coral to seed at prefseedsites
    - seedCAprop = nsiteint*1 vector of proportions of CA coral to seed at prefseedsites.
"""
function distribute_seeded_corals(total_site_area::Vector{Float64},
    prefseedsites::Vector{Int64}, available_space::Vector{Float64},
    n_to_seed::NamedTuple{(:nTA, :nCA), Tuple{Int64, Int64}}, 
    col_area_seed::NamedTuple{(:areaTA, :areaCA), Tuple{Float64, Float64}})

    # extract site area for sites selected
    site_area_seed = total_site_area[prefseedsites]

    # scale site area for sites selected by actual available space (k/100 - sum_cover)
    site_area_seed_remaining = site_area_seed .* available_space[prefseedsites]

    # proportion of available space on each site relative to total space available on these sites
    prop_area_avail = site_area_seed_remaining ./ sum(site_area_seed_remaining)

    # distribute seeded corals (as area) across sites according to available space proportions
    # proportion*(area of 1 coral * num seeded corals)
    scaled_seed_TA = prop_area_avail .* (n_to_seed.nTA * col_area_seed.areaTA)
    scaled_seed_CA = prop_area_avail .* (n_to_seed.nCA * col_area_seed.areaCA)

    # convert to relative cover proportion by dividing by site area
    scaled_seed_TA = scaled_seed_TA ./ site_area_seed
    scaled_seed_CA = scaled_seed_CA ./ site_area_seed

    return (seedTAprop = scaled_seed_TA, seedCAprop = scaled_seed_CA)
end
