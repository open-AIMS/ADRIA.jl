"""
    distribute_seeded_corals(total_site_area::Vector{Float64},
        prefseedsites::Vector{Int64}, available_space::Vector{Float64},
        seeded_area::NamedTuple{(:TA, :CA),Tuple{Float64,Float64}})::NamedTuple{(:TA, :CA),Tuple{Vector{Float64},Vector{Float64}}}

Calculate proportion of Tabular Acropora (TA) and Corymbose Acropora (CA) to
be seeded at each of the nsiteint seeding sites selected. Distributes seeded
corals according to current available space at each selected site.

# Arguments
- total_site_area : nsites*1, total area at each site in m^2.
- prefseedsites : nsiteint*1, indices for the selected seeding sites.
- available_space : nsites*1, current available space at each site in m^2.
- seeded_area : area (in m²) of each coral type to be seeded with
    - TA : colony area of a TA.
    - CA : colony area of a CA.

# Returns
scaled_seed : NamedTuple (TA = scaled_seed_TA, CA = scaled_seed_CA), where:
    - TA : nsiteint elements, proportions of TA coral to seed at prefseedsites
    - CA : nsiteint elements, proportions of CA coral to seed at prefseedsites.
"""
function distribute_seeded_corals(total_site_area::Vector{Float64},
    prefseedsites::Vector{Int64}, available_space::Vector{Float64},
    seeded_area::NamedTuple{(:TA, :CA),Tuple{Float64,Float64}})::NamedTuple{(:TA, :CA),Tuple{Vector{Float64},Vector{Float64}}}
 
    # extract site area for sites selected
    site_area_seed = total_site_area[prefseedsites]

    # proportion of available space on each site relative to total space available on these sites
    # prop_area_avail = available_space ./ sum(available_space)
    prop_area_avail = available_space[prefseedsites] ./ sum(available_space[prefseedsites])

    # distribute seeded corals (as area) across sites according to available space proportions
    # proportion*(area of 1 coral * num seeded corals)
    # convert to relative cover proportion by dividing by site area

    scaled_seed_TA = (prop_area_avail .* seeded_area.TA) ./ site_area_seed
    scaled_seed_CA = (prop_area_avail .* seeded_area.CA) ./ site_area_seed

    return (TA=scaled_seed_TA, CA=scaled_seed_CA)
end
