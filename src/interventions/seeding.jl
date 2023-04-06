"""
    distribute_seeded_corals(total_location_area::Vector{Float64},
        prefseedlocations::Vector{Int64}, available_space::Vector{Float64},
        seeded_area::NamedTuple{(:TA, :CA),Tuple{Float64,Float64}})::NamedTuple{(:TA, :CA),Tuple{Vector{Float64},Vector{Float64}}}

Calculate proportion of Tabular Acropora (TA) and Corymbose Acropora (CA) to
be seeded at each of the n_location_int seeding locations selected. Distributes seeded
corals according to current available space at each selected location.

# Arguments
- total_location_area : nlocations*1, total area at each location in m²  .
- prefseedlocations : n_location_int*1, indices for the selected seeding locations.
- available_space : nlocations*1, current available space at each location in m².
- seeded_area : area (in m²) of each coral type to be seeded with
    - TA : colony area of a TA.
    - CA : colony area of a CA.

# Returns
scaled_seed : NamedTuple (TA = scaled_seed_TA, CA = scaled_seed_CA), where:
    - TA : n_location_int elements, proportions of TA coral to seed at prefseedlocations
    - CA : n_location_int elements, proportions of CA coral to seed at prefseedlocations.
"""
function distribute_seeded_corals(total_location_area::Vector{Float64},
    prefseedlocations::Vector{Int64}, available_space::Vector{Float64},
    seeded_area::NamedTuple{(:TA, :CA),Tuple{Float64,Float64}})::NamedTuple{(:TA, :CA),Tuple{Vector{Float64},Vector{Float64}}}

    # extract location area for locations selected
    location_area_seed = total_location_area[prefseedlocations]

    # proportion of available space on each location relative to total space available on these locations
    # prop_area_avail = available_space ./ sum(available_space)
    prop_area_avail = available_space[prefseedlocations] ./ sum(available_space[prefseedlocations])

    # distribute seeded corals (as area) across locations according to available space proportions
    # proportion*(area of 1 coral * num seeded corals)
    # convert to relative cover proportion by dividing by location area

    scaled_seed_TA = (prop_area_avail .* seeded_area.TA) ./ location_area_seed
    scaled_seed_CA = (prop_area_avail .* seeded_area.CA) ./ location_area_seed

    return (TA=scaled_seed_TA, CA=scaled_seed_CA)
end
