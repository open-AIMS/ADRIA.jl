"""
    distribute_seeded_corals(seed_loc_area::Vector{Float64},
        available_space::Vector{Float64},
        seeded_area::NamedDimsArray)::NamedDimsArray

Calculate proportion of deployed corals to be seeded at each of the selected locations. 
Distributes seeded corals according to current available space at each selected site.

# Arguments
- seed_loc_area : area of locations to seed in m².
- available_space : currently available space at each seed location in m².
- seeded_area : area (in m²) of each coral type to be seeded with dim taxa.

# Returns
scaled_seed : NamedDimsArray [taxa to seed ⋅ number of seed locations]
"""
function distribute_seeded_corals(seed_loc_area::Vector{Float64},
    available_space::Vector{Float64},
    seeded_area::NamedDimsArray)::NamedDimsArray

    # Proportion of available space on each site relative to total space available on these
    # sites
    prop_area_avail = available_space ./ sum(available_space)

    # Distribute seeded corals (as area) across sites according to available space 
    # proportions:
    #     proportion * (area of 1 coral * num seeded corals)
    # Convert to relative cover proportion by dividing by site area
    scaled_seed = ((prop_area_avail .* seeded_area') ./ seed_loc_area)'

    return scaled_seed
end


function seed_corals!(Y_pstep::Matrix{Float64}, a_adapt::Vector{Float64}, total_location_area::Vector{Float64},
    prefseedsites::Vector{Int64}, leftover_space_m²::Vector{Float64}, seeded_area::NamedDimsArray,
    Yseed::Matrix{Float64}, seed_sc::Vector{Int64}, c_dist_t)

    # Calculate proportion to seed based on current available space
    scaled_seed = distribute_seeded_corals(total_location_area[prefseedsites], leftover_space_m²[prefseedsites], seeded_area)

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
        # Assume same stdev and bounds as original
        tn = TruncatedNormal.(a_adapt[seed_sc], std.(c_dist), 0.0, maximum.(c_dist))

        # Create new distributions by mixing previous and current distributions
        c_dist_t[seed_sc, loc] = map((t, t1, w) -> MixtureModel([t, t1], [w...]), c_dist, tn, wta)
    end
end
