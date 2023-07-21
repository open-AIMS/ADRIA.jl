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

"""
    seed_corals!(cover::Matrix{Float64}, total_location_area::Vector{Float64},
        leftover_space_m²::Vector{Float64}, seed_locs::BitVector,
        seeded_area::NamedDimsArray, seed_sc::BitVector, a_adapt::Vector{Float64},
        Yseed::SubArray, c_dist_t::Matrix{Distribution})::Nothing

Deploy thermally enhanced corals to indicated locations ("seeding" or "outplanting").
Increases indicated area covered by the given coral taxa and determines the modified
distribution of critical DHW thresholds.

# Arguments
- `cover` : Area currently covered by coral
- `total_location_area` : Total area for each location
- `leftover_space_m²` : Currently available area at each location
- `seed_locs` : Selected locations to seed
- `seeded_area` : Area (in m²) to seed
- `seed_sc` : Indicates in-matrix locations of the coral size classes to seed
- `a_adapt` : Mean of thermal enhancement in terms of DHW
- `Yseed` : Log of seeded locations to update
- `c_dist_t` : Critical DHW distributions of corals to update
"""
function seed_corals!(cover::Matrix{Float64}, total_location_area::Vector{Float64},
    leftover_space_m²::Vector{Float64}, seed_locs::Vector{Int64},
    seeded_area::NamedDimsArray, seed_sc::BitVector, a_adapt::Vector{Float64},
    Yseed::SubArray, c_dist_t1::Matrix{Distribution})::Nothing

    # Calculate proportion to seed based on current available space
    scaled_seed = distribute_seeded_corals(total_location_area[seed_locs], leftover_space_m²[seed_locs], seeded_area)

    # Seed each site and log
    @views cover[seed_sc, seed_locs] .+= scaled_seed
    Yseed[:, seed_locs] .= scaled_seed

    # Calculate w_taxa using vectorized operations
    w_taxa::Matrix{Float64} = scaled_seed ./ cover[seed_sc, seed_locs]

    # Update critical DHW distribution for deployed size classes
    @floop for (i, loc) in enumerate(seed_locs)
        # Previous distributions
        c_dist_t = c_dist_t1[seed_sc, loc]

        # Priors (weights based on cover for each species)
        wta = hcat(w_taxa[:, i], 1.0 .- w_taxa[:, i])

        # Truncated normal distributions for deployed corals
        # Assume same stdev and bounds as original
        tn = truncated.(Normal.(a_adapt[seed_sc], std.(c_dist_t)), 0.0, maximum.(c_dist_t))

        # Create new distributions by mixing previous and current distributions
        c_dist_t1[seed_sc, loc] = map((t, t1, w) -> MixtureModel([t, t1], [w...]), c_dist_t, tn, eachrow(wta))
    end

    return nothing
end
