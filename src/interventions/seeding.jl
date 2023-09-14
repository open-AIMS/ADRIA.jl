"""
    seed_corals!(cover::Matrix{Float64}, total_location_area::V, leftover_space::V,
        seed_locs::Vector{Int64}, seeded_area::NamedDimsArray, seed_sc::BitVector, a_adapt::V,
        Yseed::SubArray, stdev::V, c_dist_t::Matrix)::Nothing where {V<:Vector{Float64}}

Deploy thermally enhanced corals to indicated locations ("seeding" or "outplanting").
Increases indicated area covered by the given coral taxa and determines the modified
distribution of critical DHW thresholds.

Note: Units for all areas are expected to be identical, and are assumed to be in m².

# Arguments
- `cover` : Area currently covered by coral
- `total_location_area` : Total area for each location
- `leftover_space` : Currently available area at each location
- `seed_locs` : Selected locations to seed
- `seeded_area` : Area to seed
- `seed_sc` : Indicates in-matrix locations of the coral size classes to seed
- `a_adapt` : Mean of thermal enhancement in terms of DHW
- `Yseed` : Log of seeded locations to update
- `c_dist_t` : Critical DHW distributions of corals to update (i.e., for time \$t\$)
"""
function seed_corals!(cover::Matrix{Float64}, loc_k_area::V, leftover_space::V,
    seed_locs::Vector{Int64}, seeded_area::NamedDimsArray, seed_sc::BitVector, a_adapt::V,
    Yseed::SubArray, stdev::V, c_dist_t::Matrix{<:Truncated})::Nothing where {V<:Vector{Float64}}

    # Calculate proportion to seed based on current available space

    # Proportion of available space on each site relative to total space available on these
    # sites
    prop_area_avail = leftover_space[seed_locs] ./ sum(loc_k_area[seed_locs])

    # Distribute seeded corals (as area) across sites according to available space
    # proportions:
    #     proportion * (area of 1 coral * num seeded corals)
    # Convert to relative cover proportion by dividing by site area
    scaled_seed = ((prop_area_avail .* seeded_area') ./ loc_k_area[seed_locs])'
    # scaled_seed = distribute_seeded_corals(loc_k_area[seed_locs], leftover_space[seed_locs], seeded_area)

    # Seed each site and log
    @views cover[seed_sc, seed_locs] .+= scaled_seed
    Yseed[:, seed_locs] .= scaled_seed

    # Calculate distribution weights using proportion of area (used as priors for MixtureModel)
    # Note: It is entirely possible for a location to be ranked in the top N, but
    #       with no deployments (for a given species). A location with 0 cover
    #       and no deployments will therefore be NaN due to zero division.
    #       These are replaced with 1.0 so that the distribution for unseeded
    #       corals are used.
    w_taxa::Matrix{Float64} = scaled_seed ./ cover[seed_sc, seed_locs]
    replace!(w_taxa, NaN => 1.0)

    # Update critical DHW distribution for deployed size classes
    for (i, loc) in enumerate(seed_locs)
        # Previous distributions
        c_dist_ti = @view(c_dist_t[seed_sc, loc])

        # Truncated normal distributions for deployed corals
        # Assume same stdev and bounds as original
        tn::Vector{Truncated} = truncated.(Normal.(a_adapt[seed_sc], stdev[seed_sc]), minimum.(c_dist_ti), maximum.(c_dist_ti))

        # If seeding an empty location, no need to do any further calculations
        if all(isapprox.(w_taxa[:, i], 1.0))
            c_dist_t[seed_sc, loc] .= tn
            continue
        end

        # Create new distributions by mixing previous and current distributions using
        # proportional cover as the priors/weights
        # Priors (weights based on cover for each species)
        tx::Vector{MixtureModel} = MixtureModel[MixtureModel([t, t1], Float64[w, 1.0-w]) for ((t, t1), w) in zip(zip(c_dist_ti, tn), w_taxa[:, i])]
        c_dist_t[seed_sc, loc] .= truncated.(Normal.(mean.(tx), stdev[seed_sc]), minimum.(tx), maximum.(tx))
    end

    return nothing
end
