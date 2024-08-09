using StatsBase

"""
    distribute_seeded_corals(seed_loc_area::Vector{Float64}, available_space::Vector{Float64}, seeded_area::YAXArray)::YAXArray

Calculate proportion of deployed corals to be seeded at each of the selected locations.
Distributes seeded corals according to current available space at each selected site.

# Arguments
- seed_loc_k_m² : carrying capacity area of locations to seed in m².
- available_space : currently available space at each seed location in m².
- seeded_area : area (in m²) of each coral type to be seeded with dim taxa.

# Returns
YAXArray[taxa to seed ⋅ number of seed locations], area increased relative to k area.
"""
function distribute_seeded_corals(
    seed_loc_k_m²::Vector{Float64},
    available_space::Vector{Float64},
    seeded_area::YAXArray
)::YAXArray

    # Proportion of available space on each site relative to available space at these
    # locations
    prop_area_avail = available_space ./ sum(available_space)

    # Distribute seeded corals (as area) across locations according to available space
    # proportions:
    #     proportion * (area of 1 coral * num seeded corals)
    # Convert to relative cover proportion by dividing by location area
    scaled_seed = ((prop_area_avail .* seeded_area.data') ./ seed_loc_k_m²)'
    #scaled_seed = ((prop_area_avail .* seeded_area') ./ seed_loc_k_m²)'

    return DataCube(scaled_seed, taxa=caxes(seeded_area)[1].val.data, locations=1:length(available_space))
end

"""
    seed_corals!(cover::Matrix{Float64}, loc_k_area::V, leftover_space_m²::V, seed_locs::Vector{Int64}, seeded_area::YAXArray, seed_sc::BitVector, a_adapt::V, Yseed::SubArray, stdev::V, c_dist_t::Matrix)::Nothing where {V<:Vector{Float64}}

Deploy thermally enhanced corals to indicated locations ("seeding" or "outplanting").
Increases indicated area covered by the given coral taxa and determines the modified
distribution of critical DHW thresholds.

Note: Units for all areas are expected to be identical, and are assumed to be in m².

# Arguments
- `cover` : Area currently covered by coral
- `loc_k_area` : Carrying capacity of location (in m²)
- `leftover_space_m²` : Currently available area at each location
- `seed_locs` : Selected locations to seed
- `seeded_area` : Area to seed
- `seed_sc` : Indicates in-matrix locations of the coral size classes to seed
- `a_adapt` : Mean of thermal enhancement in terms of DHW
- `Yseed` : Log of seeded locations to update
- `c_dist_t` : Critical DHW distributions of corals to update (i.e., for time \$t\$)
"""
function seed_corals!(
    cover::Matrix{Float64},
    loc_k_area::V,
    leftover_space_m²::V,
    seed_locs::Vector{Int64},
    seeded_area::YAXArray,
    seed_sc::BitVector,
    a_adapt::V,
    Yseed::SubArray,
    stdev::V,
    c_dist_t::Matrix{Float64},
)::Nothing where {V<:Vector{Float64}}
    # Selected locations can fill up over time so avoid locations with no space
    seed_locs = seed_locs[findall(leftover_space_m²[seed_locs] .> 0.0)]

    # Calculate proportion to seed based on current available space
    scaled_seed = distribute_seeded_corals(
        loc_k_area[seed_locs],
        leftover_space_m²[seed_locs],
        seeded_area,
    )

    # Seed each location and log
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
        tn::Vector{Float64} = truncated_normal_mean.(
            a_adapt[seed_sc], stdev[seed_sc], 0.0, a_adapt[seed_sc] .+ HEAT_UB,
        )

        # If seeding an empty location, no need to do any further calculations
        if all(isapprox.(w_taxa[:, i], 1.0))
            c_dist_t[seed_sc, loc] .= tn
            continue
        end

        # Create new distributions by mixing previous and current distributions using
        # proportional cover as the priors/weights
        # Priors (weights based on cover for each species)
        tx::Vector{Weights} = Weights.(eachcol(vcat(w_taxa[:, i]', 1.0 .- w_taxa[:, i]')))
        c_dist_t[seed_sc, loc] = sum.(eachcol(vcat(c_dist_ti', tn')), tx)
    end

    return nothing
end
function seed_corals!(
    cover::AbstractArray{Float64, 3},
    loc_k_area::Vector{T},
    leftover_space_m²::Vector{T},
    seed_locs::Vector{Int64},
    seeded_area::YAXArray,
    seed_sc::Matrix{Bool},
    a_adapt::Matrix{T},
    Yseed::SubArray,
    stdev::Matrix{T},
    c_dist_t::Array{Float64, 3},
)::Nothing where {T<:Float64}
    # Selected locations can fill up over time so avoid locations with no space
    seed_locs = seed_locs[findall(leftover_space_m²[seed_locs] .> 0.0)]
    loc_mask::BitVector = [loc in seed_locs for loc in 1:length(loc_k_area)]

    # Calculate proportion to seed based on current available space
    scaled_seed = distribute_seeded_corals(
        loc_k_area[seed_locs],
        leftover_space_m²[seed_locs],
        seeded_area,
    )

    # Seed each location and log
    cover[seed_sc, loc_mask] .+= scaled_seed
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
        tn::Vector{Float64} = truncated_normal_mean.(
            a_adapt[seed_sc], stdev[seed_sc], 0.0, a_adapt[seed_sc] .+ HEAT_UB,
        )

        # If seeding an empty location, no need to do any further calculations
        if all(isapprox.(w_taxa[:, i], 1.0))
            c_dist_t[seed_sc, loc] .= tn
            continue
        end

        # Create new distributions by mixing previous and current distributions using
        # proportional cover as the priors/weights
        # Priors (weights based on cover for each species)
        tx::Vector{Weights} = Weights.(eachcol(vcat(w_taxa[:, i]', 1.0 .- w_taxa[:, i]')))
        c_dist_t[seed_sc, loc] = sum.(eachcol(vcat(c_dist_ti', tn')), tx)
    end

    return nothing
end
