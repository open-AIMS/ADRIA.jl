using StatsBase

"""
    distribute_seeded_corals(seed_loc_area::Vector{Float64}, available_space::Vector{Float64}, seeded_area::YAXArray)::YAXArray

Calculate proportion of deployed corals to be seeded at each of the selected locations.
Distributes seeded corals according to current available space at each selected site.

# Arguments
- seed_loc_k_m² : Carrying capacity area of locations to seed in m².
- available_space : Currently available space at each seed location in m².
- seeded_area : Area (in m²) of each coral type to be seeded with dim taxa.
- seed_volume : Absolute number of coral to deploy.

# Returns
- YAXArray[taxa to seed ⋅ number of seed locations], Proportional increase in cover relative to locations' `k` area
- Matrix[seed locations ⋅ taxa to seed], Number of coral deployed
"""
function distribute_seeded_corals(
    seed_loc_k_m²::Vector{Float64},
    available_space::Vector{Float64},
    seeded_area::YAXArray,
    seed_volume::Vector{Float64}
)::Tuple{YAXArray,Matrix{Float64}}
    total_seeded_area::Float64 = sum(seeded_area)
    total_available_space::Float64 = sum(available_space)

    # Proportion of available space on each site relative to available space at these
    # locations
    prop_area_avail = available_space ./ total_available_space
    if total_seeded_area > total_available_space
        @warn "Seeded area exceeds available space. Restricting to available space."
        seeded_area = copy(seeded_area)
        seeded_area .*= total_available_space / total_seeded_area
    end

    # Distribute seeded corals (as area) across locations according to available space
    # proportions:
    #     proportion * (area of 1 coral * num seeded corals)
    # Convert to relative cover proportion by dividing by location area
    scaled_seed = ((prop_area_avail .* seeded_area.data') ./ seed_loc_k_m²)'
    proportional_increase = DataCube(
        scaled_seed;
        taxa=caxes(seeded_area)[1].val.data,
        locations=1:length(available_space)
    )

    n_deployed_coral =
        prop_area_avail .* seed_volume' .* max(
            total_available_space / total_seeded_area, 1.0
        )

    return proportional_increase, n_deployed_coral
end

"""
    update_tolerance_distribution!(
        scaled_seed::YAXArray,
        cover::Matrix{T},
        c_dist_t::Matrix{T},
        stdev::Vector{T},
        seed_locs::Vector{Int64},
        seed_sc::BitVector,
        a_adapt::Vector{T}
    )::Nothing where {T<:Float64}

Update the thermal tolerance distribution of a population due to seeding.
Updates the store `c_dist_t` in place.

# Arguments
- `scaled_seed` : Seeding values transformed to proportion cover increase relative to k area.
- `cover` : Current coral cover state.
- `c_dist_t` : Critical DHW distributions of corals to update (i.e., for time \$t\$).
- `stdev` : Standard deviation of DHW tolerance distributions for each functional type.
- `seed_locs` : Seeding locations
- `seed_sc` : Size classes to seed
- `a_adapt` : Level of assisted adaptation (as DHW tolerance increase)
"""
function update_tolerance_distribution!(
    scaled_seed::YAXArray,
    cover::AbstractArray{T},
    c_dist_t::AbstractArray{T},
    stdev::AbstractArray{T},
    seed_locs::Vector{Int64},
    seed_sc::AbstractMatrix{Bool},
    a_adapt::AbstractMatrix{T}
)::Nothing where {T<:Float64}
    # Calculate distribution weights using proportion of area (used as priors for MixtureModel)
    # Note: It is entirely possible for a location to be ranked in the top N, but
    #       with no deployments (for a given species). A location with 0 cover
    #       and no deployments will therefore be NaN due to zero division.
    #       These are replaced with 1.0 so that the distribution for unseeded
    #       corals are used.
    w_taxa::Matrix{Float64} = scaled_seed ./ (cover[seed_sc, seed_locs] .+ scaled_seed)
    replace!(w_taxa, NaN => 1.0)

    # Update critical DHW distribution for deployed size classes
    for (i, loc) in enumerate(seed_locs)
        # Previous distributions
        c_dist_ti = @view(c_dist_t[seed_sc, loc])

        # Truncated normal distributions for deployed corals
        # Assume same stdev and bounds as original
        tn::Vector{Float64} = truncated_normal_mean.(
            a_adapt[seed_sc], stdev[seed_sc], 0.0, a_adapt[seed_sc] .+ HEAT_UB
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
