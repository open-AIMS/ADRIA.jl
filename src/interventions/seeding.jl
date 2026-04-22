using StatsBase

"""
    seed_size_groups(n_groups::Int64, n_sizes::Int64)::BitMatrix

BitMatrix with dimensions (n_groups, n_sizes) where trues represent size classes and
functional groups to seed. Seed only the smallest size class for each group.
"""
function seed_size_groups(n_groups::Int64, n_sizes::Int64)::BitMatrix
    return hcat(trues(n_groups), falses(n_groups, n_sizes - 1))
end

"""
    distribute_seeded_corals(
        seed_loc_k_m²::Vector{Float64},
        available_space::Vector{Float64},
        seed_volume::Vector{Float64},
        colony_areas::Vector{Float64},
        seeding_devices_per_m2::Float64
    )::Tuple{YAXArray,Matrix{Float64}}

Calculate proportion of deployed corals to be seeded at each of the selected locations.
Distributes seeded corals according to current available space at each selected site.

# Arguments
- `seed_loc_k_m²` : Carrying capacity area of locations to seed in m².
- `available_space` : Currently available space at each seed location in m².
- `seed_volume` : Absolute number of coral to deploy of each functional group.
- `colony_areas` : Area of single 1yo colony of each functional group.
- `seeding_devices_per_m2` : Seeding device density (number of devices per m²).

# Returns
- YAXArray[taxa to seed ⋅ number of seed locations], Proportional increase in cover relative
to locations' `k` area
- Matrix[seed locations ⋅ taxa to seed], Number of coral deployed
"""
function distribute_seeded_corals(
    seed_loc_k_m²::Vector{Float64},
    available_space::Vector{Float64},
    seed_volume::Vector{Float64},
    colony_areas::Vector{Float64},
    seeding_devices_per_m2::Float64
)::Tuple{YAXArray,Matrix{Float64}}
    total_available_space::Float64 = sum(available_space)
    seed_volume_tmp = deepcopy(seed_volume)

    # If n_devices > max_n_devices , cap seed_volume
    max_n_devices = seeding_devices_per_m2 * total_available_space
    n_devices = sum(seed_volume_tmp)            # Assuming one coral per device survives to 1-yr old
    if n_devices > max_n_devices
        seed_volume_tmp .*= (max_n_devices / n_devices)
        cap = n_devices - sum(seed_volume_tmp)
        @warn """
        Number of seeding devices exceeds available space.
        Excluding $cap devices to fit.
        Total available space: $(round(total_available_space))
        Max n devices = $max_n_devices
        N_devices = $n_devices
        """
    end

    # Extract colony areas and determine approximate seeded area in m^2
    seeded_area = colony_areas .* seed_volume_tmp
    total_seeded_area::Float64 = sum(seeded_area)

    # Proportion of available space on each site relative to available space at these
    # locations
    prop_area_avail = available_space ./ total_available_space
    if total_seeded_area > total_available_space
        if !is_test_env()
            @warn "Seeded area exceeds available space. Restricting to available space."
        end

        seeded_area .*= total_available_space / total_seeded_area

        # Update seed_volume_tmp if seeded_area is capped
        seed_volume_tmp = seeded_area ./ colony_areas
    end

    # Distribute seeded corals (as area) across locations according to available space
    # proportions:
    #     proportion available space * (area of 1 coral * num seeded corals)
    # Convert to relative cover proportion by dividing by location area
    relative_seeded_area = ((prop_area_avail .* seeded_area') ./ seed_loc_k_m²)'

    proportional_increase = DataCube(
        relative_seeded_area;
        taxa=collect(caxes(seeded_area)[1].val),
        locations=1:length(available_space)
    )

    n_deployed_coral = prop_area_avail .* seed_volume_tmp'

    @assert sum(n_deployed_coral) ≈ sum(seed_volume_tmp)

    return proportional_increase, n_deployed_coral
end

"""
    update_tolerance_distribution!(
        scaled_seed::YAXArray,
        cover::AbstractArray{T},
        c_dist_t::AbstractArray{T},
        c_mean_reference::AbstractArray{T,2},
        stdev::AbstractArray{T},
        seed_locs::Vector{Int64},
        seed_sc::AbstractMatrix{Bool},
        a_adapt::AbstractVector{T}
    )::Nothing where {T<:Float64}

Update the thermal tolerance distribution of a population due to seeding.
Updates the store `c_dist_t` in place.

# Arguments
- `scaled_seed` : Seeding values transformed to proportion cover increase relative to k area.
- `cover` : Current coral cover state.
- `c_dist_t` : Critical DHW distributions of corals to update (i.e., for time \$t\$).
- `c_mean_reference` : Matrix with c_mean to be used as reference for seeding corals heat
tolerance enhancement.
- `stdev` : Standard deviation of DHW tolerance distributions for each functional type.
- `seed_locs` : Seeding locations
- `seed_sc` : Size classes to seed
- `a_adapt` : Level of assisted adaptation (as DHW tolerance increase)
"""
function update_tolerance_distribution!(
    scaled_seed::YAXArray,
    cover::AbstractArray{T},
    c_dist_t::AbstractArray{T},
    c_mean_reference::AbstractArray{T,2},
    stdev::AbstractArray{T},
    seed_locs::Vector{Int64},
    seed_sc::AbstractMatrix{Bool},
    a_adapt::AbstractVector{T}
)::Nothing where {T<:Float64}

    # Calculate distribution weights using proportion of area (used as priors for MixtureModel)
    # Note: It is entirely possible for a location to be ranked in the top N, but
    #       with no deployments (for a given species). A location with 0 cover
    #       and no deployments will therefore be NaN due to zero division.
    #       These are replaced with 1.0 so that the distribution for unseeded
    #       corals are used.
    w_taxa::Matrix{Float64} = scaled_seed ./ (cover[seed_sc, seed_locs] .+ scaled_seed)
    # NaN occurs when both cover and scaled_seed are 0 (undeployed taxa at bare locations).
    replace!(w_taxa, NaN => 0.0)

    # Update critical DHW distribution for deployed size classes
    a_adapt_relative = copy(a_adapt)
    for (i, loc) in enumerate(seed_locs)
        a_adapt_relative .= (a_adapt .+ c_mean_reference[:, loc])
        # Previous distributions
        c_dist_ti = @view(c_dist_t[seed_sc, loc])

        # Truncated normal distributions for deployed corals
        # Assume same stdev and bounds as original
        tn::Vector{Float64} =
            truncated_normal_mean.(
                a_adapt_relative, stdev[seed_sc], 0.0, a_adapt_relative .+ HEAT_UB
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
