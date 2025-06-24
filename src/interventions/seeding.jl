using StatsBase

"""
    distribute_seeded_corals(strategy::Symbol, seed_locs::Vector{Int64}, loc_k_m²::Vector{Float64}, available_space::Vector{Float64}, max_seeded_area::YAXArray, seed_volume::Vector{Float64}, target_density::Float64)

Calculate proportion of deployed corals to be seeded at each of the selected locations.
Distributes seeded corals according to current available space at each selected site.

# Arguments
- strategy : coral seeding strategy
- seed_locs : ordered list of location indices, ordered by intervention priority.
- loc_k_m² : Carrying capacity area of locations to seed in m².
- available_space : Currently available space at each seed location in m².
- max_seeded_area : Area (in m²) of each coral type to be seeded with dim taxa.
- seed_volume : Absolute number of coral to deploy. Aligns with max_seeded_area.
- target_density : density to seed corals n_corals / m^2

# Returns
- YAXArray[taxa to seed ⋅ number of seed locations], Proportional increase in cover relative to locations' `k` area
- Matrix[seed locations ⋅ taxa to seed], Number of coral deployed
"""
function distribute_seeded_corals(
    strategy::Symbol,
    seed_locs::Vector{Int64},
    loc_k_m²::Vector{Float64},
    available_space::Vector{Float64},
    max_seeded_area::YAXArray,
    seed_volume::Vector{Float64},
    target_density::Float64
)
    available_space = available_space[seed_locs]
    n_iv_locs = length(seed_locs)
    no_dep_mask = findall(seed_volume .!= 0.0)
    maximum_density::Vector{Float64} = ones(size(seed_volume))
    maximum_density[no_dep_mask] = seed_volume[no_dep_mask] ./ max_seeded_area[no_dep_mask]

    if strategy == :VARY_LOCATIONS
        target_density, seed_volume, n_iv_locs = vary_locations(
            available_space, target_density, seed_volume, n_iv_locs
        )
    elseif strategy == :VARY_N_SEEDED
        target_density, seed_volume, n_iv_locs = vary_n_corals(
            available_space, target_density, seed_volume, n_iv_locs
        )
    elseif strategy == :VARY_SEED_DENSITY
        target_density, seed_volume, n_iv_locs = vary_seed_density(
            available_space, target_density, seed_volume, n_iv_locs
        )
    elseif strategy == :CAP_DENSITY
        target_density, seed_volume, n_iv_locs = seed_cap_density(
            available_space, target_density, seed_volume, n_iv_locs
        )
    end

    if any(maximum_density .< target_density)
        new_density = minimum(maximum_density)
        msg = "Attempting to deploy target at a density not theoretically possible."
        msg += "Density: $(target_density). Limiting to $(new_density)."
        @warn msg
        seed_volume *= new_density / target_density
        target_density = new_density
    end

    seed_locs = seed_locs[1:n_iv_locs]
    available_space = available_space[1:n_iv_locs]
    loc_k_m² = loc_k_m²[seed_locs]

    # seeded area refers to the real area covered by a seeded coral
    # seeded area = n_corals / density (/m^2)
    seeded_area = seed_volume ./ maximum_density
    # total area cover by seeded corals
    total_seeded_area = sum(seeded_area)
    total_seeded_corals = sum(seed_volume)
    total_available_space::Float64 = sum(available_space)

    # Proportion of available space on each site relative to available space at these
    # locations
    prop_area_avail = available_space ./ total_available_space
    if (total_seeded_corals / total_available_space) > target_density + 1e-7
        @warn "Seeded area exceeds available space. Restricting to 95% of available space."
        max_seeded_area = copy(max_seeded_area)
        seed_volume = copy(seed_volume)
        rescale_coef = total_available_space / (total_seeded_corals / target_density) * 0.95
        max_seeded_area .*= rescale_coef
        seed_volume .*= rescale_coef
    end

    # Distribute seeded corals (as area) across locations according to available space
    # proportions:
    #     proportion * (area of 1 coral * num seeded corals)
    # Convert to relative cover proportion by dividing by location area
    scaled_seed = ((prop_area_avail .* seeded_area') ./ loc_k_m²)'
    proportional_increase = DataCube(
        scaled_seed;
        taxa=caxes(max_seeded_area)[1].val.data,
        locations=1:length(available_space)
    )

    n_deployed_coral = prop_area_avail .* seed_volume'

    # This assumes each taxa is deployed with the same density, variable density over taxa to be added
    n_taxa_seed = size(n_deployed_coral, 2)

    return proportional_increase, n_deployed_coral, seed_locs, target_density/n_taxa_seed
end

"""Find the number of corals required to satisfy the target density and available space."""
function vary_n_corals(
    ordered_avail_areas::Vector{Float64},
    target_density::Float64,
    n_corals::Vector{Float64},
    n_iv_locs::Int64
)
    total_available_space = sum(ordered_avail_areas)
    coral_required = total_available_space * target_density
    # Preserve the ratio type of corals deplored.
    n_corals = n_corals ./ sum(n_corals) .* coral_required
    return target_density, n_corals, n_iv_locs
end

"""
Find the number of corals required to satisfy the target density and available space.

Limit the number of corals to the n_corals given. Deployments may decide to not deploy all
corals if there is not enough space, and will deploy at a lower density if there is not
enough corals.
"""
function seed_cap_density(
    ordered_avail_area::Vector{Float64},
    target_density::Float64,
    n_corals::Vector{Float64},
    n_iv_locs::Int64
)
    total_available_space = sum(ordered_avail_area)
    total_corals = sum(n_corals)
    implied_density = total_corals / total_available_space
    if implied_density == 0.0
        throw(ArgumentError("Attempting to seed 0 corals."))
    end
    if implied_density > target_density
        scale_coef = target_density / implied_density
        n_corals .*= scale_coef
        implied_density = target_density
    end
    return implied_density, n_corals, n_iv_locs
end

"""Find the outplant density required to satsfy the outplant volume and number of locations."""
function vary_seed_density(
    ordered_avail_areas::Vector{Float64},
    target_density::Float64,
    n_corals::Vector{Float64},
    n_iv_locs::Int64
)::Tuple{Float64,Vector{Float64},Int64}
    # Total available space
    sum_avail = sum(ordered_avail_areas[1:n_iv_locs])
    return sum(n_corals) / sum_avail, n_corals, n_iv_locs
end

"""Find the number of locations required to meet the target density."""
function vary_locations(
    ordered_avail_areas::Vector{Float64},
    target_density::Float64,
    n_corals::Vector{Float64},
    n_iv_locs::Int64
)::Tuple{Float64,Vector{Float64},Int64}
    n_corals_sum = sum(n_corals)
    cum_avail = cumsum(ordered_avail_areas)
    target_area = n_corals_sum / target_density
    n_iv_locs = findfirst(x -> x > target_area, cum_avail)

    if isnothing(n_iv_locs)
        new_density = n_corals_sum / cum_avail[end]
        n_iv_locs = length(ordered_avail_areas)

        # For consistency, warn if new density is updated to more than 10% greater than the original value
        if (new_density > target_density * (1.10))
            @warn "Density has been updated to accommodate coral volume to $(new_density) corals/m2."
        end
        return new_density, n_corals, n_iv_locs
    end

    if n_iv_locs!=1
        new_density = n_corals_sum ./ cum_avail[n_iv_locs-1:n_iv_locs]
        n_iv_locs = (n_iv_locs-1:n_iv_locs)[argmin(abs.(new_density.-target_density))]
    end

    target_density = n_corals_sum ./ cum_avail[n_iv_locs]

    return target_density, n_corals, n_iv_locs
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
        tn::Vector{Float64} =
            truncated_normal_mean.(
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
