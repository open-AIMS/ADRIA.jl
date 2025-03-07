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

"""
    seeding_intervention!(ctx::SimulationContext, tstep::Int64, leftover_space_m²::Vector{Float64})::Vector

Process the seeding intervention phase.
"""
function seeding_intervention!(
    ctx::SimulationContext, tstep::Int64, leftover_space_m²::Vector{Float64}
)::Vector
    selected_seed_ranks = []

    # Select locations for seeding
    if ctx.is_guided && ctx.seed_decision_years[tstep]
        # Check for locations with available space
        locs_with_space = vec(leftover_space_m²) .> 0.0
        _valid_locs = ctx.habitable_locs .& ctx.depth_criteria
        considered_locs = findall(_valid_locs .& locs_with_space)

        if length(considered_locs) > 0
            selected_seed_ranks = select_seed_locations_guided(
                ctx, considered_locs, locs_with_space, leftover_space_m²
            )

            # Log rankings
            if !isempty(selected_seed_ranks)
                ctx.log_location_ranks[tstep, At(selected_seed_ranks), At(:seed)] .=
                    1:length(selected_seed_ranks)
            end
        end
    elseif ctx.apply_seeding && ctx.seed_decision_years[tstep]
        # Unguided deployment
        selected_seed_ranks = select_locations_unguided(ctx, leftover_space_m²)

        ctx.log_location_ranks[tstep, At(selected_seed_ranks), At(:seed)] .= 1.0
    end

    # Apply seeding if locations were selected
    if !isempty(selected_seed_ranks)
        seed_locations!(ctx, tstep, selected_seed_ranks, leftover_space_m²)
    end

    return selected_seed_ranks
end

"""
    select_seed_locations_guided(ctx::SimulationContext, considered_locs::Vector{Int64},
                               locs_with_space::BitVector, leftover_space_m²::Vector{Float64})

Select locations for seeding using guided approach.
"""
function select_seed_locations_guided(ctx::SimulationContext,
    considered_locs::Vector{Int64},
    locs_with_space::BitVector, leftover_space_m²::Vector{Float64})

    # Prepare environmental projections
    tstep::Int64 = ctx.current_tstep + 1
    decay::Vector{Float64} = 0.99 .^ (1:(ctx.plan_horizon + 1)) .^ 2

    # Use modified projected DHW (may have been affected by fogging or shading)
    dhw_p = copy(ctx.dhw_scen)
    dhw_p[tstep, :] .= ctx.dhw_t

    dhw_projection = weighted_projection(dhw_p, tstep, ctx.plan_horizon, decay, ctx.tf)
    wave_projection = weighted_projection(
        ctx.wave_scen, tstep, ctx.plan_horizon, decay, ctx.tf
    )

    # Calculate connectivity strengths
    vec_abs_k = loc_k_area(ctx.domain)
    loc_coral_cover = _loc_coral_cover(ctx.C_cover_t)
    area_weighted_conn = ctx.domain.conn .* vec_abs_k
    conn_cache = similar(area_weighted_conn.data)
    in_conn, out_conn, _ = connectivity_strength(
        area_weighted_conn, vec(loc_coral_cover), conn_cache
    )

    # Update decision matrix with current criteria values
    _valid_locs = ctx.habitable_locs .& ctx.depth_criteria
    update_criteria_values!(
        ctx.decision_mat;
        heat_stress=dhw_projection[_valid_locs],
        wave_stress=wave_projection[_valid_locs],
        coral_cover=loc_coral_cover[_valid_locs],
        in_connectivity=in_conn[_valid_locs],
        out_connectivity=out_conn[_valid_locs]
    )

    # Select locations based on criteria
    return select_locations(
        ctx.seed_pref,
        ctx.decision_mat[location=locs_with_space[_valid_locs]],
        ctx.MCDA_approach,
        ctx.domain.loc_data.cluster_id,
        ctx.max_area_to_seed,
        considered_locs,
        vec(leftover_space_m²),
        ctx.min_iv_locs,
        ctx.max_members
    )
end

"""
    select_locations_unguided(ctx::SimulationContext, leftover_space_m²::Vector{Float64})

Select locations for intervention using unguided approach.
"""
function select_locations_unguided(
    ctx::SimulationContext, leftover_space_m²::Vector{Float64}
)
    return unguided_selection(
        ctx.domain.loc_ids,
        ctx.min_iv_locs,
        vec(leftover_space_m²),
        ctx.depth_criteria
    )
end

"""
    seed_locations!(
        ctx::SimulationContext, tstep::Int64, selected_seed_ranks::Vector,
        leftover_space_m²::Vector{Float64}
    )

Apply seeding to the selected locations.
"""
function seed_locations!(
    ctx::SimulationContext, tstep::Int64, selected_seed_ranks::Vector,
    leftover_space_m²::Vector{Float64}
)
    # Find selected locations
    seed_locs = findall(ctx.log_location_ranks.locations .∈ [selected_seed_ranks])

    # Filter to locations with available space
    available_space = leftover_space_m²[seed_locs]
    locs_with_space = findall(available_space .> 0.0)

    # Skip if no locations with space
    if length(locs_with_space) == 0
        return nothing
    end

    # Get locations with space
    seed_locs = seed_locs[locs_with_space]

    # Get seed volume
    taxa_names = collect(ctx.param_set.factors[occursin.("N_seed_", ctx.param_set.factors)])
    seed_volume = ctx.param_set[At(taxa_names)]

    # Calculate proportion to seed based on current available space
    vec_abs_k = loc_k_area(ctx.domain)
    proportional_increase, n_corals_seeded = distribute_seeded_corals(
        vec_abs_k[seed_locs],
        available_space,
        ctx.max_seeded_area,
        seed_volume.data
    )

    # Log estimated number of corals seeded
    ctx.seeding_log[tstep, :, seed_locs] .= n_corals_seeded'

    # Add coral seeding to recruitment
    ctx.recruitment[[1, 2, 4], seed_locs] .+= proportional_increase

    # Update tolerance distributions
    return update_tolerance_distribution!(
        proportional_increase,
        ctx.C_cover_t,
        ctx.c_mean_t,
        ctx.c_std,
        seed_locs,
        ctx.seed_sc,
        ctx.a_adapt
    )
end
