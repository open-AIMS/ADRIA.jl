"""
    _grow_corals!(C_cover, C_cover_t, recruitment, gs, tstep)

Growth phase of a timestep: take the cover of `tstep - 1` into `C_cover_t` (absolute area),
add the recruits settled last step, set the cover-dependent growth constraints per location,
run CoralBlox for each habitable location, constrain cover that outgrows the habitable
area, and write the result back as relative cover into `C_cover_t` and `C_cover[tstep]`.
`gs` holds the per-scenario growth state built before the time loop (arrays are mutated
in place, nothing is rebound).
"""
function _grow_corals!(C_cover, C_cover_t, recruitment, gs, tstep::Int64)
    (;
        # Habitable locations and their areas
        habitable_locs, habitable_loc_idxs, habitable_loc_areas, habitable_loc_areas′,
        vec_abs_k,

        # Carrying capacity (maximum projected cover)
        habitable_max_projected_cover, cache_habitable_max_projected_cover,

        # Recruitment scratch space
        _recruitment_col_sum_2d, _recruitment_col_sum,

        # Growth-constraint masks and caches
        agg_cover_above_threshold_mask, growth_threshold_mask_cache, cover_threshold_mask,
        relative_habitable_cover_cache, growth_constraints, net_growth_rates,

        # Cover thresholds for the growth transition
        cover_transition_lb, cover_transition_ub,

        # Calibration groups (per-location masks and indices) and their scaled rates
        cb_calib_group_masks, n_cb_calib_groups, loc_cb_calib_group_idxs, biogrp_lin_ext,
        biogrp_survival, _le_cache,

        # Growth acceleration curve parameters
        growth_acc_height, growth_acc_midpoint, growth_acc_steepness,

        # Functional group definitions
        functional_groups
    ) = gs

    # Convert cover to absolute values to use within CoralBlox model
    C_cover_t[:, :, habitable_locs] .=
        C_cover[tstep - 1, :, :, habitable_locs] .* habitable_loc_areas′

    # To prevent overgrowth when adding recruitment to C_cover_t, set a threshold above
    # which recruits are either set to zero or rescaled.
    constrain_recruitment!(
        recruitment, agg_cover_above_threshold_mask, C_cover_t, habitable_loc_areas
    )

    # Settlers from t-1 grow into observable sizes.
    C_cover_t[:, 1, habitable_locs] .+= recruitment

    sum!(_recruitment_col_sum_2d, recruitment)
    habitable_max_projected_cover .=
        cache_habitable_max_projected_cover .+ _recruitment_col_sum

    relative_habitable_cover_cache[habitable_locs] .=
        loc_coral_cover(C_cover_t[:, :, habitable_locs]) ./
        vec_abs_k[habitable_locs]

    # Growth constrains need to be calculated seperately for differen growth rates
    growth_threshold_mask_cache .=
        relative_habitable_cover_cache .>= cover_transition_ub
    for idx = 1:n_cb_calib_groups
        cover_threshold_mask .=
            cb_calib_group_masks[:, idx] .&& growth_threshold_mask_cache

        # Only apply linear_extension_scale_factors to locations with high cover
        @views growth_constraints[cover_threshold_mask] .= linear_extension_scale_factors(
            C_cover_t[:, :, cover_threshold_mask],
            vec_abs_k[cover_threshold_mask],
            biogrp_lin_ext[:, :, idx],
            _le_cache,
            habitable_max_projected_cover[cover_threshold_mask]
        )
    end

    # When cover is between cover_transition_lb and cover_transition_ub use a weighted mean
    @. growth_threshold_mask_cache = (
        cover_transition_lb <= relative_habitable_cover_cache <
        cover_transition_ub
    )
    if sum(growth_threshold_mask_cache) > 0
        for idx = 1:n_cb_calib_groups
            cover_threshold_mask .=
                growth_threshold_mask_cache .&& cb_calib_group_masks[:, idx]
            transition_scale =
                (
                    relative_habitable_cover_cache[cover_threshold_mask] .-
                    cover_transition_lb
                ) ./ (cover_transition_ub - cover_transition_lb)
            growth_constraints[cover_threshold_mask] .=
                (
                    transition_scale .*
                    linear_extension_scale_factors(
                        C_cover_t[:, :, cover_threshold_mask],
                        vec_abs_k[cover_threshold_mask],
                        biogrp_lin_ext[:, :, idx],
                        _le_cache,
                        habitable_max_projected_cover[cover_threshold_mask]
                    )
                ) .+
                (
                    (1 .- transition_scale) .*
                    growth_acceleration.(
                        growth_acc_height[idx],
                        growth_acc_midpoint[idx],
                        growth_acc_steepness[idx],
                        relative_habitable_cover_cache[cover_threshold_mask]
                    )
                )
        end
    end
    growth_threshold_mask_cache .=
        relative_habitable_cover_cache .< cover_transition_lb
    for idx = 1:n_cb_calib_groups
        cover_threshold_mask .=
            growth_threshold_mask_cache .&& cb_calib_group_masks[:, idx]
        growth_constraints[cover_threshold_mask] .= growth_acceleration.(
            growth_acc_height[idx],
            growth_acc_midpoint[idx],
            growth_acc_steepness[idx],
            relative_habitable_cover_cache[cover_threshold_mask]
        )
    end

    @inbounds for i in habitable_loc_idxs
        cb_idx::Int64 = loc_cb_calib_group_idxs[i]
        net_growth_rates[:, :, i] .=
            @view(biogrp_lin_ext[:, :, cb_idx]) .*
            growth_constraints[i]

        timestep!(
            functional_groups[i],
            @view(recruitment[:, i]),
            @view(net_growth_rates[:, :, i]),
            @view(biogrp_survival[:, :, cb_idx])
        )

        # Write to the cover matrix
        coral_cover(functional_groups[i], @view(C_cover_t[:, :, i]))
    end

    # Check if size classes are inappropriately out-growing habitable area
    if any(loc_coral_cover(C_cover_t)[habitable_locs] .> habitable_loc_areas)
        @warn "Cover outgrowing habitable area at tstep $tstep. Constraining."
        outgrowing_locs_mask =
            loc_coral_cover(C_cover_t)[habitable_locs] .> habitable_loc_areas
        C_cover_t[:, :, outgrowing_locs_mask] .*=
            reshape(
                vec_abs_k[outgrowing_locs_mask] ./
                loc_coral_cover(C_cover_t)[outgrowing_locs_mask],
                (1, 1, count(outgrowing_locs_mask))
            ) .* 0.999
    end

    # Convert C_cover_t to relative values after CoralBlox was run
    C_cover_t[:, :, habitable_locs] .= (
        @view(C_cover_t[:, :, habitable_locs]) ./ habitable_loc_areas′
    )
    C_cover[tstep, :, :, habitable_locs] .= @view(C_cover_t[:, :, habitable_locs])

    return nothing
end

"""
    _adapt_and_log!(
        C_cover, c_mean_t, net_growth_rates, c_mean_tol_ceil, dhw_tol_mean_log, cover_log,
        _permuted_buf, log_dhw_tols, log_cover, tstep, tf
    )

Natural adaptation of the thermal tolerance means `c_mean_t` (cover is not changed), then
the optional per-timestep tolerance and cover logs.
"""
function _adapt_and_log!(
    C_cover, c_mean_t, net_growth_rates, c_mean_tol_ceil, dhw_tol_mean_log, cover_log,
    _permuted_buf, log_dhw_tols::Bool, log_cover::Bool, tstep::Int64, tf::Int64
)
    # Natural adaptation (doesn't change C_cover_t)
    if tstep <= tf
        adjust_DHW_distribution!(
            @view(C_cover[tstep - 1, :, :, :]), c_mean_t, net_growth_rates,
            c_mean_tol_ceil
        )

        if log_dhw_tols
            # Log dhw tolerances if requested
            permutedims!(_permuted_buf, c_mean_t, (2, 1, 3))
            dhw_tol_mean_log[tstep, :, :] .= reshape(
                _permuted_buf, size(dhw_tol_mean_log)[2:3]
            )
        end

        if log_cover
            permutedims!(_permuted_buf, @view(C_cover[tstep, :, :, :]), (2, 1, 3))
            cover_log[tstep, :, :] .= reshape(_permuted_buf, size(cover_log)[2:3])
        end
    end

    return nothing
end

"""
    _reproduce!(recruitment, C_cover_t, c_mean_t, rs, rng, apply_allee_effect)

Reproduction phase: snapshot the juvenile tolerance (before coral aquaculture changes it),
compute fecundity scope and its proportions per location, the leftover space, and the
settlers (`settler_cover`, the only RNG draws of this phase) as the new `recruitment`
(relative to the location area). `leftover_space_m²` is reduced by the settlers.
`rs` holds the per-scenario buffers and constants built before the time loop.
"""
function _reproduce!(
    recruitment, C_cover_t, c_mean_t, rs, rng, apply_allee_effect::Bool
)
    (;
        # Locations and areas
        n_locs, habitable_loc_idxs, habitable_areas, vec_abs_k, _loc_coral_cover,

        # Fecundity and settlement
        fec_scope, fecundity_per_m², prop_fecundity, potential_settlers,
        basal_area_per_settler, conn, sim_params,

        # Space and recruitment scratch space
        leftover_space_m², _recruitment_col_sum_2d, _recruitment_col_sum,

        # Heat tolerance snapshot taken before coral aquaculture
        c_mean_pre_CAq
    ) = rs

    # Snapshot juvenile tolerance after natural adaptation, before coral aquaculture modifies it.
    # This is the correct baseline for a_adapt enhancement: using the post-coral-aquaculture state
    # here would cause the enhancement to compound each deployment year.
    c_mean_pre_CAq .= c_mean_t[:, JUVENILE_SIZE_CLASS, :]

    # Reproduction
    # Calculates scope for coral fedundity for each size class and at each location
    fecundity_scope!(
        fec_scope, fecundity_per_m², C_cover_t, habitable_areas;
        apply_allee_effect=apply_allee_effect
    )

    for l = 1:n_locs
        s = sum(fec_scope[:, l])
        if s > 0.0
            prop_fecundity[:, l] .= fec_scope[:, l] ./ s
        else
            prop_fecundity[:, l] .= 0.0
        end
    end

    _loc_coral_cover .= loc_coral_cover(C_cover_t)
    leftover_space_m² .= relative_leftover_space(_loc_coral_cover) .* vec_abs_k

    # Reset potential settlers to zero
    potential_settlers .= 0.0
    recruitment .= 0.0

    # habitable_areas_view = @view habitable_areas[:, habitable_locs]

    # Recruitment represents additional cover, relative to total location area
    # Recruitment/settlement occurs after the full moon in October/November
    @views recruitment[:, habitable_loc_idxs] .=
        settler_cover(
            fec_scope,
            conn,
            leftover_space_m²,
            sim_params.max_settler_density,
            sim_params.max_larval_density,
            basal_area_per_settler,
            potential_settlers;
            rng=rng
        )[
            :, habitable_loc_idxs
        ] ./ habitable_areas[:, habitable_loc_idxs]

    # Reset fecundity scope before next run
    fec_scope .= 0.0

    sum!(_recruitment_col_sum_2d, recruitment)
    leftover_space_m² .-= _recruitment_col_sum .* vec_abs_k

    return nothing
end

"""
    _apply_shading!(
        dhw_t, YShd, dhw_scen, Shd_locs_mask, Shd_decision_years, apply_Shd, Shd, tstep
    )

Take the DHW of `tstep` into `dhw_t` and, in shading years, reduce it at the shaded
locations (logged in `YShd`). This is the surface DHW used by everything that follows.
"""
function _apply_shading!(
    dhw_t, YShd, dhw_scen, Shd_locs_mask, Shd_decision_years, apply_Shd::Bool,
    Shd::Float64, tstep::Int64
)
    # Shading
    # Apply regional cooling effect before selecting locations for coral aquaculture
    dhw_t .= dhw_scen[tstep, :]  # subset of DHW for given timestep
    if apply_Shd && Shd_decision_years[tstep]
        YShd[tstep, :] .= Shd

        # Apply reduction in DHW due to shading
        dhw_t[Shd_locs_mask] .= max.(0.0, dhw_t[Shd_locs_mask] .- Shd)
    end

    return nothing
end

"""
    _update_guided_state!(
        dhw_projection, wave_projection, current_loc_cover, dhw_p, dhw_scen, dhw_t,
        wave_scen, wave_idx, area_weighted_conn, _loc_coral_cover, conn_cache, C_cover_t,
        plan_horizon, decay, tf, tstep
    )

Guided scenarios only. Update, in place, the projected DHW (with the shaded DHW of this
step) and wave stress, and the current location cover used by the location selection.
Returns the in and out connectivity strengths for the decision matrices.
"""
function _update_guided_state!(
    dhw_projection, wave_projection, current_loc_cover, dhw_p, dhw_scen, dhw_t, wave_scen,
    wave_idx, area_weighted_conn, _loc_coral_cover, conn_cache, C_cover_t, plan_horizon,
    decay, tf, tstep
)
    # Use modified projected DHW (may have been affected by shading)
    dhw_p .= dhw_scen
    dhw_p[tstep, :] .= dhw_t
    dhw_projection .= weighted_projection(dhw_p, tstep, plan_horizon, decay, tf)

    if wave_idx > 0.0
        wave_projection .= weighted_projection(
            wave_scen, tstep, plan_horizon, decay, tf
        )
    end

    # Determine connectivity strength weighting by area.
    # Accounts for strength of connectivity where there is low/no coral cover
    in_conn, out_conn, _ = connectivity_strength(
        area_weighted_conn, vec(_loc_coral_cover), conn_cache
    )

    # Calculate current location cover
    current_loc_cover .= dropdims(sum(C_cover_t; dims=(1, 2)); dims=(1, 2))

    return in_conn, out_conn
end

"""
    _apply_fogging!(fog, ivs, in_conn, out_conn, dhw_t, tstep)

Fogging: select the fogging locations of this step (guided by the decision matrix, or
unguided at random using the scenario RNG) and reduce `dhw_t` there. `fog` holds the
fogging strategy and its scenario constants, `ivs` the state shared by the intervention
phases (both built before the time loop). In and out connectivity are `nothing` unless
the scenario is guided.
"""
function _apply_fogging!(fog, ivs, in_conn, out_conn, dhw_t, tstep::Int64)
    (;
        # Strategy and location selection
        Fog_strategy, Fog_pref, Fog_decision_mat, min_iv_locs,

        # Deployment settings and tracking
        fogging, last_Fog_deployment, Yfog
    ) = fog
    (;
        # Domain, run mode and random numbers
        domain, is_guided, MCDA_approach, rng,

        # Projections and current state
        current_loc_cover, dhw_projection, wave_projection, leftover_space_m²,

        # Intervention tracking
        recent_cover_losses, log_location_ranks, iv_CAq_last_deployment
    ) = ivs

    # Fog location selection
    if !isnothing(Fog_strategy)
        state = if is_guided
            build_state(
                domain, Fog_strategy,
                (
                    current_cover=current_loc_cover,
                    recent_cover_losses=recent_cover_losses,
                    last_deployment=iv_CAq_last_deployment
                )
            )
        else
            nothing
        end

        # Get candidate locations from strategy
        candidate_locs = filter_candidate_locations(Fog_strategy, tstep, state)
        candidate_loc_indices = findall(
            in.(domain.loc_ids, Ref(candidate_locs))
        )

        if !isempty(candidate_locs)
            if is_guided
                # Update decision matrix with current conditions
                update_criteria_values!(
                    Fog_decision_mat[location = At(candidate_locs)];
                    heat_stress=dhw_projection[candidate_loc_indices],
                    wave_stress=wave_projection[candidate_loc_indices],
                    coral_cover=current_loc_cover[candidate_loc_indices],
                    in_connectivity=in_conn[candidate_loc_indices],
                    out_connectivity=out_conn[candidate_loc_indices]
                )

                # Build state for target locations only
                selected_Fog_ranks = select_locations(
                    Fog_pref,
                    Fog_decision_mat[location = At(candidate_locs)],
                    MCDA_approach,
                    min_iv_locs
                )
            else
                selected_Fog_ranks = unguided_selection(
                    candidate_locs,
                    min_iv_locs,
                    vec(leftover_space_m²[candidate_loc_indices]);
                    rng=rng
                )
            end
            if !isempty(selected_Fog_ranks)
                log_val = is_guided ? (1:length(selected_Fog_ranks)) : 1.0
                log_location_ranks[tstep, At(selected_Fog_ranks), At(:fog)] .=
                    log_val
                selected_Fog_loc_idx = findall(
                    in.(domain.loc_ids, Ref(selected_Fog_ranks))
                )
                last_Fog_deployment[selected_Fog_loc_idx] .= tstep

                # Fog selected locations
                Fog_locs = findall(log_location_ranks.locations .∈ [selected_Fog_ranks])
                Fog_locations!(@view(Yfog[tstep, :]), Fog_locs, dhw_t, fogging)
            end
        end
    end

    return nothing
end

"""
    _apply_larval_methods!(
        lvm, ivs, in_conn, out_conn, recruitment, prop_fecundity, param_set, tstep
    )

Larval methods: for each target share, select the locations (guided, or unguided at
random using the scenario RNG) and add the settlers to `recruitment`, reducing the leftover
space. `lvm` holds the strategy and scenario constants, `ivs` the state shared by the
intervention phases (both built before the time loop). In and out connectivity are
`nothing` unless the scenario is guided.
"""
function _apply_larval_methods!(
    lvm, ivs, in_conn, out_conn, recruitment, prop_fecundity, param_set, tstep::Int64
)
    (;
        # Strategy and location selection
        LvM_strategy, LvM_pref, LvM_decision_mat, LvM_min_iv_locs, depth_criteria,

        # Deployment settings and tracking
        iv_LvM_N_settlers, last_LvM_deployment, YLvM,

        # Colony sizes and areas
        colony_areas, _CAq_size_groups, vec_abs_k
    ) = lvm
    (;
        # Domain, run mode and random numbers
        domain, is_guided, MCDA_approach, rng,

        # Projections and current state
        current_loc_cover, dhw_projection, wave_projection, leftover_space_m²,

        # Intervention tracking
        recent_cover_losses, log_location_ranks
    ) = ivs

    # Larval methods intervention
    if !isnothing(LvM_strategy)
        state = if is_guided
            build_state(
                domain, LvM_strategy,
                (
                    current_cover=current_loc_cover,
                    recent_cover_losses=recent_cover_losses,
                    last_deployment=last_LvM_deployment
                )
            )
        else
            nothing
        end

        if is_decision_year(LvM_strategy, tstep)
            # Get candidate locations from strategy
            candidate_locs = filter_candidate_locations(LvM_strategy, tstep, state)

            # LvM_share is a @NamedTuple{weight::Float64, target_locs::Vector{String}}
            for LvM_share in domain.LvM_target_locations
                share_candidate_locs = intersect(
                    candidate_locs, LvM_share.target_locs
                )
                if isempty(share_candidate_locs)
                    if ADRIA.decision.strategy_type(param_set, "iv_LvM") ==
                        PeriodicStrategy
                        @warn """
                            tstep $tstep: Deployment with PeriodicStrategy on
                            $(length(LvM_share.target_locs)) reefs with weight
                            $(LvM_share.weight) skipped because no candidate reefs
                            were found.
                        """
                    end
                    continue
                end

                share_candidate_loc_idx = findall(
                    in.(domain.loc_ids, Ref(share_candidate_locs))
                )

                if is_guided
                    # Update decision matrix with current conditions
                    update_criteria_values!(
                        LvM_decision_mat[location = At(share_candidate_locs)];
                        heat_stress=dhw_projection[share_candidate_loc_idx],
                        wave_stress=wave_projection[share_candidate_loc_idx],
                        coral_cover=current_loc_cover[share_candidate_loc_idx],
                        in_connectivity=in_conn[share_candidate_loc_idx],
                        out_connectivity=out_conn[share_candidate_loc_idx]
                    )

                    selected_LvM_ranks = select_locations(
                        LvM_pref,
                        LvM_decision_mat[location = At(share_candidate_locs)],
                        MCDA_approach,
                        LvM_min_iv_locs
                    )
                else
                    # Unguided deployment, anywhere with available space > 0
                    selected_LvM_ranks = unguided_selection(
                        share_candidate_locs,
                        LvM_min_iv_locs,
                        vec(leftover_space_m²[share_candidate_loc_idx]),
                        depth_criteria[share_candidate_loc_idx];
                        rng=rng
                    )
                end
                if !isempty(selected_LvM_ranks)
                    log_val = is_guided ? (1:length(selected_LvM_ranks)) : 1.0
                    log_location_ranks[tstep, At(selected_LvM_ranks), At(:lvm)] .=
                        log_val
                    selected_LvM_loc_idx = findall(
                        in.(domain.loc_ids, Ref(selected_LvM_ranks))
                    )
                    last_LvM_deployment[selected_LvM_loc_idx] .= tstep
                end

                # Check if locations are selected
                has_LvM_locs::Bool = !isempty(selected_LvM_ranks)

                # Apply Larval Methods (assumed to occur after spawning)
                if has_LvM_locs
                    LvM_loc_idx = findall(domain.loc_ids .∈ [selected_LvM_ranks])

                    @views available_space = leftover_space_m²[LvM_loc_idx]

                    locs_with_space = findall(available_space .> 0.0)

                    # If there are locations with space to select from, deploy what we can
                    if length(locs_with_space) > 0
                        LvM_loc_idx = LvM_loc_idx[locs_with_space]
                        available_space = available_space[locs_with_space]

                        @views LvM_proportional_increase, n_LvM_corals = distribute_LvM_corals(
                            vec_abs_k[LvM_loc_idx],
                            available_space,
                            iv_LvM_N_settlers * LvM_share.weight,
                            colony_areas[_CAq_size_groups],
                            prop_fecundity[:, LvM_loc_idx]
                        )

                        @views recruitment[:, LvM_loc_idx] .+=
                            LvM_proportional_increase
                        leftover_space_m²[LvM_loc_idx] .-=
                            vec(sum(Array(LvM_proportional_increase); dims=1)) .*
                            vec_abs_k[LvM_loc_idx]

                        # Log estimated number of corals moved
                        YLvM[tstep, :, LvM_loc_idx] .= n_LvM_corals
                    end
                end
            end
        end
    end

    return nothing
end

"""
    _apply_coral_aquaculture!(
        caq, ivs, in_conn, out_conn, recruitment, C_cover_t, c_mean_t, param_set, tstep
    )

Coral aquaculture: for each target share, select the locations (guided, or unguided at
random using the scenario RNG), add the deployed corals to `recruitment`, reduce the
leftover space (after larval methods) and mix the deployed tolerance into `c_mean_t`.
`caq` holds the strategy and scenario constants, `ivs` the state shared by the
intervention phases (both built before the time loop). In and out connectivity are
`nothing` unless the scenario is guided.
"""
function _apply_coral_aquaculture!(
    caq, ivs, in_conn, out_conn, recruitment, C_cover_t, c_mean_t, param_set, tstep::Int64
)
    (;
        # Strategy and location selection
        CAq_strategy, CAq_pref, CAq_decision_mat, min_iv_locs, depth_criteria,

        # Deployment settings and tracking
        CAq_volume, CAq_devices_per_m2, YCAq,

        # Colony sizes and areas
        colony_areas, _CAq_size_groups, vec_abs_k,

        # Heat tolerance of the outplanted corals
        c_mean_reference, c_std, a_adapt, c_mean_tol_ceil
    ) = caq
    (;
        # Domain, run mode and random numbers
        domain, is_guided, MCDA_approach, rng,

        # Projections and current state
        current_loc_cover, dhw_projection, wave_projection, leftover_space_m²,

        # Intervention tracking
        recent_cover_losses, log_location_ranks, iv_CAq_last_deployment
    ) = ivs

    # Coral aquaculture location selection
    if !isnothing(CAq_strategy)
        state = if is_guided
            build_state(
                domain,
                CAq_strategy,
                (
                    current_cover=current_loc_cover,
                    recent_cover_losses=recent_cover_losses,
                    last_deployment=iv_CAq_last_deployment
                )
            )
        else
            nothing
        end

        if is_decision_year(CAq_strategy, tstep)
            # Get candidate locations from strategy
            candidate_locs = filter_candidate_locations(CAq_strategy, tstep, state)

            # CAq_share is a @NamedTuple{weight::Float64, target_locs::Vector{String}}
            for CAq_share in domain.CAq_target_locations
                share_candidate_locs = intersect(
                    candidate_locs, CAq_share.target_locs
                )
                if isempty(share_candidate_locs)
                    if ADRIA.decision.strategy_type(param_set, "iv_CAq") ==
                        PeriodicStrategy
                        @warn """
                            tstep $tstep: Deployment with PeriodicStrategy on
                            $(length(CAq_share.target_locs)) reefs with weight
                            $(CAq_share.weight) skipped because no candidate reefs
                            were found.
                        """
                    end
                    continue
                end

                share_candidate_loc_idx = findall(
                    in.(domain.loc_ids, Ref(share_candidate_locs))
                )

                if is_guided
                    # Update decision matrix with current conditions
                    update_criteria_values!(
                        CAq_decision_mat[location = At(share_candidate_locs)];
                        heat_stress=dhw_projection[share_candidate_loc_idx],
                        wave_stress=wave_projection[share_candidate_loc_idx],
                        coral_cover=current_loc_cover[share_candidate_loc_idx],
                        in_connectivity=in_conn[share_candidate_loc_idx],
                        out_connectivity=out_conn[share_candidate_loc_idx]
                    )

                    # Build state for target locations only
                    selected_CAq_ranks = select_locations(
                        CAq_pref,
                        CAq_decision_mat[location = At(share_candidate_locs)],
                        MCDA_approach,
                        min_iv_locs
                    )
                else
                    # Unguided deployment, coral-aquaculture/fog corals anywhere, so long as available space > 0
                    selected_CAq_ranks = unguided_selection(
                        share_candidate_locs,
                        min_iv_locs,
                        vec(leftover_space_m²[share_candidate_loc_idx]),
                        depth_criteria[share_candidate_loc_idx];
                        rng=rng
                    )
                end
                if !isempty(selected_CAq_ranks)
                    log_val = is_guided ? (1:length(selected_CAq_ranks)) : 1.0
                    log_location_ranks[tstep, At(selected_CAq_ranks), At(:caq)] .=
                        log_val
                    selected_CAq_loc_idx = findall(
                        in.(domain.loc_ids, Ref(selected_CAq_ranks))
                    )
                    iv_CAq_last_deployment[selected_CAq_loc_idx] .= tstep
                end

                # Check if locations are selected (can reuse previous selection)
                has_CAq_locs::Bool = !isempty(selected_CAq_ranks)

                # Apply coral aquaculture (assumed to occur after spawning)
                if has_CAq_locs  # CAq_decision_years[tstep] &&
                    # Deploy coral aquaculture at selected locations
                    # Selected locations can fill up over time so avoid locations with no space´
                    CAq_loc_idx = findall(domain.loc_ids .∈ [selected_CAq_ranks])

                    available_space = leftover_space_m²[CAq_loc_idx]
                    locs_with_space = findall(available_space .> 0.0)

                    # If there are locations with space to select from, then deploy what we can
                    # Otherwise, do nothing.
                    if length(locs_with_space) > 0
                        # Calculate proportion for coral aquaculture based on current available space
                        CAq_loc_idx = CAq_loc_idx[locs_with_space]
                        available_space = available_space[locs_with_space]

                        proportional_increase, n_corals_CAq = distribute_CAq_corals(
                            vec_abs_k[CAq_loc_idx],
                            available_space,
                            CAq_volume .* CAq_share.weight,
                            colony_areas[_CAq_size_groups],
                            CAq_devices_per_m2
                        )

                        # Log estimated number of corals deployed via coral aquaculture
                        YCAq[tstep, :, CAq_loc_idx] .= n_corals_CAq'

                        # Add coral aquaculture to recruitment
                        recruitment[:, CAq_loc_idx] .+= proportional_increase
                        leftover_space_m²[CAq_loc_idx] .-=
                            vec(sum(Array(proportional_increase); dims=1)) .*
                            vec_abs_k[CAq_loc_idx]

                        update_tolerance_distribution!(
                            proportional_increase,
                            C_cover_t,
                            c_mean_t,
                            c_mean_reference[:, :, end],
                            c_std,
                            CAq_loc_idx,
                            _CAq_size_groups,
                            a_adapt,
                            c_mean_tol_ceil
                        )
                    end
                end
            end
        end
    end

    return nothing
end

"""
    _apply_disturbances!(
        ds, C_cover_t, ΔC_cover_t, c_mean_t, c_mean_t_1, dhw_t, eff_dhw_t, tstep
    )

Environmental disturbances of a timestep: snapshot cover (`ΔC_cover_t`) and tolerance
(`c_mean_t_1`) after all seeding and before any mortality, then bleaching mortality (on the
depth-attenuated DHW), the tolerance reference update and cyclone mortality.
`ds` holds the per-scenario constants and buffers built before the time loop.
"""
function _apply_disturbances!(
    ds, C_cover_t, ΔC_cover_t, c_mean_t, c_mean_t_1, dhw_t, eff_dhw_t, tstep::Int64
)
    (;
        # Location data and cyclone mortality
        loc_data, cyclone_mortality_scen,

        # DHW and bleaching
        eff_dhw_base, eff_dhw_mix, bleach_dhw, low_dhw_bleach_reduction,
        low_dhw_bleach_threshold,

        # Heat tolerance
        c_std, c_mean_tol_ceil, c_mean_reference, c_mean_pre_CAq, a_adapt_ref
    ) = ds

    # Apply disturbances

    # ΔC_cover_t should only hold changes in cover due to env. disturbances
    ΔC_cover_t .= C_cover_t

    # Store current c_mean before bleaching. This will be used in the next timestep to
    # update the settler's DHW tolerance dists.
    c_mean_t_1 .= c_mean_t

    # Bleaching mortality
    # Calculate and apply bleaching mortality
    # Bleaching typically occurs in the warmer months (November - February)
    # but mortality may occur within this timeframe or some time after it.
    # An earlier attempt to account for the cooling effect of storms / high wave
    # activity involves a linear reduction:
    #    `dhw_t .* (1.0 .- wave_scen[tstep, :])`
    # `wave_scen` is normalized to the maximum value found for the given wave scenario
    # so what causes 100% mortality can differ between runs.
    # Wave activity is said to also promote nutrient cycling, and improve growth.
    # Removed as I could not find any justification behind the implementation.
    eff_dhw_t .= effective_dhw_at_depth.(
        dhw_t, loc_data.depth_med, eff_dhw_base, eff_dhw_mix
    )
    bleaching_mortality!(
        C_cover_t,
        eff_dhw_t,
        c_std,
        c_mean_t,
        @view(bleach_dhw[(tstep - 1):tstep, :, :, :]),
        c_mean_tol_ceil;
        low_dhw_reduction=low_dhw_bleach_reduction,
        low_dhw_threshold=low_dhw_bleach_threshold
    )

    # Store current means to be used in future timesteps.
    # Use pre-coral-aquaculture snapshot so that a_adapt enhancement is not compounded
    # into the reference year-over-year.
    if a_adapt_ref > 0
        c_mean_reference[:, :, 2:end] .= c_mean_reference[:, :, 1:(end - 1)]
    end
    c_mean_reference[:, :, 1] .= c_mean_pre_CAq

    # Coral deaths due to selected cyclone scenario
    # Peak cyclone period is January to March
    # TODO: Update cyclone data to hold data for relevant functional groups
    cyclone_mortality!(C_cover_t, cyclone_mortality_scen[tstep, :, :]')

    return nothing
end

"""
    _step_cots!(cots, Ycots, Ycots_bc, C_cover_t, tstep, year)

One year of the built-in COTS life cycle and its grazing on coral (reduces `C_cover_t`),
then larval dispersal between locations and larvae arriving from outside the domain.
Does nothing when COTS is off (`cots` is `nothing`).
"""
function _step_cots!(cots, Ycots, Ycots_bc, C_cover_t, tstep::Int64, year::Int64)
    # COTS: one year of the COTS life cycle and its grazing on coral, then larval
    # dispersal between locations and larvae arriving from outside the domain
    cots_step!(cots, C_cover_t, tstep, year)
    log_cots!(Ycots, Ycots_bc, cots, tstep)

    return nothing
end

"""
    _step_submodel_hooks!(slots, cover_before, C_cover_t, dhw_t, vec_abs_k, tstep, year)

Submodels attached to the domain, in attach order, at the `AfterDisturbances` phase:
they may only reduce `C_cover_t`, checked per cell against the snapshot `cover_before`
taken just before the phase. Only called when submodels are attached. `@noinline` keeps
the hooks (specialised per submodel combination) out of the compiled time loop.
"""
@noinline function _step_submodel_hooks!(
    slots::Tuple, cover_before, C_cover_t, dhw_t, vec_abs_k, tstep::Int64, year::Int64
)
    copyto!(cover_before, C_cover_t)
    _step_submodels!(
        slots, AfterDisturbances(), CoralStepView(C_cover_t, dhw_t, vec_abs_k), tstep, year
    )
    _check_decrease_only(cover_before, C_cover_t, tstep)
    _record_submodels!(slots, tstep)

    return nothing
end

"""
    _propagate_survival!(
        survival_rate_cache, ΔC_cover_t, C_cover_t, functional_groups, recruitment,
        habitable_areas, C_cover, n_locs, tstep
    )

Survival propagation: the survival ratio of this step's disturbances (`C_cover_t` over its
pre-disturbance snapshot `ΔC_cover_t`, must be <= 1) is applied to CoralBlox and to the
pending `recruitment` (size class 1), which is also converted to absolute area (m²) by
multiplying by `habitable_areas`. The cover is stored in `C_cover[tstep]`.
"""
function _propagate_survival!(
    survival_rate_cache, ΔC_cover_t, C_cover_t, functional_groups, recruitment,
    habitable_areas, C_cover, n_locs::Int64, tstep::Int64
)
    # Calculate survival_rate due to env. disturbances
    @inbounds for i in eachindex(survival_rate_cache)
        d = ΔC_cover_t[i]
        survival_rate_cache[i] = d == 0.0 ? 1.0 : C_cover_t[i] / d
    end
    @assert !any(>(1.0), survival_rate_cache) "Survival rate should be <= 1"

    for loc = 1:n_locs
        apply_mortality!(
            functional_groups[loc], @view(survival_rate_cache[:, :, loc])
        )
    end

    recruitment .*= (view(survival_rate_cache,:,1,:) .* habitable_areas)

    C_cover[tstep, :, :, :] .= C_cover_t

    return nothing
end

"""
    _track_cover_loss!(
        recent_cover_losses, Δcover_loss_proportion, ΔC_cover_t, C_cover_t, n_locs
    )

Proportional cover loss at each location due to this timestep's disturbances, pushed on the
rolling buffer used by reactive strategies. Note that the same `Δcover_loss_proportion`
vector is pushed every step (kept as is).
"""
function _track_cover_loss!(
    recent_cover_losses, Δcover_loss_proportion, ΔC_cover_t, C_cover_t, n_locs::Int64
)
    # Calculate proportional cover loss at each location
    # due to disturbances this timestep
    Δcover_loss_proportion .= 0.0
    for loc = 1:n_locs
        cover_before = sum(@view(ΔC_cover_t[:, :, loc]))
        cover_after = sum(@view(C_cover_t[:, :, loc]))

        if cover_before > 0.0
            Δcover_loss_proportion[loc] = (cover_before - cover_after) / cover_before
        end
    end

    # Store in rolling buffer
    push!(recent_cover_losses, Δcover_loss_proportion)

    return nothing
end
