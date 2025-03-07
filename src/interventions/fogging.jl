function fog_locations!(Yfog, locs, dhw_t, fogging)
    dhw_t[locs] .= dhw_t[locs] .* (1.0 .- fogging)
    return Yfog[locs] .= fogging
end

"""
    select_fog_locations_guided(ctx::SimulationContext)

Select locations for fogging using guided approach.
"""
function select_fog_locations_guided(ctx::SimulationContext)
    return select_locations(
        ctx.fog_pref,
        ctx.decision_mat,
        ctx.MCDA_approach,
        ctx.min_iv_locs
    )
end

"""
    fogging_intervention!(ctx::SimulationContext, tstep::Int64, leftover_space_m²::Vector{Float64})

Process the fogging intervention phase.
"""
function fogging_intervention!(
    ctx::SimulationContext, tstep::Int64, leftover_space_m²::Vector{Float64}
)
    selected_fog_ranks = []

    # Select locations for fogging
    if ctx.is_guided
        if ctx.fog_decision_years[tstep] && (ctx.fogging > 0.0)
            selected_fog_ranks = select_fog_locations_guided(ctx)

            if !isempty(selected_fog_ranks)
                ctx.log_location_ranks[tstep, At(selected_fog_ranks), At(:fog)] .=
                    1:length(selected_fog_ranks)
            end
        end
    elseif ctx.apply_fogging && ctx.fog_decision_years[tstep]
        selected_fog_ranks = select_locations_unguided(ctx, leftover_space_m²)

        ctx.log_location_ranks[tstep, At(selected_fog_ranks), At(:fog)] .= 1.0
    end

    # Apply fogging to selected locations
    has_fog_locs = !isempty(selected_fog_ranks)
    if has_fog_locs
        fog_locs = findall(ctx.log_location_ranks.locations .∈ [selected_fog_ranks])
        fog_locations!(@view(ctx.fogging_log[tstep, :]), fog_locs, ctx.dhw_t, ctx.fogging)
    end

    return selected_fog_ranks
end