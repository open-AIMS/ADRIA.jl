"""
    simulation_step!(ctx::SimulationContext)::SimulationContext

Perform one time step of the simulation.
"""
function simulation_step!(ctx::SimulationContext)::SimulationContext
    tstep = ctx.current_tstep + 1

    # Run growth stage
    growth_phase!(ctx, tstep)

    # Run adaptation phase
    adaptation_phase!(ctx, tstep)

    # Run reproduction phase
    leftover_space_m² = reproduction_phase!(ctx)

    # Run intervention phases (shading, fogging, seeding)
    intervention_phases!(ctx, tstep, leftover_space_m²)

    # Run disturbance phases (bleaching, cyclones)
    disturbance_phases!(ctx, tstep)

    # Update current time step
    ctx.current_tstep = tstep

    return ctx
end

"""
    growth_phase!(ctx::SimulationContext, tstep::Int64)::Nothing

Perform the coral growth phase of the simulation.
"""
function growth_phase!(ctx::SimulationContext, tstep::Int64)::Nothing
    # Convert cover to absolute values to use within CoralBlox model
    ctx.C_cover_t[:, :, ctx.habitable_locs] .=
        ctx.C_cover[tstep - 1, :, :, ctx.habitable_locs] .* ctx.habitable_loc_areas′

    # IMPORTANT: Update functional groups with current cover state
    cover_view = [@view ctx.C_cover[tstep-1, :, :, loc] for loc in 1:(ctx.n_locs)]
    vec_abs_k = ctx.habitable_area[1, :]
    for i in 1:length(ctx.functional_groups)
        ctx.functional_groups[i] = reuse_buffers!(
            ctx.functional_groups[i],
            cover_view[i] * vec_abs_k[i]
        )
    end

    # Settlers from t-1 grow into observable sizes
    ctx.C_cover_t[:, 1, ctx.habitable_locs] .+= ctx.recruitment

    # Calculate linear extension scale factors
    lin_ext_scale_factors = calculate_linear_extension_factors(ctx)

    # Perform timestep for each location
    coral_growth!(ctx, lin_ext_scale_factors)

    # Convert C_cover_t to relative values after CoralBlox was run
    ctx.C_cover_t[:, :, ctx.habitable_locs] .= (
        @view(ctx.C_cover_t[:, :, ctx.habitable_locs]) ./ ctx.habitable_loc_areas′
    )

    return nothing
end

"""
    calculate_linear_extension_factors(ctx::SimulationContext)

Calculate linear extension scale factors for coral growth.
"""
function calculate_linear_extension_factors(ctx::SimulationContext)

    # Calculate max projected cover
    habitable_max_projected_cover = max_projected_cover(
        ctx.linear_extensions,
        ctx.bin_edges,
        ctx.habitable_loc_areas
    )

    # Calculate scale factors
    lin_ext_scale_factors = linear_extension_scale_factors(
        ctx.C_cover_t[:, :, ctx.habitable_locs],
        ctx.habitable_loc_areas,
        ctx.linear_extensions,
        ctx.bin_edges,
        habitable_max_projected_cover
    )

    # Apply correction for locations with low coral cover
    lin_ext_scale_factors[_loc_coral_cover(ctx.C_cover_t)[ctx.habitable_locs] .< (0.7 .* ctx.habitable_loc_areas)] .=
        1

    return lin_ext_scale_factors
end

"""
    coral_growth!(ctx::SimulationContext, lin_ext_scale_factors::Vector{Float64})

Perform growth calculations for all locations.
"""
function coral_growth!(ctx::SimulationContext, lin_ext_scale_factors::Vector{Float64})
    survival_rate = 1.0 .- _to_group_size(ctx.domain.coral_details, ctx.corals.mb_rate)

    @floop for i in ctx.habitable_loc_idxs
        # Perform timestep
        timestep!(
            ctx.functional_groups[i],
            ctx.recruitment[:, i],
            ctx.linear_extensions .* lin_ext_scale_factors[i],
            survival_rate
        )

        # Write to the cover matrix
        coral_cover(ctx.functional_groups[i], @view(ctx.C_cover_t[:, :, i]))
    end

    # Check if size classes are inappropriately out-growing habitable area
    @assert (
        sum(
            _loc_coral_cover(ctx.C_cover_t)[ctx.habitable_locs] .> ctx.habitable_loc_areas
        ) == 0
    ) "Cover outgrowing habitable area"
end

"""
    adaptation_phase!(ctx::SimulationContext, tstep::Int64)

Perform the natural adaptation phase of the simulation.
"""
function adaptation_phase!(ctx::SimulationContext, tstep::Int64)
    # Process natural adaptation
    adjust_DHW_distribution!(
        @view(ctx.C_cover_t[:, :, :]), ctx.c_mean_t,
        ctx.growth_rate  # THIS IS AN ISSUE, IT MUST COME FROM CORALBLOX.
    )

    # Set values for t to t-1 (for next timestep)
    ctx.c_mean_t_1 .= ctx.c_mean_t

    # Log DHW tolerances if in debug mode
    if parse(Bool, ENV["ADRIA_DEBUG"]) == true
        ctx.dhw_tol_mean_log[tstep, :, :] .= reshape(
            mean.(ctx.c_mean_t), size(ctx.interim_dhw_tol_log)
        )
    end
end

"""
    reproduction_phase!(ctx::SimulationContext)::Vector{Float64}

Perform the coral reproduction and settlement phase.
"""
function reproduction_phase!(ctx::SimulationContext)::Vector{Float64}
    # Calculate fecundity scope
    fecundity_scope!(
        ctx.fec_scope, ctx.fec_all, ctx.fec_params_per_m², ctx.C_cover_t, ctx.habitable_area
    )

    # Calculate leftover space
    loc_coral_cover = _loc_coral_cover(ctx.C_cover_t)
    vec_abs_k = loc_k_area(ctx.domain)
    leftover_space_m² = relative_leftover_space(loc_coral_cover) .* vec_abs_k

    # Reset matrices
    ctx.potential_settlers .= 0.0
    ctx.recruitment .= 0.0

    # Calculate recruitment
    recruitment_phase!(ctx, leftover_space_m², loc_coral_cover, vec_abs_k)

    # Update leftover space after recruitment
    leftover_space_m² .-= dropdims(sum(ctx.recruitment; dims=1); dims=1) .* vec_abs_k

    return leftover_space_m²
end

"""
    recruitment_phase!(
        ctx::SimulationContext, leftover_space_m²::T, loc_coral_cover::T, vec_abs_k::T
    ) where {T<:Vector{Float64}}

Calculate coral recruitment and update tolerances.
"""
function recruitment_phase!(
    ctx::SimulationContext, leftover_space_m²::T, loc_coral_cover::T, vec_abs_k::T
)::Nothing where {T<:Vector{Float64}}
    # Set up parameters
    sim_params = ctx.domain.sim_constants
    conn = ctx.domain.conn
    habitable_area = ctx.habitable_area

    # Preallocate cache for source/sink locations
    ctx.valid_sources = falses(size(conn, 2))
    ctx.valid_sinks = falses(size(conn, 1))

    # Calculate basal area per settler
    basal_area_per_settler = colony_mean_area(
        ctx.corals.mean_colony_diameter_m[ctx.corals.class_id .== 1]
    )

    # Calculate recruitment
    @views ctx.recruitment[:, ctx.habitable_locs] .=
        settler_cover(
            ctx.fec_scope,
            conn,
            leftover_space_m²,
            sim_params.max_settler_density,
            sim_params.max_larval_density,
            basal_area_per_settler,
            ctx.potential_settlers,
            ctx.valid_sources,
            ctx.valid_sinks
        )[
            :, ctx.habitable_locs
        ] ./ habitable_area[:, ctx.habitable_locs]

    # Calculate transition probability matrix
    TP_data::YAXArray{Float64,2} = conn ./ sum(conn; dims=1)
    replace!(TP_data, NaN => 0.0)

    # Update settler DHW tolerance
    settler_DHW_tolerance!(
        ctx.c_mean_t_1,
        ctx.c_mean_t,
        vec_abs_k,
        TP_data,  # Pass transition probability matrix, not connectivity
        ctx.recruitment,
        ctx.fec_params_per_m²,
        ctx.param_set[At("heritability")]
    )

    return nothing
end

"""
    intervention_phases!(ctx::SimulationContext, tstep::Int64, leftover_space_m²::Vector{Float64})::Vector

Process all intervention phases (shading, fogging, seeding).
"""
function intervention_phases!(
    ctx::SimulationContext, tstep::Int64, leftover_space_m²::Vector{Float64}
)
    # Process DHW for this timestep
    ctx.dhw_t .= ctx.dhw_scen[tstep, :]

    # Process shading intervention
    shading_intervention!(ctx, tstep)

    # Process fogging intervention
    fogging_intervention!(ctx, tstep, leftover_space_m²)

    # Process seeding intervention
    return seeding_intervention!(ctx, tstep, leftover_space_m²)
end

"""
    disturbance_phases!(ctx::SimulationContext, tstep::Int64)::Nothing

Process all disturbance phases (bleaching, cyclones).
"""
function disturbance_phases!(ctx::SimulationContext, tstep::Int64)::Nothing
    # Copy current cover for calculating survival rates
    ctx.ΔC_cover_t .= copy(ctx.C_cover_t)

    # Apply cyclone mortality
    cyclone_mortality!(
        ctx.C_cover_t,
        ctx.cyclone_mortality_scen[tstep, :, :]'
    )

    # Apply bleaching mortality
    bleaching_mortality!(
        ctx.C_cover_t,
        ctx.dhw_t,
        ctx.depth_coeff,
        ctx.c_std,
        ctx.c_mean_t_1,
        ctx.c_mean_t,
        @view(ctx.bleaching_mort[(tstep - 1):tstep, :, :, :])
    )

    # Determine survival and apply to functional groups
    apply_survival_rates!(ctx, tstep)

    # Update C_cover after all disturbances are applied
    ctx.C_cover[tstep, :, :, :] .= ctx.C_cover_t

    return nothing
end

"""
    apply_survival_rates!(ctx::SimulationContext, tstep::Int64)::Nothing

Calculate survival rates from disturbances and apply to functional groups.
"""
function apply_survival_rates!(ctx::SimulationContext, tstep::Int64)::Nothing
    # Calculate survival rates
    ctx.ΔC_cover_t[ctx.ΔC_cover_t .== 0.0] .= 1.0
    ctx.survival_rate .= ctx.C_cover_t ./ ctx.ΔC_cover_t
    @assert sum(ctx.survival_rate .> 1) == 0 "Survival rate should be <= 1"

    # Apply mortality to functional groups
    survival_rate_slices = [
        @view ctx.survival_rate[:, :, loc] for loc in 1:(ctx.n_locs)
    ]
    apply_mortality!.(ctx.functional_groups, survival_rate_slices)

    # Apply mortality to recruitment
    ctx.recruitment .*= reshape(
        ctx.survival_rate[:, 1, :], (ctx.n_groups, ctx.n_locs)
    )

    ctx.recruitment .*= ctx.habitable_area

    return nothing
end

"""
    run_model(domain::Domain, param_set::Union{DataFrameRow,YAXArray})::NamedTuple
    run_model(domain::Domain, param_set::DataFrameRow, functional_groups::Vector{Vector{FunctionalGroup}})::NamedTuple
    run_model(domain::Domain, param_set::YAXArray, functional_groups::Vector{Vector{FunctionalGroup}})::NamedTuple

Core scenario running function.
"""
function run_model(domain::Domain, param_set::Union{DataFrameRow,YAXArray})::NamedTuple
    n_locs = domain.coral_details.n_locs
    n_sizes = domain.coral_details.n_sizes
    n_groups = domain.coral_details.n_groups
    _bin_edges = bin_edges()  # TODO: Move to CoralDetails
    functional_groups = Vector{FunctionalGroup}[
        FunctionalGroup.(
            eachrow(_bin_edges[:, 1:(end - 1)]),
            eachrow(_bin_edges[:, 2:end]),
            eachrow(zeros(n_groups, n_sizes))
        ) for _ in 1:n_locs
    ]
    return run_model(domain, param_set, functional_groups)
end

function run_model(
    domain::Domain,
    param_set::DataFrameRow,
    functional_groups::Vector{Vector{FunctionalGroup}}
)::NamedTuple
    setup()
    ps = DataCube(Vector(param_set); factors=names(param_set))
    return run_model(domain, ps, functional_groups)
end

function run_model(
    domain::Domain,
    param_set::YAXArray,
    functional_groups::Vector{Vector{FunctionalGroup}}
)::NamedTuple
    # Create and initialize the simulation context
    ctx = SimulationContext(domain, param_set, functional_groups)
    initialize_context!(ctx)

    # Run simulation for all time steps
    for _ in 2:(ctx.tf)
        simulation_step!(ctx)
    end

    # Process and return results
    return process_results(ctx)
end

"""
    process_results(ctx::SimulationContext)::NamedTuple

Process and return final simulation results.
"""
function process_results(ctx::SimulationContext)::NamedTuple
    # Final reshape of cover data
    raw = reshape(
        permutedims(ctx.C_cover, (1, 3, 2, 4)),
        (ctx.tf, ctx.n_group_and_size, ctx.n_locs)
    )

    # Process DHW tolerance log if in debug mode
    if parse(Bool, ENV["ADRIA_DEBUG"]) == true
        coral_dhw_log = DataCube(
            ctx.dhw_tol_mean_log;
            timesteps=timesteps(ctx.domain),
            species=ctx.corals.coral_id,
            sites=1:(ctx.n_locs)
        )
    else
        coral_dhw_log = false
    end

    # Return results
    return (
        raw=raw,
        seed_log=ctx.seeding_log,
        fog_log=ctx.fogging_log,
        shade_log=ctx.shading_log,
        site_ranks=ctx.log_location_ranks,
        bleaching_mortality=ctx.bleaching_mort,
        coral_dhw_log=coral_dhw_log
    )
end
