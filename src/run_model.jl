import .decision: DecisionPreference

"""
    SimulationContext

Holds the state during a model run, separating the data flow from processing logic.
"""
mutable struct SimulationContext
    # Domain data
    domain::Domain

    # Scenario factors/settings
    param_set::YAXArray

    # Coral parameters
    corals::DataFrame
    functional_groups::Vector{Vector{FunctionalGroup}}
    fec_params_per_m²::Matrix{Float64}
    linear_extensions::Matrix{Float64}
    bin_edges::Matrix{Float64}
    growth_rate::Matrix{Float64}

    # Environmental scenarios
    dhw_scen::YAXArray
    wave_scen::YAXArray
    cyclone_mortality_scen::YAXArray

    # Configuration parameters
    is_guided::Bool
    is_unguided::Bool
    MCDA_approach::Union{Nothing,Function,DataType}
    plan_horizon::Int64
    seeding::Bool
    apply_seeding::Bool
    apply_fogging::Bool
    apply_shading::Bool
    depth_criteria::BitVector

    # Intervention parameters
    seed_decision_years::Vector{Bool}
    fog_decision_years::Vector{Bool}
    shade_decision_years::Vector{Bool}
    min_iv_locs::Int64
    max_members::Int64
    max_area_to_seed::Float64

    # Simulation constants
    tf::Int64
    n_locs::Int64
    n_groups::Int64
    n_sizes::Int64
    n_group_and_size::Int64
    habitable_locs::BitVector
    habitable_loc_areas::Vector{Float64}
    habitable_loc_areas′::Array{Float64,3}
    habitable_loc_idxs::Vector{Int64}
    habitable_area::Matrix{Float64}

    # Connectivity
    valid_sources::BitVector
    valid_sinks::BitVector

    # Derived parameters
    fogging::Float64
    srm::Float64
    seed_sc::BitMatrix
    max_seeded_area::YAXArray
    a_adapt::Matrix{Float64}

    # Work/intermediate matrices
    c_mean_t_1::Array{Float64,3}
    c_mean_t::Array{Float64,3}
    c_std::Matrix{Float64}

    # Results
    C_cover::Array{Float64,4}  # Timestep, groups, sizes, locations
    log_location_ranks::YAXArray
    shading_log::Array{Float64,2}  # Timestep, locations
    fogging_log::Array{Float64,2}
    seeding_log::Array{Float64,3}  # Timestep, groups, locations
    bleaching_mort::Array{Float64,4}
    current_tstep::Int64
    dhw_tol_mean_log::Array{Float64,3}

    # Caches and temporary variables
    dhw_t::Vector{Float64}
    C_cover_t::Array{Float64,3}
    ΔC_cover_t::Array{Float64,3}
    depth_coeff::Vector{Float64}
    fec_scope::Matrix{Float64}
    fec_all::Array{Float64,3}
    recruitment::Matrix{Float64}
    potential_settlers::Matrix{Float64}
    survival_rate::Array{Float64,3}
    interim_dhw_tol_log::Matrix{Float64}

    # Decision making
    decision_mat::Union{Nothing,YAXArray}
    seed_pref::Union{Nothing,DecisionPreference}
    fog_pref::Union{Nothing,DecisionPreference}

    # Constructor (with defaults for optional params)
    function SimulationContext(
        domain::Domain, param_set::YAXArray,
        functional_groups::Vector{Vector{FunctionalGroup}}
    )
        # Initialize with base properties
        return new(domain, param_set, to_coral_spec(param_set), functional_groups)
    end
end

"""
    initialize_context!(ctx::SimulationContext)

Set up all the context variables needed for a simulation run.
"""
function initialize_context!(ctx::SimulationContext)
    # Set random seed
    rnd_seed_val = floor(Int64, sum(ctx.param_set[Where(x -> x != "RCP")]))
    Random.seed!(rnd_seed_val)

    # Initialize environmental scenarios
    initialize_environmental_scenarios!(ctx)

    # Set up simulation constants
    ctx.tf = size(ctx.dhw_scen, 1)
    ctx.n_locs = n_locations(ctx.domain)
    ctx.n_groups = ctx.domain.coral_details.n_groups
    ctx.n_sizes = ctx.domain.coral_details.n_sizes
    ctx.n_group_and_size = ctx.domain.coral_details.n_group_and_size

    # Initialize habitat information
    initialize_habitat_data!(ctx)

    # Set up intervention parameters
    initialize_intervention_parameters!(ctx)

    # Create results storage
    initialize_result_matrices!(ctx)

    # Empty the old contents of the buffers and add the new blocks
    cover_view = [@view ctx.C_cover[1, :, :, loc] for loc in 1:(ctx.n_locs)]
    ctx.functional_groups = reuse_buffers!.(
        ctx.functional_groups, (cover_view .* loc_k_area(ctx.domain))
    )

    # Initialize coral fecundity parameters
    ctx.fec_params_per_m² = _to_group_size(
        ctx.domain.coral_details, ctx.corals.fecundity
    )

    ctx.linear_extensions = _to_group_size(
        ctx.domain.coral_details, ctx.corals.linear_extension
    )

    ctx.bin_edges = bin_edges()
    ctx.growth_rate = growth_rate(ctx.linear_extensions, bin_widths())

    # Set up derived parameters
    ctx.current_tstep = 1

    return ctx
end

"""
    initialize_environmental_scenarios!(ctx::SimulationContext)

Extract and process environmental scenarios from the domain based on configuration.
"""
function initialize_environmental_scenarios!(ctx::SimulationContext)
    # Extract DHW scenario
    dhw_idx = Int64(ctx.param_set[At("dhw_scenario")])
    if dhw_idx > 0.0
        ctx.dhw_scen = @view(ctx.domain.dhw_scens[:, :, dhw_idx])
    else
        # Run with no DHW disturbances
        ctx.dhw_scen = copy(ctx.domain.dhw_scens[:, :, 1])
        ctx.dhw_scen .= 0.0
    end

    # Extract wave scenario
    wave_idx = Int64(ctx.param_set[At("wave_scenario")])
    if wave_idx > 0.0
        ctx.wave_scen = copy(ctx.domain.wave_scens[:, :, wave_idx])
        ctx.wave_scen .= ctx.wave_scen ./ maximum(ctx.wave_scen)
        replace!(ctx.wave_scen, Inf => 0.0, NaN => 0.0)
    else
        ctx.wave_scen = copy(ctx.domain.wave_scens[:, :, 1])
        ctx.wave_scen .= 0.0
    end

    # Extract cyclone mortality scenario
    cyclone_mortality_idx = Int64(ctx.param_set[At("cyclone_mortality_scenario")])
    if cyclone_mortality_idx > 0.0
        ctx.cyclone_mortality_scen = @view(
            ctx.domain.cyclone_mortality_scens[:, :, :, cyclone_mortality_idx]
        )
    else
        ctx.cyclone_mortality_scen = copy(ctx.domain.cyclone_mortality_scens[:, :, :, 1])
        ctx.cyclone_mortality_scen .= 0.0
    end
end

"""
    initialize_habitat_data!(ctx::SimulationContext)::Nothing

Initialize habitat-related data in the simulation context.
"""
function initialize_habitat_data!(ctx::SimulationContext)::Nothing
    # Set up habitat data
    vec_abs_k = loc_k_area(ctx.domain)
    ctx.habitable_locs = location_k(ctx.domain) .> 0.0
    ctx.habitable_loc_areas = vec_abs_k[ctx.habitable_locs]
    ctx.habitable_loc_areas′ = reshape(
        ctx.habitable_loc_areas, (1, 1, length(ctx.habitable_locs))
    )
    ctx.habitable_loc_idxs = findall(ctx.habitable_locs)
    ctx.habitable_area = Matrix{Float64}(vec_abs_k')

    # Depth criteria for interventions
    ctx.depth_criteria = identify_within_depth_bounds(
        ctx.domain.loc_data.depth_med,
        ctx.param_set[At("depth_min")],
        ctx.param_set[At("depth_offset")]
    )

    # Cache depth coefficients
    ctx.depth_coeff = depth_coefficient.(ctx.domain.loc_data.depth_med)

    return nothing
end

"""
    initialize_intervention_parameters!(ctx::SimulationContext)::Nothing

Set up parameters related to interventions (seeding, fogging, shading).
"""
function initialize_intervention_parameters!(ctx::SimulationContext)::Nothing
    # Determine intervention strategy
    ctx.is_guided = ctx.param_set[At("guided")] > 0
    ctx.is_unguided = ctx.param_set[At("guided")] == 0.0

    if ctx.is_guided
        ctx.MCDA_approach = mcda_methods()[Int64(ctx.param_set[At("guided")])]
        ctx.plan_horizon = Int64(ctx.param_set[At("plan_horizon")])
    else
        ctx.MCDA_approach = nothing
        ctx.plan_horizon = 0
    end

    # Locations to intervene
    ctx.min_iv_locs = ctx.param_set[At("min_iv_locations")]
    ctx.max_members = ctx.param_set[At("cluster_max_member")]

    # Intervention flags
    taxa_names = collect(ctx.param_set.factors[occursin.("N_seed_", ctx.param_set.factors)])
    ctx.seeding = any(ctx.param_set[At(taxa_names)] .> 0.0)
    ctx.apply_seeding = ctx.is_unguided && ctx.seeding
    ctx.fogging = ctx.param_set[At("fogging")]
    ctx.apply_fogging = ctx.is_unguided && (ctx.fogging > 0.0)
    ctx.srm = ctx.param_set[At("SRM")]
    ctx.apply_shading = ctx.srm > 0.0

    # Years to start interventions
    seed_start_year = Int64(ctx.param_set[At("seed_year_start")])
    shade_start_year = Int64(ctx.param_set[At("shade_year_start")])
    fog_start_year = Int64(ctx.param_set[At("fog_year_start")])

    # Intervention durations
    seed_years = Int64(ctx.param_set[At("seed_years")])
    shade_years = Int64(ctx.param_set[At("shade_years")])
    fog_years = Int64(ctx.param_set[At("fog_years")])

    # Decision frequencies
    ctx.seed_decision_years = decision_frequency(
        seed_start_year, ctx.tf, seed_years,
        Int64(ctx.param_set[At("seed_deployment_freq")])
    )
    ctx.fog_decision_years = decision_frequency(
        fog_start_year, ctx.tf, fog_years, Int64(ctx.param_set[At("fog_deployment_freq")])
    )
    ctx.shade_decision_years = decision_frequency(
        shade_start_year, ctx.tf, shade_years,
        Int64(ctx.param_set[At("shade_deployment_freq")])
    )

    # Set up seeding parameters
    initialize_seeding_parameters!(ctx, taxa_names)

    # Initialize decision matrices if using guided approach
    if ctx.is_guided
        initialize_decision_matrices!(ctx)
    end
end

"""
    initialize_seeding_parameters!(ctx::SimulationContext, taxa_names::Vector{String})

Set up parameters specific to coral seeding.
"""
function initialize_seeding_parameters!(ctx::SimulationContext, taxa_names::Vector{String})
    # Define taxa and size class to seed
    taxa_to_seed = [2, 3, 5]
    target_class_id = ctx.corals.class_id .== 1
    ctx.seed_sc = _to_group_size(
        ctx.domain.coral_details,
        (ctx.corals.taxa_id .∈ [taxa_to_seed]) .& target_class_id
    )

    # Extract colony areas and determine approximate seeded area in m^2
    seed_volume = ctx.param_set[At(taxa_names)]
    colony_areas = _to_group_size(
        ctx.domain.coral_details,
        colony_mean_area(ctx.corals.mean_colony_diameter_m)
    )
    ctx.max_seeded_area = colony_areas[ctx.seed_sc] .* seed_volume

    # Set up assisted adaptation values
    ctx.a_adapt = zeros(ctx.n_groups, ctx.n_sizes)
    ctx.a_adapt[ctx.seed_sc] .= ctx.param_set[At("a_adapt")]

    # Enhance from mean of "natural" DHW tolerance
    ctx.a_adapt[ctx.a_adapt .> 0.0] .+= _to_group_size(
        ctx.domain.coral_details, ctx.corals.dist_mean
    )[ctx.a_adapt .> 0.0]

    # Calculate total area to seed
    return ctx.max_area_to_seed = sum(ctx.max_seeded_area)
end

"""
    initialize_decision_matrices!(ctx::SimulationContext)::Nothing

Initialize decision matrices for guided intervention approaches.
"""
function initialize_decision_matrices!(ctx::SimulationContext)::Nothing
    ctx.seed_pref = SeedPreferences(ctx.domain, ctx.param_set)
    ctx.fog_pref = FogPreferences(ctx.domain, ctx.param_set)

    # Create shared decision matrix
    ctx.decision_mat = decision_matrix(
        ctx.domain.loc_ids,
        ctx.seed_pref.names;
        depth=ctx.domain.loc_data.depth_med
    )

    # Remove locations that cannot support corals or are out of depth bounds
    _valid_locs = ctx.habitable_locs .& ctx.depth_criteria

    ctx.decision_mat = ctx.decision_mat[_valid_locs, :]

    return nothing
end

"""
    initialize_result_matrices!(ctx::SimulationContext)::Nothing

Initialize matrices to store simulation results.
"""
function initialize_result_matrices!(ctx::SimulationContext)::Nothing
    n_tf, n_groups, n_sizes, n_locs = ctx.tf, ctx.n_groups, ctx.n_sizes, ctx.n_locs

    # Coral cover relative to available area
    ctx.C_cover = zeros(n_tf, n_groups, n_sizes, n_locs)
    ctx.C_cover[1, :, :, :] .= _reshape_init_cover(
        ctx.domain.init_coral_cover, (n_sizes, n_groups, n_locs)
    )

    # Location ranks for interventions
    ctx.log_location_ranks = ZeroDataCube(;
        T=Float64,
        timesteps=timesteps(ctx.domain),
        locations=ctx.domain.loc_ids,
        intervention=interventions()
    )

    # Intervention logs
    ctx.shading_log = zeros(n_tf, n_locs)
    ctx.fogging_log = zeros(n_tf, n_locs)
    ctx.seeding_log = zeros(n_tf, 3, n_locs)

    # Bleaching mortality log
    ctx.bleaching_mort = zeros(n_tf, n_groups, n_sizes, n_locs)

    # Set up distributions for natural adaptation/heritability
    ctx.c_mean_t_1 = repeat(
        _to_group_size(ctx.domain.coral_details, ctx.corals.dist_mean),
        1, 1, n_locs
    )
    ctx.c_std = _to_group_size(ctx.domain.coral_details, ctx.corals.dist_std)
    ctx.c_mean_t = copy(ctx.c_mean_t_1)

    # Initialize intermediate/tracking arrays
    ctx.C_cover_t = zeros(n_groups, n_sizes, n_locs)
    ctx.survival_rate = ones(n_groups, n_sizes, n_locs)
    ctx.ΔC_cover_t = zeros(n_groups, n_sizes, n_locs)

    # Set up cache for fecundity calculations
    ctx.fec_scope = zeros(n_groups, n_locs)
    ctx.fec_all = zeros(n_groups, n_sizes, n_locs)
    ctx.recruitment = zeros(n_groups, n_locs)
    ctx.dhw_t = zeros(n_locs)

    # Cache for potential settlers
    ctx.potential_settlers = zeros(size(ctx.fec_scope)...)

    return nothing
end

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
    ctx.functional_groups = reuse_buffers!.(
        ctx.functional_groups, (cover_view .* vec_abs_k)
    )

    # Settlers from t-1 grow into observable sizes
    ctx.C_cover_t[:, 1, ctx.habitable_locs] .+= ctx.recruitment

    # Calculate linear extension scale factors
    lin_ext_scale_factors = calculate_linear_extension_factors(ctx)

    # Perform timestep for each location
    coral_details!(ctx, tstep, lin_ext_scale_factors)

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
    coral_details!(ctx::SimulationContext, tstep::Int64, lin_ext_scale_factors::Vector{Float64})

Perform growth calculations for all locations.
"""
function coral_details!(ctx::SimulationContext, tstep::Int64, lin_ext_scale_factors::Vector{Float64})
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
    TP_data = conn ./ sum(conn; dims=1)
    replace!(TP_data, NaN => 0)

    # Calculate fecundity parameters
    fec_params_per_m² = _to_group_size(
        ctx.domain.coral_growth, ctx.corals.fecundity
    )

    # Update settler DHW tolerance
    settler_DHW_tolerance!(
        ctx.c_mean_t_1,
        ctx.c_mean_t,
        vec_abs_k,
        TP_data,  # Pass transition probability matrix, not connectivity
        ctx.recruitment,
        fec_params_per_m²,
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
    shading_intervention!(ctx::SimulationContext, tstep::Int64)

Process the shading intervention phase.
"""
function shading_intervention!(ctx::SimulationContext, tstep::Int64)
    if ctx.apply_shading && ctx.shade_decision_years[tstep]
        # Apply regional cooling effect
        ctx.shading_log[tstep, :] .= ctx.srm

        # Apply reduction in DHW due to SRM
        ctx.dhw_t .= max.(0.0, ctx.dhw_t .- ctx.srm)
    end
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
    seed_locations!(ctx, tstep, selected_seed_ranks, leftover_space_m²)

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
    # Skip if no locations selected
    has_seed_locs = !isempty(selected_seed_ranks)
    if !has_seed_locs
        return nothing
    end

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
