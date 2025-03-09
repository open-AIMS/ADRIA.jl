import CoralBlox: FunctionalGroup

import .decision: DecisionPreference


"""
    SimulationContext

Holds the state during a model run, separating the data flow from processing logic.
"""
mutable struct SimulationContext

    # Could create types/structs for each conceptual area of concern and compose them here
    # but for now we use a flat structure.

    domain::Domain

    # Base factors/settings
    param_set::YAXArray

    # Coral parameters
    corals::DataFrame
    functional_groups::Vector{Vector{FunctionalGroup}}
    fec_params_per_m²::Matrix{Float64}
    linear_extensions::Matrix{Float64}
    bin_edges::Matrix{Float64}
    growth_rate::Matrix{Float64}

    # Environmental scenarios
    dhw_scen::YAXArray  # The specific DHW scenario from domain.dhw_scen
    wave_scen::YAXArray  # The specific wave scenario from domain.wave_scen
    cyclone_mortality_scen::YAXArray  # As above, but for domain.cyclone_mortality_scen

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

    # Bioregion-related fields
    unique_biogroups::Vector{Int64}
    n_biogroups::Int64
    biogroup_masks::BitMatrix
    loc_biogrp_idxs::Vector{Int64}

    # Scale factors
    growth_accel_parameters::Matrix{Float64}
    growth_acc_steepness::Vector{Float64}
    growth_acc_height::Vector{Float64}
    growth_acc_midpoint::Vector{Float64}
    scale_factors::Array{Float64, 3}
    biogrp_lin_ext::Array{Float64, 3}
    biogrp_survival::Array{Float64, 3}
    growth_constraints::Vector{Float64}

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

    # Initialize bioregion-related properties
    initialize_bioregions!(ctx)

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

function initialize_bioregions!(ctx::SimulationContext)
    if !("GROUPED_BIOREGION" in names(ctx.domain.loc_data))
        ctx.domain.loc_data.GROUPED_BIOREGION .= 1
    end

    loc_data = ctx.domain.loc_data

    # Extract unique biogroups
    ctx.unique_biogroups = unique(loc_data.GROUPED_BIOREGION)
    ctx.n_biogroups = length(ctx.unique_biogroups)

    # Create biogroup masks
    ctx.biogroup_masks = falses(ctx.n_locs, ctx.n_biogroups)
    for (idx, biogroup) in enumerate(ctx.unique_biogroups)
        ctx.biogroup_masks[:, idx] .= loc_data.GROUPED_BIOREGION .== biogroup
    end

    # Create biogroup indices for each location
    ctx.loc_biogrp_idxs = [
        findfirst(x -> x == biogrp, ctx.unique_biogroups)
        for biogrp in loc_data.GROUPED_BIOREGION
    ]

    if length(ctx.unique_biogroups) > 1
        # Extract growth acceleration parameters
        growth_acc_names = accel_params_array_to_vec(
            generate_growth_accel_names(ctx.unique_biogroups)
        )
        ctx.growth_accel_parameters = accel_params_vec_to_array(
            ctx.param_set[factors=At(growth_acc_names)], ctx.n_biogroups
        )

        ctx.growth_acc_steepness = ctx.growth_accel_parameters[:, 1]
        ctx.growth_acc_height = ctx.growth_accel_parameters[:, 2]
        ctx.growth_acc_midpoint = ctx.growth_accel_parameters[:, 3]

        # Extract scale factors
        sf_col_names = scale_factor_array_to_vec(
            generate_scale_factor_names(ctx.unique_biogroups)
        )
        ctx.scale_factors = scale_factor_vec_to_array(
            ctx.param_set[factors=At(sf_col_names)],
            ctx.n_groups,
            ctx.n_biogroups,
            2
        )

        # Initialize bioregion-specific extensions and survival rates
        ctx.biogrp_lin_ext = repeat(ctx.linear_extensions, 1, 1, ctx.n_biogroups)
        ctx.biogrp_survival = repeat(ctx.survival_rate, 1, 1, ctx.n_biogroups)
        for i in 1:ctx.n_biogroups
            ctx.biogrp_lin_ext[:, :, i] .*= ctx.scale_factors[:, 1, i]
            ctx.biogrp_survival[:, :, i] .= apply_mortality_scaling(
                ctx.biogrp_survival[:, :, i], ctx.scale_factors[:, 2, i]
            )
        end
    end

    # Preallocate vector for growth constraints
    ctx.growth_constraints = zeros(Float64, ctx.n_locs)
end
