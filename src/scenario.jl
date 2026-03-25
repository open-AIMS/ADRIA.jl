"""Scenario running functions"""

using CoralBlox

import CoralBlox:
    SizeClass,
    FunctionalGroup,
    reuse_buffers!,
    apply_mortality!,
    timestep!,
    coral_cover,
    max_projected_cover,
    linear_extension_scale_factors

using .metrics:
    relative_cover,
    relative_loc_taxa_cover,
    total_absolute_cover,
    absolute_shelter_volume,
    relative_shelter_volume,
    relative_juveniles,
    relative_taxa_cover,
    juvenile_indicator,
    coral_evenness

using SparseArrays
using DataStructures: CircularBuffer

using .decision

"""
    setup_cache(domain::Domain)::NamedTuple

Establish tuple of matrices/vectors for use as reusable data stores to avoid repeated memory allocations.
"""
function setup_cache(domain::Domain)::NamedTuple

    # Simulation constants
    n_locs::Int64 = domain.coral_growth.n_locs
    n_sizes::Int64 = domain.coral_growth.n_sizes
    n_groups::Int64 = domain.coral_growth.n_groups
    tf = length(timesteps(domain))

    cache = (
        # sf=zeros(n_groups, n_locs),  # stressed fecundity, commented out as it is disabled
        fec_scope=zeros(n_groups, n_locs),  # fecundity scope
        recruitment=zeros(n_groups, n_locs),  # coral recruitment
        dhw_step=zeros(n_locs),  # DHW for each time step
        C_cover_t=zeros(n_groups, n_sizes, n_locs),  # Cover for previous timestep
        depth_coeff=zeros(n_locs),  # store for depth coefficient
        loc_area=Matrix{Float64}(loc_area(domain)'),  # area of locations
        habitable_area=Matrix{Float64}(loc_k_area(domain)'),  # location carrying capacity
        wave_damage=zeros(tf, n_sizes * n_groups, n_locs),  # damage coefficient for each size class
        dhw_tol_mean_log=zeros(tf, n_sizes * n_groups, n_locs)  # tmp log for mean dhw tolerances
    )

    return cache
end

"""
    _reshape_init_cover(data::AbstractMatrix{<:Union{Float32, Float64}})

Reshape initial coral cover of shape [n_groups * n_sizes ⋅ n_locations] to shape
[groups ⋅ sizes ⋅ locations]
"""
function _reshape_init_cover(
    C_cover::AbstractMatrix{T},
    dims::NTuple{3,Int64}
)::Array{T} where {T<:Union{Float32,Float64}}
    return permutedims(
        reshape(
            C_cover,
            dims
        ),
        (2, 1, 3)
    )
end

"""
    apply_survival_scaling(survival_rate::AbstractMatrix, scaling_param::AbstractVector)::AbstractMatrix
"""
function apply_survival_scaling(
    survival_rate::AbstractMatrix,
    scaling_param::AbstractVector
)::AbstractMatrix
    # If `scaling_param` is a very low negative number this function will return
    # negative numbers for lower values of `survival_rate`. For instance, when
    # `scaling_param = -1` that will happen whenever `survival_rate < 1/3`.
    # In practice, the lowest `survival_rate` value we have is 0.716, so that issue never
    # happens.
    return survival_rate .+ 0.5 .* (1 .- survival_rate) .* scaling_param
end

"""
    _to_group_size(growth_spec::CoralGrowth, data::AbstractVector{<:Union{Float32, Float64, Tuple}})::Matrix{<:Union{Float32, Float64}}

Reshape vector to shape [functional_groups ⋅ sizes]
"""
function _to_group_size(
    growth_spec::CoralGrowth, data::AbstractVector{T}
)::Matrix{T} where {T<:Union{Float32,Float64,Bool}}
    # Data is reshaped to size ⋅ groups then transposed to maintain expected order
    return Matrix(reshape(data, (growth_spec.n_sizes, growth_spec.n_groups))')
end

"""
    run_scenarios(dom::Domain, scens::DataFrame, RCP::String; show_progress=true, remove_workers=true)
    run_scenarios(dom::Domain, scens::DataFrame, RCP::Vector{String}; show_progress=true, remove_workers=true)

Run scenarios defined by the parameter table storing results to disk.
Scenarios are run in parallel where the number of scenarios > 256 for base Domain data
packages and > 4 for GBR-scale datasets.

# Notes
- Returned `Domain` holds scenario invoke time used as unique result set identifier.

# Examples
```julia-repl
...
julia> rs_45 = ADRIA.run_scenarios(dom, scens, "45")
julia> rs_45_60 = ADRIA.run_scenarios(dom, scens, ["45", "60"])
```

# Arguments
- `dom` : Domain, to run scenarios with
- `scens` : DataFrame of scenarios to run
- `RCP` : ID or list of RCP(s) to run scenarios under.
- `show_progress` : Display progress
- `remove_workers` : If running in parallel, removes workers after completion

# Returns
ResultSet
"""
function run_scenarios(
    dom::Domain, scens::DataFrame, RCP::String; show_progress=true, remove_workers=true
)::ResultSet
    return run_scenarios(dom, scens, [RCP]; show_progress, remove_workers)
end
function run_scenarios(
    dom::Domain,
    scens::DataFrame,
    RCP::Vector{String};
    show_progress=true,
    remove_workers=true
)::ResultSet
    # Initialize ADRIA configuration options
    setup()

    # Sort RCPs so the dataframe order match the output filepath
    RCP = sort(RCP)

    @info "Running $(nrow(scens)) scenarios over $(length(RCP)) RCPs: $RCP"

    # Cross product between rcps and param_df to have every row of param_df for each rcp
    rcps_df = DataFrame(; RCP=parse.(Int64, RCP))
    scenarios_df = crossjoin(scens, rcps_df)
    sort!(scenarios_df, :RCP)

    @info "Setting up Result Set"
    dom, data_store = ADRIA.setup_result_store!(dom, scenarios_df)

    # Convert DataFrame to named matrix for faster iteration
    scenarios_matrix::YAXArray = DataCube(
        Matrix(scenarios_df); scenarios=1:nrow(scenarios_df), factors=names(scenarios_df)
    )

    # The idea to initialise the functional groups here is so that each scenario run can
    # reuse the memory. After a few scenarios runs there will be very few new
    # allocations as buffers grow larger and are not exceeded
    #
    # The problem is that even if a subsequent run only has to reallocate once, the
    # buffer is so large it takes very long regardless.
    #
    # This is also forced me to added an argument to run_model which breaks
    # encapsulation, so its quite ugly
    n_locs::Int64 = dom.coral_growth.n_locs
    n_sizes::Int64 = dom.coral_growth.n_sizes
    n_groups::Int64 = dom.coral_growth.n_groups
    _bin_edges::Matrix{Float64} = bin_edges(; unit=:m)
    functional_groups = [
        FunctionalGroup.(
            eachrow(_bin_edges[:, 1:(end - 1)]),
            eachrow(_bin_edges[:, 2:end]),
            eachrow(zeros(n_groups, n_sizes))
        ) for _ in 1:n_locs
    ]

    para_threshold::Int64 =
        ((typeof(dom) == RMEDomain) || (typeof(dom) == ReefModDomain)) ? 8 : 256
    active_cores::Int64 = parse(Int64, ENV["ADRIA_NUM_CORES"])
    parallel::Bool =
        (parse(Bool, ENV["ADRIA_DEBUG"]) == false) && (active_cores > 1) &&
        (nrow(scens) >= para_threshold)
    if parallel && nworkers() == 1
        @info "Setting up parallel processing..."
        spinup_time = @elapsed begin
            _setup_workers()

            # Load ADRIA on workers and define helper function
            # Note: Workers do not share the same cache in parallel workloads.
            #       Previously, cache would be reserialized so each worker has access to
            #       a separate cache.
            #       Using CachingPool() resolves the the repeated reserialization but it
            #       seems each worker was then attempting to use the same cache, causing the
            #       Julia kernel to crash in multi-processing contexts.
            #       Getting each worker to create its own cache reduces serialization time
            #       (at the cost of increased run time) but resolves the kernel crash issue.
            @sync @async @everywhere @eval begin
                using ADRIA
                func = (dfx) -> run_scenario(dfx..., functional_groups, data_store)
            end
        end

        @info "Time taken to spin up workers: $(round(spinup_time; digits=2)) seconds"
    end

    if parallel
        # Define local helper
        func = (dfx) -> run_scenario(dfx..., functional_groups, data_store)

        try
            for rcp in RCP
                run_msg = "Running $(nrow(scens)) scenarios for RCP $rcp"

                # Switch RCPs so correct data is loaded
                dom = switch_RCPs!(dom, rcp)
                target_rows = findall(
                    scenarios_matrix[factors=At("RCP")] .== parse(Float64, rcp)
                )
                scen_args = _scenario_args(dom, scenarios_matrix, rcp, length(target_rows))

                if show_progress
                    @showprogress desc = run_msg dt = 4 pmap(
                        func, CachingPool(workers()), scen_args
                    )
                else
                    pmap(func, CachingPool(workers()), scen_args)
                end
            end
        catch err
            # Remove hanging workers if any error occurs
            _remove_workers()
            rethrow(err)
        end
    else
        # Define local helper
        func = dfx -> run_scenario(dfx..., functional_groups, data_store)

        for rcp in RCP
            run_msg = "Running $(nrow(scens)) scenarios for RCP $rcp"

            # Switch RCPs so correct data is loaded
            dom = switch_RCPs!(dom, rcp)
            scen_args = _scenario_args(
                dom, scenarios_matrix, rcp, size(scenarios_matrix, 1)
            )

            if show_progress
                @showprogress desc = run_msg dt = 4 map(func, scen_args)
            else
                map(func, scen_args)
            end
        end
    end

    if parallel && remove_workers
        _remove_workers()
    end

    return load_results(_result_location(dom, RCP))
end

function growth_acceleration(
    height::Float64, midpoint::Float64, steepness::Float64, available_space::Float64
)
    return height / (1 + exp(-steepness * (available_space - midpoint))) + 1.0
end

function _scenario_args(dom, scenarios_matrix::YAXArray, rcp::String, n::Int)
    target_rows = findall(
        collect(scenarios_matrix[factors=At("RCP")]) .== parse(Float64, rcp)
    )
    rep_doms = Iterators.repeated(dom, n)
    return zip(
        rep_doms, target_rows, eachrow(scenarios_matrix[target_rows, :])
    )
end

"""
    run_scenario(domain::Domain, idx::Int64, scenario::Union{AbstractVector,DataFrameRow}, functional_groups::Vector{Vector{FunctionalGroup}}, data_store::NamedTuple)::Nothing
    run_scenario(domain::Domain, scenario::Union{AbstractVector,DataFrameRow})::NamedTuple
    run_scenario(domain::Domain, scenario::Union{AbstractVector,DataFrameRow}, RCP::String)::NamedTuple

Run individual scenarios for a given domain, saving results to a Zarr data store.
Results are stored in Zarr format at a pre-configured location.
Sets up a new `cache` if not provided.

# Arguments
- `domain` : Simulation domain (may be modified via `switch_RCPs!`).
- `idx` : Scenario index, to store results into `data_store`.
- `scenario` : Parameter row describing the scenario.
- `functional_groups` : Preallocated functional group buffers.
- `data_store` : Pre-opened store with arrays to write results into.

# Returns
Nothing
"""
function run_scenario(
    domain::Domain,
    idx::Int64,
    scenario::Union{AbstractVector,DataFrameRow},
    functional_groups::Vector{Vector{FunctionalGroup}},
    data_store::NamedTuple
)::Nothing
    if domain.RCP == ""
        local rcp
        try
            rcp = string(Int64(scenario[At("RCP")]))
        catch err
            if !(err isa MethodError)
                rethrow(err)
            end

            rcp = scenario.RCP  # Extract from dataframe
        end

        domain = switch_RCPs!(domain, rcp)
    end

    result_set = run_model(domain, scenario, functional_groups)

    # Capture results to disk
    # Set values below threshold to 0 to save space
    threshold = parse(Float32, ENV["ADRIA_THRESHOLD"])

    # rs_raw has dimensions [timesteps ⋅ group ⋅ sizes ⋅ locations]
    rs_raw::Array{Float64} = result_set.raw
    vals = relative_cover(rs_raw)
    vals[vals .< threshold] .= 0.0
    data_store.relative_cover[:, :, idx] .= vals

    vals = absolute_shelter_volume(rs_raw, loc_k_area(domain), scenario)
    vals[vals .< threshold] .= 0.0
    data_store.absolute_shelter_volume[:, :, idx] .= vals
    vals = relative_shelter_volume(rs_raw, loc_k_area(domain), scenario)
    vals[vals .< threshold] .= 0.0
    data_store.relative_shelter_volume[:, :, idx] .= vals

    coral_spec::DataFrame = to_coral_spec(scenario)
    vals = relative_juveniles(rs_raw)
    vals[vals .< threshold] .= 0.0
    data_store.relative_juveniles[:, :, idx] .= vals

    vals = juvenile_indicator(rs_raw, coral_spec, loc_k_area(domain))
    vals[vals .< threshold] .= 0.0
    data_store.juvenile_indicator[:, :, idx] .= vals

    vals = relative_taxa_cover(rs_raw, loc_k_area(domain))
    vals[vals .< threshold] .= 0.0
    data_store.relative_taxa_cover[:, :, idx] .= vals

    vals = relative_loc_taxa_cover(rs_raw)

    vals = coral_evenness(vals.data)
    vals[vals .< threshold] .= 0.0
    data_store.coral_evenness[:, :, idx] .= vals

    # Store logs
    c_dim = Base.ndims(result_set.raw) + 1
    log_stores = (:site_ranks, :mc_log, :seed_log, :fog_log, :shade_log, :coral_dhw_log)
    for k in log_stores
        if k == :seed_log || k == :site_ranks
            concat_dim = c_dim
        else
            concat_dim = c_dim - 1
        end

        vals = getfield(result_set, k)

        try
            vals[vals .< threshold] .= Float32(0.0)
        catch err
            err isa MethodError ? nothing : rethrow(err)
        end

        if k == :seed_log || k == :mc_log
            getfield(data_store, k)[:, :, :, idx] .= vals
        elseif k == :site_ranks
            if !isnothing(data_store.site_ranks)
                # Squash site ranks down to average rankings over environmental repeats
                data_store.site_ranks[:, :, :, idx] .= vals
            end
        elseif k == :coral_dhw_log
            # Only log coral DHW tolerances if in debug mode
            if parse(Bool, ENV["ADRIA_DEBUG"]) == true
                getfield(data_store, k)[:, :, :, idx] .= vals
            end
        else
            getfield(data_store, k)[:, :, idx] .= vals
        end
    end

    return nothing
end
function run_scenario(
    domain::Domain, scenario::Union{AbstractVector,DataFrameRow}
)::NamedTuple
    return run_model(domain, scenario)
end
function run_scenario(
    domain::Domain, scenario::Union{AbstractVector,DataFrameRow}, RCP::String
)::NamedTuple
    domain = switch_RCPs!(domain, RCP)
    return run_scenario(domain, scenario)
end

"""
    run_model(domain::Domain, param_set::Union{DataFrameRow,YAXArray})::NamedTuple
    run_model(domain::Domain, param_set::DataFrameRow, functional_groups::Vector{Vector{FunctionalGroup}})::NamedTuple
    run_model(domain::Domain, param_set::YAXArray, functional_groups::Vector{Vector{FunctionalGroup}})::NamedTuple

Core scenario running function. When called with only `domain` and `param_set` arguments,
`ADRIA.setup()` must be run beforehand.

# Returns
NamedTuple of collated results
- `raw` : Array, Coral cover relative to k area
- `seed_log` : Array, Log of seeded locations
- `fog_log` : Array, Log of fogged locations
- `shade_log` : Array, Log of shaded locations
- `site_ranks` : Array, Log of location rankings
- `bleaching_mortality` : Array, Log of mortalities caused by bleaching
- `coral_dhw_log` : Array, Log of DHW tolerances / adaptation over time (only logged in debug mode)
"""
function run_model(
    domain::Domain,
    param_set::Union{DataFrameRow,YAXArray}
)
    n_locs::Int64 = domain.coral_growth.n_locs
    n_sizes::Int64 = domain.coral_growth.n_sizes
    n_groups::Int64 = domain.coral_growth.n_groups

    edges = bin_edges()
    bin_start = edges[:, 1:(end - 1)]
    bin_end = edges[:, 2:end]
    functional_groups::Vector{Vector{FunctionalGroup}} = Vector(undef, n_locs)
    group_info = zeros(n_groups, n_sizes)
    for i in 1:n_locs
        fg::Vector{FunctionalGroup} = Vector(undef, n_groups)
        for g in 1:n_groups
            fg[g] = FunctionalGroup(
                bin_start[g, :],
                bin_end[g, :],
                group_info[g, :]
            )
        end

        functional_groups[i] = fg
    end

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
    corals = to_coral_spec(param_set)
    cache = setup_cache(domain)

    factor_names::Vector{String} = collect(param_set.factors.val)

    # Set random seed using intervention values
    # TODO: More robust way of getting intervention/criteria values
    rnd_seed_val::Int64 = floor(Int64, sum(param_set[Where(x -> x != "RCP")]))  # select everything except RCP
    Random.seed!(rnd_seed_val)

    # Extract environmental data
    dhw_idx::Int64 = Int64(param_set[At("dhw_scenario")])
    if dhw_idx > 0.0
        dhw_scen = @view(domain.dhw_scens[:, :, dhw_idx])
    else
        # Run with no DHW disturbances
        dhw_scen = copy(domain.dhw_scens[:, :, 1])
        dhw_scen .= 0.0
    end

    wave_idx::Int64 = Int64(param_set[At("wave_scenario")])
    if wave_idx > 0.0
        # TODO: Better conversion of Ub to wave mortality
        #       Currently scaling significant wave height by its max to non-dimensionalize values
        wave_scen = copy(domain.wave_scens[:, :, wave_idx])
        wave_scen .= wave_scen ./ maximum(wave_scen)
        replace!(wave_scen, Inf => 0.0, NaN => 0.0)
    else
        # Run with no wave disturbances
        wave_scen = copy(domain.wave_scens[:, :, 1])
        wave_scen .= 0.0
    end

    cyclone_mortality_idx::Int64 = Int64(param_set[At("cyclone_mortality_scenario")])
    if cyclone_mortality_idx > 0.0
        cyclone_mortality_scen = @view(
            domain.cyclone_mortality_scens[:, :, :, cyclone_mortality_idx]
        )
    else
        # Run with no cyclone disturbances
        cyclone_mortality_scen = copy(domain.cyclone_mortality_scens[:, :, :, 1])
        cyclone_mortality_scen .= 0.0
    end
    # Environment variables are stored as strings, so convert to bool for use
    in_debug_mode = parse(Bool, get(ENV, "ADRIA_DEBUG", "false")) == true

    # Sim constants
    sim_params = domain.sim_constants
    tf::Int64 = size(dhw_scen, 1)
    n_locs::Int64 = domain.coral_growth.n_locs
    n_groups::Int64 = domain.coral_growth.n_groups
    n_sizes::Int64 = domain.coral_growth.n_sizes

    # Initialize cover loss tracking for reactive strategies
    max_lookback = Int64(param_set[At("reactive_response_delay")])
    if max_lookback < 1
        recent_cover_losses = CircularBuffer{Vector{Float64}}(1)
    else
        recent_cover_losses = CircularBuffer{Vector{Float64}}(max_lookback)
    end
    push!(recent_cover_losses, zeros(n_locs))

    # Locations to intervene
    min_iv_locs::Int64 = param_set[At("min_iv_locations")]

    fogging::Float64 = param_set[At("fogging")]  # proportion of bleaching mortality reduction through fogging
    srm::Float64 = param_set[At("SRM")]  # DHW equivalents reduced by some shading mechanism
    shade_years::Int64 = param_set[At("shade_years")]  # number of years to shade

    # Years to start seeding/shading/fogging
    shade_start_year::Int64 = param_set[At("shade_year_start")]

    colony_areas = _to_group_size(
        domain.coral_growth, colony_mean_area(corals.mean_colony_diameter_m)
    )

    habitable_areas::Matrix{Float64} = cache.habitable_area
    fecundity_per_m²::Matrix{Float64} = _to_group_size(
        domain.coral_growth, corals.fecundity
    ) # number of larvae produced per m²

    # Caches
    conn = sparse(domain.conn.data)

    # Determine contribution of each source to a sink location
    # i.e., columns should sum to 1!
    TP_data = conn ./ sum(conn; dims=1)
    replace!(nonzeros(TP_data), NaN => 0.0)
    dropzeros!(TP_data)

    # sf = cache.sf  # unused as it is currently deactivated
    fec_scope = cache.fec_scope
    recruitment = cache.recruitment
    dhw_t = cache.dhw_step
    C_cover_t = cache.C_cover_t
    depth_coeff = cache.depth_coeff

    # Used to distribute moving corals settlers
    prop_fecundity = copy(cache.fec_scope)

    loc_data = domain.loc_data
    depth_coeff .= depth_coefficient.(loc_data.depth_med)

    # Coral cover relative to available area (i.e., 1.0 == site is filled to max capacity)
    C_cover::Array{Float64,4} = zeros(tf, n_groups, n_sizes, n_locs)
    C_cover[1, :, :, :] .= _reshape_init_cover(
        domain.init_coral_cover, (n_sizes, n_groups, n_locs)
    )

    # Locations that can support corals
    vec_abs_k = loc_k_area(domain)
    habitable_locs::BitVector = location_k(domain) .> 0.0
    habitable_loc_areas = vec_abs_k[habitable_locs]
    habitable_loc_areas′ = reshape(habitable_loc_areas, (1, 1, length(habitable_locs)))
    n_habitable_locs::Int64 = length(habitable_locs)

    # Avoid placing importance on sites that were not considered
    # Lower values are higher importance/ranks.
    # Values of n_locs+1 indicate locations that were not considered in rankings.
    log_location_ranks = ZeroDataCube(;     # log seeding/fogging ranks
        T=Float64,
        timesteps=1:tf,
        locations=domain.loc_ids,
        intervention=interventions()
    )

    Yshade = spzeros(tf, n_locs)
    Yfog = spzeros(tf, n_locs)
    Yseed = zeros(tf, n_groups, n_locs)  # 3 = the number of seeded coral types
    Ymc = zeros(tf, n_groups, n_locs)

    # Prep scenario-specific flags/values
    # Intervention strategy: < 0 is no intervention, 0 is random location selection, > 0 is guided
    is_guided = param_set[At("guided")] > 0
    if is_guided
        MCDA_approach = mcda_methods()[Int64(param_set[At("guided")])]
    end

    # Number of time steps in environmental layers to look ahead when making decisions
    plan_horizon::Int64 = Int64(param_set[At("plan_horizon")])

    # Decisions should place more weight on environmental conditions
    # closer to the decision point
    α = 0.99
    decay = α .^ (1:(plan_horizon + 1)) .^ 2

    # Years at which intervention locations are re-evaluated and deployed
    # seed_decision_years = decision_frequency(
    #     seed_start_year, tf, seed_years, param_set[At("seed_deployment_freq")]
    # )
    # fog_decision_years = decision_frequency(
    #     fog_start_year, tf, fog_years, param_set[At("fog_deployment_freq")]
    # )
    last_seed_deployment = zeros(Int64, n_locs)
    last_fog_deployment = zeros(Int64, n_locs)
    last_mc_deployment = zeros(Int64, n_locs)
    shade_decision_years = decision_frequency(
        shade_start_year, tf, shade_years, param_set[At("shade_deployment_freq")]
    )

    # Define taxa and size class to seed, and identify their factor names
    # TODO: Seed 1-year old corals!!! If this is the 1st size class, that's fine but needs
    # to be confirmed with ecoRRAP
    _seed_size_groups::BitMatrix = seed_size_groups(n_groups, n_sizes)

    # Set up assisted adaptation values
    a_adapt::Vector{Float64} = fill(param_set[At("a_adapt")], n_groups)

    seed_volume = param_set[At(factor_names[contains.(factor_names, "N_seed")])]
    seeding_devices_per_m2::Float64 = param_set[At("seeding_devices_per_m2")]

    is_unguided = param_set[At("guided")] == 0.0
    is_seeding = any(seed_volume .> 0)

    # Flag indicating whether to seed or not to seed when unguided
    unguided_seeding = is_unguided && is_seeding

    # Flag indicating whether to fog or not fog
    is_fogging = fogging > 0.0
    unguided_fogging = is_unguided && is_fogging

    # Moving corals flag
    n_mc_settlers = param_set[At("N_mc_settlers")]
    is_mc = n_mc_settlers > 0.0
    unguided_moving_corals = is_unguided && is_mc

    # Flag indicating whether to apply shading
    apply_shading = srm > 0.0

    depth_criteria = identify_within_depth_bounds(
        loc_data.depth_med, param_set[At("depth_min")], param_set[At("depth_offset")]
    )

    if is_guided
        seed_pref, seed_decision_mat, seed_strategy = setup_guided_intervention(
            domain, param_set, depth_criteria, SeedPreferences,
            domain.seed_target_locations, is_seeding, build_seed_strategy
        )
        fog_pref, fog_decision_mat, fog_strategy = setup_guided_intervention(
            domain, param_set, depth_criteria, FogPreferences, domain.fog_target_locations,
            is_fogging, build_fog_strategy
        )
        mc_pref, mc_decision_mat, mc_strategy = setup_guided_intervention(
            domain, param_set, depth_criteria, MCPreferences, domain.mc_target_locations,
            is_mc, build_mc_strategy
        )
    else
        seed_strategy =
            unguided_seeding ?
            build_seed_strategy(param_set, domain, domain.seed_target_locations) :
            nothing
        fog_strategy =
            unguided_fogging ?
            build_fog_strategy(param_set, domain, domain.fog_target_locations) :
            nothing
        mc_strategy =
            unguided_moving_corals ?
            build_mc_strategy(param_set, domain, domain.mc_target_locations) :
            nothing
    end

    dhw_projection::Vector{Float64} = zeros(Float64, n_locs)
    wave_projection::Vector{Float64} = zeros(Float64, n_locs)

    # Set up distributions for natural adaptation/heritability
    c_mean_t_1::Array{Float64,3} = repeat(
        _to_group_size(domain.coral_growth, corals.dist_mean),
        1,
        1,
        n_locs
    )
    c_std::Array{Float64,2} = _to_group_size(
        domain.coral_growth, corals.dist_std
    )

    c_mean_t = copy(c_mean_t_1)

    # Log of distributions
    dhw_tol_mean_log = cache.dhw_tol_mean_log  # tmp log for mean dhw tolerances

    # Cache for proportional mortality and coral population increases
    bleaching_mort = zeros(tf, n_groups, n_sizes, n_locs)
    #### End coral constants

    ## Update ecological parameters based on intervention option
    # Pre-calculate proportion of survivers from wave stress
    # Sw_t = wave_damage!(cache.wave_damage, wave_scen, corals.wavemort90, n_species)

    area_weighted_conn = sparse(conn .* vec_abs_k)
    conn_cache = copy(area_weighted_conn)

    # basal_area_per_settler is the area in m^2 of a size class one coral
    basal_area_per_settler = colony_mean_area(
        corals.mean_colony_diameter_m[corals.class_id .== 1]
    )

    # Dummy vars to fill/replace with ranks of selected locations
    selected_seed_ranks = []
    selected_fog_ranks = []
    selected_mc_ranks = []

    # Cache matrix to store potential settlers
    potential_settlers = zeros(size(fec_scope)...)
    _linear_extensions = _to_group_size(domain.coral_growth, corals.linear_extension)
    _bin_edges = bin_edges()
    survival_rate = 1.0 .- _to_group_size(domain.coral_growth, corals.mb_rate)

    # Empty the old contents of the buffers and add the new blocks
    cover_view = [@view C_cover[1, :, :, loc] for loc in 1:n_locs]
    functional_groups = reuse_buffers!.(
        functional_groups, (cover_view .* vec_abs_k)
    )

    # Preallocate memory for temporaries
    survival_rate_cache = ones(n_groups, n_sizes, n_locs)
    C_cover_t::Array{Float64,3} = zeros(n_groups, n_sizes, n_locs)
    ΔC_cover_t = zeros(n_groups, n_sizes, n_locs)

    # Max projected cover is used on linear extensions scale factors
    habitable_max_projected_cover = max_projected_cover(
        _linear_extensions,
        _bin_edges,
        habitable_loc_areas
    )

    # Extract unique cb_calib_groups and create location masks for each biogroup
    unique_cb_calib_groups::Vector{Int64} = sort(unique(domain.loc_data.CB_CALIB_GROUPS))
    n_cb_calib_groups::Int64 = length(unique_cb_calib_groups)
    cb_calib_group_masks::BitMatrix = falses(n_locs, n_cb_calib_groups)
    for (idx, cb_calib_group) in enumerate(unique_cb_calib_groups)
        cb_calib_group_masks[:, idx] .= domain.loc_data.CB_CALIB_GROUPS .== cb_calib_group
    end

    # Index into unique_biogroups for each location
    loc_cb_calib_group_idxs::Vector{Int64} = [
        findfirst(x -> x == cb_group, unique_cb_calib_groups)
        for cb_group in domain.loc_data.CB_CALIB_GROUPS
    ]

    is_growth_acc_mask = occursin.("growth_acceleration", factor_names)
    growth_acc_steepness::Vector{Float64} =
        param_set[is_growth_acc_mask .&& occursin.("steepness", factor_names)].data
    growth_acc_height::Vector{Float64} =
        param_set[is_growth_acc_mask .&& occursin.("height", factor_names)].data
    growth_acc_midpoint::Vector{Float64} =
        param_set[is_growth_acc_mask .&& occursin.("midpoint", factor_names)].data

    # linear_extension scale factors with dimensions (cb_calib_groups ⋅ functional_groups)
    _linear_extension_scale_factors = reshape(
        param_set[occursin.("linear_extension_scale", factor_names)].data,
        n_cb_calib_groups,
        n_groups
    )

    # mb_rate scale factors with dimensions (cb_calib_groups ⋅ functional_groups)
    _mb_rate_scale_factors = reshape(
        param_set[occursin.("mb_rate_scale", factor_names)].data,
        n_cb_calib_groups,
        n_groups
    )

    biogrp_lin_ext::Array{Float64,3} = repeat(_linear_extensions, 1, 1, n_cb_calib_groups)
    biogrp_survival::Array{Float64,3} = repeat(survival_rate, 1, 1, n_cb_calib_groups)
    for i in 1:n_cb_calib_groups
        biogrp_lin_ext[:, :, i] .*= _linear_extension_scale_factors[i, :]
        biogrp_survival[:, :, i] .= apply_survival_scaling(
            biogrp_survival[:, :, i], _mb_rate_scale_factors[i, :]
        )
    end

    # Preallocate vector for growth constraints
    growth_constraints::Vector{Float64} = zeros(Float64, n_locs)

    FLoops.assistant(false)
    habitable_loc_idxs = findall(habitable_locs)

    # Growth constraints and acceleration caches and masks
    LIN_EXT_SCALE_FACTOR_THRESHOLD = 0.5
    relative_habitable_cover_cache = zeros(n_locs)
    growth_threshold_mask_cache::BitVector = trues(n_locs)
    cover_threshold_mask::BitVector = falses(n_locs)
    cache_habitable_max_projected_cover = copy(habitable_max_projected_cover)
    agg_cover_above_threshold_mask::BitVector = falses(n_habitable_locs)
    net_growth_rates = zeros(n_groups, n_sizes, n_locs)

    # Assume heat tolerance enhancement is based on population from `a_adapt_ref` years ago
    # If a_adapt == 0, use always first year as reference for a_adapt
    a_adapt_ref::Int64 = param_set[At("a_adapt_ref")]

    # c_mean tolerance of the "natural" population used as the reference for seeding corals
    # heat tolerance distribution
    c_mean_reference::Array{Float64,3} = if a_adapt_ref == 0
        # If a_adapt_ref == 0, 3rd dim holds c_mean at t-1 (idx 1), updated yearly, and
        # c_mean at t=1
        zeros(n_groups, n_locs, 2)
    else
        # If a_adapt_ref > 0, updates c_mean_reference yearly moving elements 1:end-1 to
        # 2:end and adds c_mean at t-1 (idx 1).
        zeros(n_groups, n_locs, a_adapt_ref)
    end

    # Keeps track of the heat tolerance means of juvenilles to use as a base for thermal
    # tolerance enhancement
    c_mean_reference .= copy(c_mean_t[:, 1, :])

    for tstep::Int64 in 2:tf
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

        habitable_max_projected_cover =
            cache_habitable_max_projected_cover .+
            dropdims(sum(recruitment; dims=1); dims=1)

        relative_habitable_cover_cache[habitable_locs] .=
            loc_coral_cover(C_cover_t[:, :, habitable_locs]) ./
            vec_abs_k[habitable_locs]

        # Growth constrains need to be calculated seperately for differen growth rates
        growth_threshold_mask_cache .=
            relative_habitable_cover_cache .>= LIN_EXT_SCALE_FACTOR_THRESHOLD
        for idx in 1:n_cb_calib_groups
            cover_threshold_mask .=
                cb_calib_group_masks[:, idx] .&& growth_threshold_mask_cache

            # Only apply linear_extension_scale_factors to locations with high cover
            growth_constraints[cover_threshold_mask] .= linear_extension_scale_factors(
                C_cover_t[:, :, cover_threshold_mask],
                vec_abs_k[cover_threshold_mask],
                biogrp_lin_ext[:, :, idx],
                _bin_edges,
                habitable_max_projected_cover[cover_threshold_mask]
            )
        end
        growth_threshold_mask_cache .=
            relative_habitable_cover_cache .< LIN_EXT_SCALE_FACTOR_THRESHOLD
        for idx in 1:n_cb_calib_groups
            cover_threshold_mask .=
                growth_threshold_mask_cache .&& cb_calib_group_masks[:, idx]
            growth_constraints[cover_threshold_mask] .=
                growth_acceleration.(
                    growth_acc_height[idx],
                    growth_acc_midpoint[idx],
                    growth_acc_steepness[idx],
                    relative_habitable_cover_cache[cover_threshold_mask]
                )
        end

        @inbounds for i in habitable_loc_idxs
            net_growth_rates[:, :, i] .=
                biogrp_lin_ext[:, :, loc_cb_calib_group_idxs[i]] .*
                growth_constraints[i]
            @views timestep!(
                functional_groups[i],
                recruitment[:, i],
                net_growth_rates[:, :, i],
                biogrp_survival[:, :, loc_cb_calib_group_idxs[i]]
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

        # Natural adaptation (doesn't change C_cover_t)
        if tstep <= tf
            adjust_DHW_distribution!(
                @view(C_cover[tstep - 1, :, :, :]), c_mean_t, net_growth_rates
            )

            # Set values for t to t-1
            c_mean_t_1 .= c_mean_t

            if in_debug_mode
                # Log dhw tolerances if in debug mode
                dhw_tol_mean_log[tstep, :, :] .= reshape(
                    mean.(c_mean_t), size(dhw_tol_mean_log)[2:3]
                )
            end
        end

        # Reproduction
        # Calculates scope for coral fedundity for each size class and at each location
        fecundity_scope!(fec_scope, fecundity_per_m², C_cover_t, habitable_areas)

        for l in 1:n_locs
            prop_fecundity[:, l] .= fec_scope[:, l] ./ sum(fec_scope[:, l])
        end

        _loc_coral_cover = loc_coral_cover(C_cover_t)
        leftover_space_m² = relative_leftover_space(_loc_coral_cover) .* vec_abs_k

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
                potential_settlers
            )[
                :, habitable_loc_idxs
            ] ./ habitable_areas[:, habitable_loc_idxs]

        # Reset fecundity scope before next run
        fec_scope .= 0.0

        leftover_space_m² .-= dropdims(sum(recruitment; dims=1); dims=1) .* vec_abs_k

        # Determine intervention locations whose deployment is assumed to occur
        # between November to February.
        #
        # Nominal sequence of events is conceptualised as:
        # - Corals spawn and are recruited
        # - SRM is applied
        # - Fogging is applied next
        # - Seeding interventions occur
        # - bleaching then occurs
        # - then cyclones hit
        #
        # Seeding occurs before bleaching as current tech only allows deployments shortly
        # after spawning. If bio-banking or similar tech comes up to speed then we could
        # potentially deploy at alternate times.

        # Shading
        # Apply regional cooling effect before selecting locations to seed
        dhw_t .= dhw_scen[tstep, :]  # subset of DHW for given timestep
        if apply_shading && shade_decision_years[tstep]
            shade_locs_mask = domain.loc_ids .∈ [domain.shade_target_locations]
            Yshade[tstep, :] .= srm

            # Apply reduction in DHW due to SRM
            dhw_t[shade_locs_mask] .= max.(0.0, dhw_t[shade_locs_mask] .- srm)
        end

        if is_guided
            # Use modified projected DHW (may have been affected by shading)
            dhw_p = copy(dhw_scen)
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
            current_loc_cover = dropdims(sum(C_cover_t; dims=(1, 2)); dims=(1, 2))
        end

        # Fogging
        # if is_guided
        #     if fog_decision_years[tstep] && (fogging .> 0.0)
        #         selected_fog_ranks = select_locations(
        #             fog_pref,
        #             decision_mat,
        #             MCDA_approach,
        #             min_iv_locs
        #         )

        #         if !isempty(selected_fog_ranks)
        #             log_location_ranks[tstep, At(selected_fog_ranks), At(:fog)] .=
        #                 1:length(selected_fog_ranks)
        #         end
        #     end
        # elseif apply_fogging && fog_decision_years[tstep]
        #     selected_fog_ranks = unguided_selection(
        #         domain.loc_ids,
        #         min_iv_locs,
        #         vec(leftover_space_m²)
        #         # depth_criteria
        #     )

        #     log_location_ranks[tstep, At(selected_fog_ranks), At(:fog)] .= 1.0
        # end

        # Fog location selection
        if !isnothing(fog_strategy)
            state = if is_guided
                build_state(
                    domain, fog_strategy,
                    (
                        current_cover=current_loc_cover,
                        recent_cover_losses=recent_cover_losses,
                        last_deployment=last_seed_deployment
                    )
                )
            else
                nothing
            end

            # Get candidate locations from strategy
            candidate_locs = filter_candidate_locations(fog_strategy, tstep, state)
            candidate_loc_indices = findall(
                in.(domain.loc_ids, Ref(candidate_locs))
            )

            if !isempty(candidate_locs)
                selected_fog_ranks = []
                if is_guided
                    # Update decision matrix with current conditions
                    update_criteria_values!(
                        fog_decision_mat[location=At(candidate_locs)];
                        heat_stress=dhw_projection[candidate_loc_indices],
                        wave_stress=wave_projection[candidate_loc_indices],
                        coral_cover=current_loc_cover[candidate_loc_indices],
                        in_connectivity=in_conn[candidate_loc_indices],
                        out_connectivity=out_conn[candidate_loc_indices]
                    )

                    # Build state for target locations only
                    selected_fog_ranks = select_locations(
                        fog_pref,
                        fog_decision_mat[location=At(candidate_locs)],
                        MCDA_approach,
                        min_iv_locs
                    )
                else
                    selected_fog_ranks = unguided_selection(
                        candidate_locs,
                        min_iv_locs,
                        vec(leftover_space_m²[candidate_loc_indices])
                    )
                end
                if !isempty(selected_fog_ranks)
                    log_val = is_guided ? (1:length(selected_fog_ranks)) : 1.0
                    log_location_ranks[tstep, At(selected_fog_ranks), At(:fog)] .=
                        log_val
                    last_fog_deployment[candidate_loc_indices] .= tstep
                end
            end
        end

        has_fog_locs::Bool = !isempty(selected_fog_ranks)

        # Fog selected locations
        if has_fog_locs  # fog_decision_years[tstep] &&
            fog_locs = findall(log_location_ranks.locations .∈ [selected_fog_ranks])
            fog_locations!(@view(Yfog[tstep, :]), fog_locs, dhw_t, fogging)

            # Empty selected_fog_ranks before the next iteration
            selected_fog_ranks = []
        end

        # Moving corals intervention
        if !isnothing(mc_strategy)
            state = if is_guided
                build_state(
                    domain, mc_strategy,
                    (
                        current_cover=current_loc_cover,
                        recent_cover_losses=recent_cover_losses,
                        last_deployment=last_seed_deployment
                    )
                )
            else
                nothing
            end

            # Get candidate locations from strategy
            candidate_locs = filter_candidate_locations(mc_strategy, tstep, state)
            candidate_loc_indices = findall(
                in.(domain.loc_ids, Ref(candidate_locs))
            )
            if !isempty(candidate_locs)
                if is_guided
                    # Update decision matrix with current conditions
                    update_criteria_values!(
                        mc_decision_mat[location=At(candidate_locs)];
                        heat_stress=dhw_projection[candidate_loc_indices],
                        wave_stress=wave_projection[candidate_loc_indices],
                        coral_cover=current_loc_cover[candidate_loc_indices],
                        in_connectivity=in_conn[candidate_loc_indices],
                        out_connectivity=out_conn[candidate_loc_indices]
                    )

                    # Build state for target locations only
                    selected_mc_ranks = select_locations(
                        mc_pref,
                        mc_decision_mat[location=At(candidate_locs)],
                        MCDA_approach,
                        min_iv_locs
                    )
                else
                    # Unguided deployment, seed/fog corals anywhere, so long as available space > 0
                    selected_mc_ranks = unguided_selection(
                        candidate_locs,
                        min_iv_locs,
                        vec(leftover_space_m²[candidate_loc_indices]),
                        depth_criteria[candidate_loc_indices]
                    )
                end
                if !isempty(selected_mc_ranks)
                    log_val = is_guided ? (1:length(selected_mc_ranks)) : 1.0
                    # log_location_ranks[tstep, At(selected_mc_ranks), At(:seed)] .= log_val
                    last_mc_deployment[candidate_loc_indices] .= tstep
                end
            end
        end

        # Check if locations are selected (can reuse previous selection)
        has_mc_locs::Bool = !isempty(selected_mc_ranks)

        # Apply Moving Corals (assumed to occur after spawning)
        if has_mc_locs
            # Selected locations can fill up over time so avoid locations with no space´
            mc_loc_idx = findall(domain.loc_ids .∈ [selected_mc_ranks])

            # Discount recruits added via seeding from leftover space
            @views available_space =
                leftover_space_m²[mc_loc_idx] .-
                dropdims(sum(recruitment[:, mc_loc_idx]; dims=1); dims=1)

            # available_space = leftover_space_m²[mc_loc_idx]
            locs_with_space = findall(available_space .> 0.0)

            # If there are locations with space to select from, then deploy what we can
            # Otherwise, do nothing.
            if length(locs_with_space) > 0
                # Calculate proportion to seed based on current available space
                mc_loc_idx = mc_loc_idx[locs_with_space]
                available_space = available_space[locs_with_space]

                @views mc_proportional_increase, n_mc_corals = distribute_moving_corals(
                    vec_abs_k[mc_loc_idx],
                    available_space,
                    n_mc_settlers,
                    colony_areas[_seed_size_groups],
                    prop_fecundity[:, mc_loc_idx]
                )

                @views recruitment[:, mc_loc_idx] .+= mc_proportional_increase

                # TODO log moving corals
                # Log estimated number of corals seeded
                Ymc[tstep, :, mc_loc_idx] .= n_mc_corals
            end

            # Empty selected_mc_ranks before the next iteration
            selected_mc_ranks = []
        end

        # Set next settlers DHW tolerance after reproducton and moving corals intervention
        settler_DHW_tolerance!(
            c_mean_t_1,
            c_mean_t,
            vec_abs_k,
            TP_data,  # ! IMPORTANT: Pass in transition probability matrix, not connectivity!
            recruitment,
            fecundity_per_m²,
            param_set[At("heritability")]
        )

        # Seeding
        # IDs of valid locations considering locations that have space for corals
        # locs_with_space = vec(leftover_space_m²) .> 0.0

        # if is_guided && seed_decision_years[tstep] && (length(considered_locs) > 0)
        #     considered_locs = findall(_valid_locs .& locs_with_space)

        #     # Use modified projected DHW (may have been affected by fogging or shading)
        #     dhw_p = copy(dhw_scen)
        #     dhw_p[tstep, :] .= dhw_t

        #     dhw_projection = weighted_projection(dhw_p, tstep, plan_horizon, decay, tf)
        #     wave_projection = weighted_projection(wave_scen, tstep, plan_horizon, decay, tf)

        #     # Determine connectivity strength weighting by area.
        #     # Accounts for strength of connectivity where there is low/no coral cover
        #     in_conn, out_conn, _ = connectivity_strength(
        #         area_weighted_conn, vec(_loc_coral_cover), conn_cache
        #     )

        #     update_criteria_values!(
        #         decision_mat;
        #         heat_stress=dhw_projection[_valid_locs],
        #         wave_stress=wave_projection[_valid_locs],
        #         coral_cover=_loc_coral_cover[_valid_locs],  # Coral cover relative to `k`
        #         in_connectivity=in_conn[_valid_locs],  # area weighted connectivities for time `t`
        #         out_connectivity=out_conn[_valid_locs]
        #     )

        #     selected_seed_ranks = select_locations(
        #         seed_pref,
        #         decision_mat[location=locs_with_space[_valid_locs]],
        #         MCDA_approach,
        #         considered_locs,
        #         min_iv_locs
        #     )

        #     # Log rankings as appropriate
        #     if !isempty(selected_seed_ranks)
        #         log_location_ranks[tstep, At(selected_seed_ranks), At(:seed)] .=
        #             1:length(selected_seed_ranks)
        #     end
        # elseif apply_seeding && seed_decision_years[tstep]
        #     # Unguided deployment, seed/fog corals anywhere, so long as available space > 0
        #     selected_seed_ranks = unguided_selection(
        #         domain.loc_ids,
        #         min_iv_locs,
        #         vec(leftover_space_m²),
        #         depth_criteria
        #     )

        #     log_location_ranks[tstep, At(selected_seed_ranks), At(:seed)] .= 1.0

        #     # Estimate proportional change in cover to apply to cubes
        # end

        # Seeding location selection
        if !isnothing(seed_strategy)
            state = if is_guided
                build_state(
                    domain, seed_strategy,
                    (
                        current_cover=current_loc_cover,
                        recent_cover_losses=recent_cover_losses,
                        last_deployment=last_seed_deployment
                    )
                )
            else
                nothing
            end

            # Get candidate locations from strategy
            candidate_locs = filter_candidate_locations(seed_strategy, tstep, state)
            candidate_loc_indices = findall(
                in.(domain.loc_ids, Ref(candidate_locs))
            )
            if !isempty(candidate_locs)
                if is_guided
                    # Update decision matrix with current conditions
                    update_criteria_values!(
                        seed_decision_mat[location=At(candidate_locs)];
                        heat_stress=dhw_projection[candidate_loc_indices],
                        wave_stress=wave_projection[candidate_loc_indices],
                        coral_cover=current_loc_cover[candidate_loc_indices],
                        in_connectivity=in_conn[candidate_loc_indices],
                        out_connectivity=out_conn[candidate_loc_indices]
                    )

                    # Build state for target locations only
                    selected_seed_ranks = select_locations(
                        seed_pref,
                        seed_decision_mat[location=At(candidate_locs)],
                        MCDA_approach,
                        min_iv_locs
                    )
                else
                    # Unguided deployment, seed/fog corals anywhere, so long as available space > 0
                    selected_seed_ranks = unguided_selection(
                        candidate_locs,
                        min_iv_locs,
                        vec(leftover_space_m²[candidate_loc_indices]),
                        depth_criteria[candidate_loc_indices]
                    )
                end
                if !isempty(selected_seed_ranks)
                    log_val = is_guided ? (1:length(selected_seed_ranks)) : 1.0
                    log_location_ranks[tstep, At(selected_seed_ranks), At(:seed)] .=
                        log_val
                    last_seed_deployment[candidate_loc_indices] .= tstep
                end
            end
        end

        # Check if locations are selected (can reuse previous selection)
        has_seed_locs::Bool = !isempty(selected_seed_ranks)

        # Apply seeding (assumed to occur after spawning)
        if has_seed_locs  # seed_decision_years[tstep] &&
            # Seed selected locations
            # Selected locations can fill up over time so avoid locations with no space´
            # Main.@infiltrate
            seed_loc_idx = findall(domain.loc_ids .∈ [selected_seed_ranks])

            available_space = leftover_space_m²[seed_loc_idx]
            locs_with_space = findall(available_space .> 0.0)

            # If there are locations with space to select from, then deploy what we can
            # Otherwise, do nothing.
            if length(locs_with_space) > 0
                # Calculate proportion to seed based on current available space
                seed_loc_idx = seed_loc_idx[locs_with_space]
                available_space = available_space[locs_with_space]

                # Extract colony areas and determine approximate seeded area in m^2
                seeded_area = colony_areas[_seed_size_groups] .* seed_volume

                proportional_increase, n_corals_seeded = distribute_seeded_corals(
                    vec_abs_k[seed_loc_idx],
                    available_space,
                    seed_volume.data,
                    seeded_area,
                    seeding_devices_per_m2
                )

                # Log estimated number of corals seeded
                Yseed[tstep, :, seed_loc_idx] .= n_corals_seeded'

                # Add coral seeding to recruitment
                recruitment[:, seed_loc_idx] .+= proportional_increase

                update_tolerance_distribution!(
                    proportional_increase,
                    C_cover_t,
                    c_mean_t,
                    c_mean_reference[:, :, end],
                    c_std,
                    seed_loc_idx,
                    _seed_size_groups,
                    a_adapt
                )
            end

            # Empty selected_seed_ranks before the next iteration
            selected_seed_ranks = []
        end

        # Apply disturbances

        # ΔC_cover_t should only hold changes in cover due to env. disturbances
        ΔC_cover_t .= copy(C_cover_t)

        # Bleaching
        # Calculate and apply bleaching mortality
        # Bleaching typically occurs in the warmer months (November - February)
        #    This: `dhw_t .* (1.0 .- wave_scen[tstep, :])`
        #    attempts to account for the cooling effect of storms / high wave activity
        # `wave_scen` is normalized to the maximum value found for the given wave scenario
        # so what causes 100% mortality can differ between runs.
        bleaching_mortality!(
            C_cover_t,
            dhw_t,  # collect(dhw_t .* (1.0 .- @view(wave_scen[tstep, :]))),
            depth_coeff,
            c_std,
            c_mean_t_1,
            c_mean_t,
            @view(bleaching_mort[(tstep - 1):tstep, :, :, :])
        )

        if a_adapt_ref > 0
            c_mean_reference[:, :, 2:end] .= c_mean_reference[:, :, 1:(end - 1)]
        end
        c_mean_reference[:, :, 1] .= c_mean_t[:, 1, :]

        # Coral deaths due to selected cyclone scenario
        # Peak cyclone period is January to March
        # TODO: Update cyclone data to hold data for relevant functional groups
        cyclone_mortality!(C_cover_t, cyclone_mortality_scen[tstep, :, :]')

        # Calculate survival_rate due to env. disturbances
        no_mortality_mask = ΔC_cover_t .== 0.0
        survival_rate_cache[.!no_mortality_mask] .= (C_cover_t ./ ΔC_cover_t)[.!no_mortality_mask]
        survival_rate_cache[no_mortality_mask] .= 1.0
        @assert sum(survival_rate_cache .> 1) == 0 "Survival rate should be <= 1"

        for loc in 1:n_locs
            apply_mortality!(
                functional_groups[loc], @view(survival_rate_cache[:, :, loc])
            )
        end

        recruitment .*= (view(survival_rate_cache, :, 1, :) .* habitable_areas)

        C_cover[tstep, :, :, :] .= C_cover_t

        # Track cover loss for reactive strategies
        _is_reactive = any(
            is_reactive(param_set[At(["seed_strategy", "fog_strategy", "mc_strategy"])])
        )
        if is_guided && _is_reactive && (tstep > 1)
            # Calculate proportional cover loss at each location
            # due to disturbances this timestep
            Δcover_loss_proportion = zeros(n_locs)
            for loc in 1:n_locs
                cover_before = sum(@view(ΔC_cover_t[:, :, loc]))
                cover_after = sum(@view(C_cover_t[:, :, loc]))

                if cover_before > 0.0
                    Δcover_loss_proportion[loc] =
                        (cover_before - cover_after) / cover_before
                end
            end

            # Store in rolling buffer
            push!(recent_cover_losses, Δcover_loss_proportion)
        end
    end

    # Could collate critical DHW threshold log for corals to reduce disk space...
    # dhw_tol_mean = dropdims(mean(dhw_tol_mean_log, dims=3), dims=3)
    # dhw_tol_mean_std = dropdims(mean(dhw_tol_std_log, dims=3), dims=3)
    # collated_dhw_tol_log = DataCube(cat(dhw_tol_mean, dhw_tol_mean_std, dims=3);
    #     timesteps=1:tf, species=corals.coral_id, stat=[:mean, :stdev])
    if in_debug_mode
        collated_dhw_tol_log = DataCube(
            dhw_tol_mean_log; timesteps=1:tf, species=corals.coral_id, sites=1:n_locs
        )
    else
        collated_dhw_tol_log = false
    end

    # Set variables to nothing so garbage collector clears them
    # Leads to memory leak issues in multiprocessing contexts without these.
    wave_scen = nothing
    dhw_tol_mean_log = nothing

    return (
        raw=C_cover,
        seed_log=Yseed,
        mc_log=Ymc,
        fog_log=Yfog,
        shade_log=Yshade,
        site_ranks=log_location_ranks,
        bleaching_mortality=bleaching_mort,
        coral_dhw_log=collated_dhw_tol_log
    )
end
