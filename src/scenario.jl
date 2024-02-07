"""Scenario running functions"""

using ADRIA.metrics: relative_cover, relative_loc_taxa_cover, total_absolute_cover, absolute_shelter_volume, relative_shelter_volume
using ADRIA.metrics: relative_juveniles, relative_taxa_cover, juvenile_indicator
using ADRIA.metrics: coral_evenness
using ADRIA.decision:
    within_depth_bounds,
    summary_stat_env,
    DMCDA_vars,
    guided_site_selection,
    unguided_site_selection

"""
    setup_cache(domain::Domain)::NamedTuple

Establish tuple of matrices/vectors for use as reusable data stores to avoid repeated memory allocations.
"""
function setup_cache(domain::Domain)::NamedTuple

    # Simulation constants
    n_sites::Int64 = domain.coral_growth.n_sites
    n_species::Int64 = domain.coral_growth.n_species
    n_groups::Int64 = domain.coral_growth.n_groups
    tf = length(timesteps(domain))

    cache = (
        # sf=zeros(n_groups, n_sites),  # stressed fecundity, commented out as it is disabled
        fec_all=zeros(n_species, n_sites),  # all fecundity
        fec_scope=zeros(n_groups, n_sites),  # fecundity scope
        recruitment=zeros(n_groups, n_sites),  # coral recruitment
        dhw_step=zeros(n_sites),  # DHW for each time step
        cov_tmp=zeros(n_species, n_sites),  # Cover for previous timestep
        depth_coeff=zeros(n_sites),  # store for depth coefficient
        site_area=Matrix{Float64}(site_area(domain)'),  # area of locations
        site_k_area=Matrix{Float64}(site_k_area(domain)'),  # location carrying capacity
        wave_damage=zeros(tf, n_species, n_sites),  # damage coefficient for each size class
        dhw_tol_mean_log=zeros(tf, n_species, n_sites),  # tmp log for mean dhw tolerances
    )

    return cache
end


function run_scenarios(param_df::DataFrame, domain::Domain, RCP::String; show_progress=true, remove_workers=true)::ResultSet
    msg = """
    `run_scenarios(param_df, domain, RCP)` is now deprecated and will be removed in
    ADRIA v1.0

    Instead, use:
        `run_scenarios(dom, scens, RCP)`
    """
    @warn msg
    return run_scenarios(domain, param_df, [RCP]; show_progress, remove_workers)
end
function run_scenarios(param_df::DataFrame, domain::Domain, RCP::Vector{String}; show_progress=true, remove_workers=true)::ResultSet
    msg = """
    `run_scenarios(param_df, domain, RCP)` is now deprecated and will be removed in
    ADRIA v1.0

    Instead, use:
        `run_scenarios(dom, scens, RCP)`
    """
    @warn msg
    return run_scenarios(domain, param_df, RCP; show_progress, remove_workers)
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
    dom::Domain,
    scens::DataFrame,
    RCP::String;
    show_progress=true,
    remove_workers=true
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
    rcps_df = DataFrame(RCP=parse.(Int64, RCP))
    scenarios_df = crossjoin(scens, rcps_df)
    sort!(scenarios_df, :RCP)

    @info "Setting up Result Set"
    dom, data_store = ADRIA.setup_result_store!(dom, scenarios_df)

    # Convert DataFrame to named matrix for faster iteration
    scenarios_matrix = NamedDimsArray(
        Matrix(scenarios_df);
        scenarios=1:nrow(scenarios_df),
        factors=names(scenarios_df)
    )

    para_threshold = typeof(dom) == RMEDomain ? 8 : 256
    parallel = (parse(Bool, ENV["ADRIA_DEBUG"]) == false) && (nrow(scens) >= para_threshold)
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
                func = (dfx) -> run_scenario(dfx..., data_store)
            end
        end

        @info "Time taken to spin up workers: $(round(spinup_time; digits=2)) seconds"

        # Define local helper
        func = (dfx) -> run_scenario(dfx..., data_store)
    end

    if parallel
        for rcp in RCP
            run_msg = "Running $(nrow(scens)) scenarios for RCP $rcp"

            # Switch RCPs so correct data is loaded
            target_rows = findall(scenarios_matrix("RCP") .== parse(Float64, rcp))
            rep_doms = Iterators.repeated(dom, length(target_rows))
            scenario_args = zip(rep_doms, target_rows, eachrow(scenarios_matrix[target_rows, :]))
            if show_progress
                @showprogress desc=run_msg dt=4 pmap(func, CachingPool(workers()), scenario_args)
            else
                pmap(func, CachingPool(workers()), scenario_args)
            end
        end
    else
        # Define local helper
        func = dfx -> run_scenario(dfx..., data_store)

        for rcp in RCP
            run_msg = "Running $(nrow(scens)) scenarios for RCP $rcp"

            # Switch RCPs so correct data is loaded
            dom = switch_RCPs!(dom, rcp)
            target_rows = findall(scenarios_matrix("RCP") .== parse(Float64, rcp))
            rep_doms = Iterators.repeated(dom, size(scenarios_matrix, 1))
            scenario_args = zip(rep_doms, target_rows, eachrow(scenarios_matrix[target_rows, :]))
            if show_progress
                @showprogress desc=run_msg dt=4 map(func, scenario_args)
            else
                map(func, scenario_args)
            end
        end
    end

    if parallel && remove_workers
        _remove_workers()
    end

    return load_results(_result_location(dom, RCP))
end

"""
    run_scenario(domain::Domain, idx::Int64, scenario::Union{AbstractVector, DataFrameRow}, data_store::NamedTuple)::Nothing
    run_scenario(domain::Domain, idx::Int64, scenario::Union{AbstractVector, DataFrameRow}, domain::Domain, data_store::NamedTuple)::Nothing
    run_scenario(domain::Domain, scenario::Union{AbstractVector, DataFrameRow})::NamedTuple
    run_scenario(domain::Domain, scenario::NamedTuple)::NamedTuple

Run individual scenarios for a given domain, saving results to a Zarr data store.
Results are stored in Zarr format at a pre-configured location.
Sets up a new `cache` if not provided.

# Notes
Logs of site ranks only store the mean site rankings over all environmental scenarios.
This is to reduce the volume of data stored.

# Returns
Nothing
"""
function run_scenario(
    domain::Domain,
    idx::Int64,
    scenario::Union{AbstractVector,DataFrameRow},
    data_store::NamedTuple
)::Nothing
    if domain.RCP == ""
        local rcp
        try
            rcp = scenario("RCP")  # Try extracting from NamedDimsArray
        catch err
            if !(err isa MethodError)
                rethrow(err)
            end

            rcp = scenario.RCP  # Extract from dataframe
        end

        domain = switch_RCPs!(domain, string(Int64(rcp)))
    end

    result_set = run_model(domain, scenario)

    # Capture results to disk
    # Set values below threshold to 0 to save space
    threshold = parse(Float32, ENV["ADRIA_THRESHOLD"])

    rs_raw = result_set.raw
    vals = relative_cover(rs_raw)
    vals[vals.<threshold] .= 0.0
    data_store.relative_cover[:, :, idx] .= vals

    vals .= absolute_shelter_volume(rs_raw, site_k_area(domain), scenario)
    vals[vals.<threshold] .= 0.0
    data_store.absolute_shelter_volume[:, :, idx] .= vals

    vals .= relative_shelter_volume(rs_raw, site_k_area(domain), scenario)
    vals[vals.<threshold] .= 0.0
    data_store.relative_shelter_volume[:, :, idx] .= vals

    coral_spec::DataFrame = to_coral_spec(scenario)
    vals .= relative_juveniles(rs_raw, coral_spec)
    vals[vals.<threshold] .= 0.0
    data_store.relative_juveniles[:, :, idx] .= vals

    vals .= juvenile_indicator(rs_raw, coral_spec, site_k_area(domain))
    vals[vals.<threshold] .= 0.0
    data_store.juvenile_indicator[:, :, idx] .= vals

    vals = relative_taxa_cover(rs_raw, site_k_area(domain))
    vals[vals.<threshold] .= 0.0
    data_store.relative_taxa_cover[:, :, idx] .= vals

    vals = relative_loc_taxa_cover(rs_raw, site_k_area(domain))
    vals = coral_evenness(NamedDims.unname(vals))
    vals[vals.<threshold] .= 0.0
    data_store.coral_evenness[:, :, idx] .= vals

    # Store raw results if no metrics specified
    # if length(metrics) == 0
    #     data_store.raw[:, :, :, idx] .= r.raw
    # end

    # Store logs
    c_dim = Base.ndims(result_set.raw) + 1
    log_stores = (:site_ranks, :seed_log, :fog_log, :shade_log, :coral_dhw_log)
    for k in log_stores
        if k == :seed_log || k == :site_ranks
            concat_dim = c_dim
        else
            concat_dim = c_dim - 1
        end

        vals = getfield(result_set, k)

        try
            vals[vals.<threshold] .= 0.0
        catch err
            err isa MethodError ? nothing : rethrow(err)
        end

        if k == :seed_log
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

    if (idx % 256) == 0
        @everywhere GC.gc()
    end

    return nothing
end
function run_scenario(
    domain::Domain,
    scenario::Union{AbstractVector,DataFrameRow}
)::NamedTuple
    return run_model(domain, scenario)
end
function run_scenario(
    domain::Domain,
    scenario::Union{AbstractVector,DataFrameRow},
    RCP::String
)::NamedTuple
    domain = switch_RCPs!(domain, RCP)
    return run_scenario(domain, scenario)
end

"""
    run_model(domain::Domain, param_set::Union{NamedTuple,DataFrameRow})::NamedTuple

Core scenario running function.

# Notes
Only the mean site rankings are kept

# Returns
NamedTuple of collated results
"""
function run_model(domain::Domain, param_set::DataFrameRow)::NamedTuple
    ps = NamedDimsArray(Vector(param_set), factors=names(param_set))
    return run_model(domain, ps)
end
function run_model(domain::Domain, param_set::NamedDimsArray)::NamedTuple
    p = domain.coral_growth.ode_p
    corals = to_coral_spec(param_set)
    cache = setup_cache(domain)

    # Set random seed using intervention values
    # TODO: More robust way of getting intervention/criteria values
    rnd_seed_val::Int64 = floor(Int64, sum(param_set(!=("RCP"))))  # select everything except RCP
    Random.seed!(rnd_seed_val)

    dhw_idx::Int64 = Int64(param_set("dhw_scenario"))
    wave_idx::Int64 = Int64(param_set("wave_scenario"))
    cyclone_mortality_idx::Int64 = Int64(param_set("cyclone_mortality_scenario"))

    # Extract environmental data
    dhw_scen = @view(domain.dhw_scens[:, :, dhw_idx])

    # TODO: Better conversion of Ub to wave mortality
    #       Currently scaling significant wave height by its max to non-dimensionalize values
    wave_scen = copy(domain.wave_scens[:, :, wave_idx])
    wave_scen .= wave_scen ./ maximum(wave_scen)
    replace!(wave_scen, Inf=>0.0, NaN=>0.0)

    cyclone_mortality_scen = @view(
        domain.cyclone_mortality_scens[:, :, :, cyclone_mortality_idx]
    )

    tspan::Tuple = (0.0, 1.0)
    solver::Euler = Euler()

    MCDA_approach::Int64 = param_set("guided")

    # Environment variables are stored as strings, so convert to bool for use
    in_debug_mode = parse(Bool, ENV["ADRIA_DEBUG"]) == true

    # Sim constants
    sim_params = domain.sim_constants
    tf::Int64 = size(dhw_scen, 1)
    n_site_int::Int64 = sim_params.n_site_int
    n_locs::Int64 = domain.coral_growth.n_sites
    n_species::Int64 = domain.coral_growth.n_species
    n_groups::Int64 = domain.coral_growth.n_groups

    # Years to start seeding/shading/fogging
    seed_start_year::Int64 = param_set("seed_year_start")
    shade_start_year::Int64 = param_set("shade_year_start")
    fog_start_year::Int64 = param_set("fog_year_start")

    fogging::Real = param_set("fogging")  # proportion of bleaching mortality reduction through fogging
    srm::Real = param_set("SRM")  # DHW equivalents reduced by some shading mechanism
    seed_years::Int64 = param_set("seed_years")  # number of years to seed
    shade_years::Int64 = param_set("shade_years")  # number of years to shade
    fog_years::Int64 = param_set("fog_years")  # number of years to fog

    loc_k_area::Matrix{Float64} = cache.site_k_area
    fec_params_per_m²::Vector{Float64} = corals.fecundity  # number of larvae produced per m²

    # Caches
    TP_data = domain.TP_data
    # sf = cache.sf  # unused as it is currently deactivated
    fec_all = cache.fec_all
    fec_scope = cache.fec_scope
    recruitment = cache.recruitment
    dhw_t = cache.dhw_step
    C_t = cache.cov_tmp
    depth_coeff = cache.depth_coeff

    site_data = domain.site_data
    depth_coeff .= depth_coefficient.(site_data.depth_med)

    # Coral cover relative to available area (i.e., 1.0 == site is filled to max capacity)
    C_cover::Array{Float64,3} = zeros(tf, n_species, n_locs)
    C_cover[1, :, :] .= domain.init_coral_cover
    cover_tmp = zeros(n_locs)

    # Locations that can support corals
    valid_locs::BitVector = location_k(domain) .> 0.0
    site_ranks = SparseArray(zeros(tf, n_locs, 2))  # log seeding/fogging/shading ranks
    Yshade = SparseArray(spzeros(tf, n_locs))
    Yfog = SparseArray(spzeros(tf, n_locs))
    Yseed = SparseArray(zeros(tf, 3, n_locs))  # 3 = the number of seeded coral types

    # Prep scenario-specific flags/values
    # Intervention strategy: < 0 is no intervention, 0 is random location selection, > 0 is guided
    is_guided = param_set("guided") > 0

    # Decisions should place more weight on environmental conditions
    # closer to the decision point
    α = 0.95
    decay = α .^ (1:Int64(param_set("plan_horizon"))+1)

    # Years at which intervention locations are re-evaluated
    seed_decision_years = fill(false, tf)
    shade_decision_years = fill(false, tf)
    fog_decision_years = fill(false, tf)

    seed_start_year = max(seed_start_year, 2)
    if param_set("seed_freq") > 0
        max_consider = min(seed_start_year + seed_years - 1, tf)
        seed_decision_years[seed_start_year:Int64(param_set("seed_freq")):max_consider] .= true
    else
        # Start at year 2 or the given specified seed start year
        seed_decision_years[seed_start_year] = true
    end

    shade_start_year = max(shade_start_year, 2)
    if param_set("shade_freq") > 0
        max_consider = min(shade_start_year + shade_years - 1, tf)
        shade_decision_years[shade_start_year:Int64(param_set("shade_freq")):max_consider] .= true
    else
        # Start at year 2 or the given specified shade start year
        shade_decision_years[shade_start_year] = true
    end

    fog_start_year = max(fog_start_year, 2)
    if param_set("fog_freq") > 0
        max_consider = min(fog_start_year + fog_years - 1, tf)
        fog_decision_years[fog_start_year:Int64(param_set("fog_freq")):max_consider] .= true
    else
        # Start at year 2 or the given specified fog start year
        fog_decision_years[fog_start_year] = true
    end

    seed_locs::Vector{Int64} = zeros(Int64, n_site_int)
    fog_locs::Vector{Int64} = zeros(Int64, n_site_int)

    # Define taxa and size class to seed, and identify their factor names
    taxa_to_seed = [2, 3, 5]
    target_class_id::BitArray = corals.class_id .== 2  # seed second smallest size class
    taxa_names = param_set.factors[occursin.("N_seed_", param_set.factors)]

    # Identify taxa and size class to be seeded
    seed_sc = (corals.taxa_id .∈ [taxa_to_seed]) .& target_class_id

    # Extract colony areas for sites selected in m^2 and add adaptation values
    colony_areas = colony_mean_area(corals.mean_colony_diameter_m)
    seeded_area = colony_areas[seed_sc] .* param_set(taxa_names)

    # Set up assisted adaptation values
    a_adapt = zeros(n_species)
    a_adapt[seed_sc] .= param_set("a_adapt")

    # Flag indicating whether to seed or not to seed
    seed_corals = any(param_set(taxa_names) .> 0.0)
    # Flag indicating whether to fog or not fog
    apply_fogging = fogging > 0.0
    # Flag indicating whether to apply shading
    apply_shading = srm > 0.0

    # Defaults to considering all sites if depth cannot be considered.
    depth_priority = collect(1:n_locs)

    # Calculate total area to seed respecting tolerance for minimum available space to still
    # seed at a site
    area_to_seed = sum(seeded_area)

    # Filter out sites outside of desired depth range
    if .!all(site_data.depth_med .== 0)
        max_depth::Float64 = param_set("depth_min") + param_set("depth_offset")
        depth_criteria::BitArray{1} = within_depth_bounds(
            site_data.depth_med, max_depth, param_set("depth_min")
        )

        if any(depth_criteria .> 0)
            # If sites can be filtered based on depth, do so.
            depth_priority = depth_priority[depth_criteria]
        else
            # Otherwise if no sites can be filtered, remove depth as a criterion.
            @warn "No sites within provided depth range of $(param_set("depth_min")) - $(max_depth) meters. Considering all sites."
        end
    end

    if is_guided
        # Pre-allocate rankings
        rankings = [depth_priority zeros(Int, length(depth_priority)) zeros(Int, length(depth_priority))]

        # Prep site selection
        mcda_vars = DMCDA_vars(domain, param_set, depth_priority, sum(C_cover[1, :, :], dims=1), area_to_seed)

        # Number of time steps in environmental layers to look ahead when making decisions
        plan_horizon::Int64 = Int64(param_set("plan_horizon"))
    end

    # Set up distributions for natural adaptation/heritability
    c_mean_t_1::Matrix{Float64} = repeat(corals.dist_mean, 1, n_locs)
    c_mean_t = copy(c_mean_t_1)

    # Log of distributions
    dhw_tol_mean_log = cache.dhw_tol_mean_log  # tmp log for mean dhw tolerances

    # Cache for proportional mortality and coral population increases
    bleaching_mort = zeros(tf, n_species, n_locs)

    #### End coral constants

    ## Update ecological parameters based on intervention option

    # Treat as enhancement from mean of "natural" DHW tolerance
    a_adapt[a_adapt.>0.0] .+= corals.dist_mean[a_adapt.>0.0]

    # Pre-calculate proportion of survivers from wave stress
    # Sw_t = wave_damage!(cache.wave_damage, wave_scen, corals.wavemort90, n_species)

    p.r .= corals.growth_rate
    p.mb .= corals.mb_rate
    ode_u = zeros(n_species, n_locs)
    growth::ODEProblem = ODEProblem{true}(growthODE, ode_u, tspan, p)
    alg_hint = Symbol[:nonstiff]

    area_weighted_TP = TP_data .* site_k_area(domain)
    TP_cache = similar(area_weighted_TP)

    # basal_area_per_settler is the area in m^2 of a size class one coral
    basal_area_per_settler = colony_mean_area(corals.mean_colony_diameter_m[corals.class_id.==1])

    # Cache matrix to store potential settlers
    potential_settlers = zeros(size(fec_scope)...)
    for tstep::Int64 in 2:tf
        # Copy cover for previous timestep as basis for current timestep
        C_t .= C_cover[tstep - 1, :, :]

        # Coral deaths due to selected cyclone scenario
        # Peak cyclone period is January to March
        cyclone_mortality!(@views(C_t), p, cyclone_mortality_scen[tstep, :, :]')

        # Calculates scope for coral fedundity for each size class and at each location
        fecundity_scope!(fec_scope, fec_all, fec_params_per_m², C_t, loc_k_area)

        loc_coral_cover = sum(C_t, dims=1)  # dims: 1 * nsites
        leftover_space_m² = relative_leftover_space(loc_coral_cover) .* loc_k_area

        # Reset potential settlers to zero
        potential_settlers .= 0.0
        recruitment .= 0.0

        # Recruitment represents additional cover, relative to total site area
        # Recruitment/settlement occurs after the full moon in October/November
        recruitment[:, valid_locs] .= settler_cover(
            fec_scope,
            TP_data,
            leftover_space_m²,
            sim_params.max_settler_density,
            sim_params.max_larval_density,
            basal_area_per_settler,
            potential_settlers
        )[:, valid_locs] ./ loc_k_area[:, valid_locs]

        settler_DHW_tolerance!(
            c_mean_t_1,
            c_mean_t,
            site_k_area(domain),
            TP_data,
            recruitment,
            fec_params_per_m²,
            param_set("heritability")
        )

        # Add recruits to current cover
        C_t[p.small, :] .= recruitment

        # Check whether current timestep is in deployment period for each intervention
        in_fog_timeframe = fog_start_year <= tstep <= (fog_start_year + fog_years - 1)
        in_shade_timeframe =
            shade_start_year <= tstep <= (shade_start_year + shade_years - 1)
        in_seed_timeframe = seed_start_year <= tstep <= (seed_start_year + seed_years - 1)

        # Apply regional cooling effect before selecting locations to seed
        dhw_t .= dhw_scen[tstep, :]  # subset of DHW for given timestep
        if apply_shading && in_shade_timeframe
            Yshade[tstep, :] .= srm

            # Apply reduction in DHW due to SRM
            dhw_t .= max.(0.0, dhw_t .- srm)
        end

        # Determine intervention locations whose deployment is assumed to occur
        # between November to February.
        # - SRM is applied first
        # - Fogging is applied next
        # - Bleaching then occurs
        # - Then intervention locations are seeded
        if is_guided && (in_seed_timeframe || in_fog_timeframe)
            # Update dMCDA values

            # Determine subset of data to select data for planning horizon
            horizon::UnitRange{Int64} = tstep:min(tstep + plan_horizon, tf)
            d_s::UnitRange{Int64} = 1:length(horizon)

            # Put more weight on projected conditions closer to the decision point
            @views env_horizon = decay[d_s] .* dhw_scen[horizon, :]
            mcda_vars.heat_stress_prob .= summary_stat_env(env_horizon, :timesteps)

            @views env_horizon = decay[d_s] .* wave_scen[horizon, :]
            mcda_vars.dam_prob .= summary_stat_env(env_horizon, 1)
            mcda_vars.leftover_space .= leftover_space_m²

            # Determine connectivity strength
            # Account for cases where there is no coral cover
            in_conn, out_conn, strong_pred = connectivity_strength(area_weighted_TP, vec(loc_coral_cover), TP_cache)
            (seed_locs, fog_locs, rankings) = guided_site_selection(
                mcda_vars,
                MCDA_approach,
                seed_decision_years[tstep],
                fog_decision_years[tstep],
                seed_locs,
                fog_locs,
                rankings,
                in_conn[mcda_vars.site_ids],
                out_conn[mcda_vars.site_ids],
                strong_pred[mcda_vars.site_ids],
            )

            # Log site ranks
            # First col only holds site index ids so skip (with 2:end)
            site_ranks[tstep, rankings[:, 1], :] = rankings[:, 2:end]
        elseif (seed_corals && in_seed_timeframe) || (apply_fogging && in_fog_timeframe)
            # Unguided deployment, seed/fog corals anywhere, so long as available space > 0
            seed_locs, fog_locs = unguided_site_selection(
                seed_locs,
                fog_locs,
                seed_decision_years[tstep],
                fog_decision_years[tstep],
                n_site_int, vec(leftover_space_m²), depth_priority)

            site_ranks[tstep, seed_locs, 1] .= 1.0
            site_ranks[tstep, fog_locs, 2] .= 1.0
        end

        has_fog_locs::Bool = !all(fog_locs .== 0)

        # Check if locations are selected, and selected locations have space,
        # otherwise no valid locations were selected for seeding.
        has_seed_locs::Bool = true
        if !all(seed_locs .== 0)
            has_seed_locs = !all(leftover_space_m²[seed_locs] .== 0.0)
        else
            has_seed_locs = false
        end

        # Fog selected locations
        if apply_fogging && in_fog_timeframe && has_fog_locs
            fog_locations!(@view(Yfog[tstep, :]), fog_locs, dhw_t, fogging)
        end

        # Calculate and apply bleaching mortality
        # Bleaching typically occurs in the warmer months (November - February)
        #    This: `dhw_t .* (1.0 .- wave_scen[tstep, :])`
        #    attempts to account for the cooling effect of storms / high wave activity
        # `wave_scen` is normalized to the maximum value found for the given wave scenario
        # so what causes 100% mortality can differ between runs.
        bleaching_mortality!(
            C_t,
            collect(dhw_t .* (1.0 .- @view(wave_scen[tstep, :]))),
            depth_coeff,
            corals.dist_std,
            c_mean_t_1,
            c_mean_t,
            @view(bleaching_mort[tstep-1:tstep, :, :])
        )

        # Apply seeding
        # Assumes coral seeding occurs in the months after disturbances
        # (such as cyclones/bleaching).
        if seed_corals && in_seed_timeframe && has_seed_locs
            # Seed selected locations
            seed_corals!(C_t, vec(loc_k_area), vec(leftover_space_m²),
                seed_locs,
                seeded_area,
                seed_sc,
                a_adapt,
                @view(Yseed[tstep, :, :]),
                corals.dist_std,
                c_mean_t,
            )
        end

        # Update initial condition
        growth.u0 .= C_t
        sol::ODESolution = solve(growth, solver, save_everystep=false, save_start=false,
            alg_hints=alg_hint, dt=1.0)

        # Assign results
        C_cover[tstep, :, valid_locs] .= sol.u[end][:, valid_locs]

        # TODO:
        # Check if size classes are inappropriately out-growing available space
        proportional_adjustment!(
            @view(C_cover[tstep, :, valid_locs]),
            cover_tmp[valid_locs]
        )

        if tstep <= tf
            # Natural adaptation
            adjust_DHW_distribution!(
                @view(C_cover[tstep-1:tstep, :, :]),
                n_groups,
                c_mean_t,
                p.r
            )

            if in_debug_mode
                # Log dhw tolerances if in debug mode
                dhw_tol_mean_log[tstep, :, :] .= mean.(c_mean_t)
            end

            c_mean_t_1 .= c_mean_t
        end
    end

    # Could collate critical DHW threshold log for corals to reduce disk space...
    # dhw_tol_mean = dropdims(mean(dhw_tol_mean_log, dims=3), dims=3)
    # dhw_tol_mean_std = dropdims(mean(dhw_tol_std_log, dims=3), dims=3)
    # collated_dhw_tol_log = NamedDimsArray(cat(dhw_tol_mean, dhw_tol_mean_std, dims=3),
    #     timesteps=1:tf, species=corals.coral_id, stat=[:mean, :stdev])
    if in_debug_mode
        collated_dhw_tol_log = NamedDimsArray(
            dhw_tol_mean_log,
            timesteps=1:tf,
            species=corals.coral_id,
            sites=1:n_locs,
        )
    else
        collated_dhw_tol_log = false
    end

    # Set variables to nothing so garbage collector clears them
    # Leads to memory leak issues in multiprocessing contexts without these.
    wave_scen = nothing
    dhw_tol_mean_log = nothing

    # Avoid placing importance on sites that were not considered
    # (lower values are higher importance)
    site_ranks[site_ranks.==0.0] .= n_locs + 1
    return (raw=C_cover, seed_log=Yseed, fog_log=Yfog, shade_log=Yshade, site_ranks=site_ranks, bleaching_mortality=bleaching_mort, coral_dhw_log=collated_dhw_tol_log)
end

function cyclone_mortality!(coral_cover, coral_params, cyclone_mortality)::Nothing
    # Small class coral mortality
    coral_deaths_small = coral_cover[coral_params.small, :] .* cyclone_mortality
    coral_cover[coral_params.small, :] -= coral_deaths_small

    # Mid class coral mortality
    coral_mid = hcat(collect(Iterators.partition(coral_params.mid, 4))...)
    for i in size(coral_mid, 1)
        coral_deaths_mid = coral_cover[coral_mid[i, :], :] .* cyclone_mortality
        coral_cover[coral_mid[i, :], :] -= coral_deaths_mid
    end

    # Large class coral mortality
    coral_deaths_large = coral_cover[coral_params.large, :] .* cyclone_mortality
    coral_cover[coral_params.large, :] -= coral_deaths_large
    return nothing
end
