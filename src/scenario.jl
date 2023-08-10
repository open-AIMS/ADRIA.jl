"""Scenario running functions"""

using ADRIA.metrics: relative_cover, relative_loc_taxa_cover, total_absolute_cover, absolute_shelter_volume, relative_shelter_volume
using ADRIA.metrics: relative_juveniles, relative_taxa_cover, juvenile_indicator
using ADRIA.metrics: coral_evenness

"""
    setup_cache(domain::Domain)::NamedTuple

Establish tuple of matrices/vectors for use as reusable data stores to avoid repeated memory allocations.
"""
function setup_cache(domain::Domain)::NamedTuple

    # sim constants
    n_sites::Int64 = domain.coral_growth.n_sites
    n_species::Int64 = domain.coral_growth.n_species
    n_groups::Int64 = domain.coral_growth.n_groups

    init_cov = domain.init_coral_cover
    cache = (
        sf=zeros(n_groups, n_sites),  # stressed fecundity
        fec_all=zeros(size(init_cov)...),  # all fecundity
        fec_scope=zeros(n_groups, n_sites),  # fecundity scope
        dhw_step=zeros(n_sites),  # DHW each time step
        cov_tmp=zeros(size(init_cov)...),  # Cover for previous timestep
        depth_coeff=zeros(n_sites),  # store for depth coefficient
        site_area=Matrix{Float64}(site_area(domain)'),  # site areas
        TP_data=Matrix{Float64}(domain.TP_data),  # transition probabilities
        waves=zeros(length(timesteps(domain)), n_species, n_sites)
    )

    return cache
end


"""
    run_scenarios(param_df::DataFrame, domain::Domain, RCP::String; show_progress=true, remove_workers=true)
    run_scenarios(param_df::DataFrame, domain::Domain, RCP::Vector{String}; show_progress=true, remove_workers=true)

Run scenarios defined by the parameter table storing results to disk.
Scenarios are run in parallel where the number of scenarios > 256.

# Notes
- Returned `domain` holds scenario invoke time used as unique result set identifier.
- If multiple RCPs are specified, this method will temporarily use double the disk space
  to consolidate results into a single ResultSet.

# Examples
```julia-repl
...
julia> rs_45 = ADRIA.run_scenarios(p_df, dom, "45")
julia> rs_45_60 = ADRIA.run_scenarios(p_df, dom, ["45", "60"])
```

# Arguments
- param_df : DataFrame of scenarios to run
- domain : Domain, to run scenarios with
- RCP : ID or list of of RCP(s) to run scenarios under.
- show_progress : Display progress
- remove_workers : If running in parallel, removes workers after completion

# Returns
ResultSet
"""
function run_scenarios(param_df::DataFrame, domain::Domain, RCP::String; show_progress=true, remove_workers=true)::ResultSet
    return run_scenarios(param_df, domain, [RCP]; show_progress, remove_workers)
end
function run_scenarios(param_df::DataFrame, domain::Domain, RCP::Vector{String}; show_progress=true, remove_workers=true)::ResultSet
    # Initialize ADRIA configuration options
    setup()

    # Sort RCPs so the dataframe order match the output filepath
    RCP = sort(RCP)

    @info "Running $(nrow(param_df)) scenarios over $(length(RCP)) RCPs: $RCP"

    # Cross product between rcps and param_df to have every row of param_df for each rcp
    rcps_df = DataFrame(RCP=parse.(Int64, RCP))
    scenarios_df = crossjoin(param_df, rcps_df)
    sort!(scenarios_df, :RCP)

    @info "Setting up Result Set"
    domain, data_store = ADRIA.setup_result_store!(domain, scenarios_df)

    # Convert DataFrame to named matrix for faster iteration
    scenarios_matrix = NamedDimsArray(
        Matrix(scenarios_df);
        scenarios=1:nrow(scenarios_df),
        factors=names(scenarios_df)
    )

    parallel = (nrow(param_df) >= 256) && (parse(Bool, ENV["ADRIA_DEBUG"]) == false)
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
            run_msg = "Running $(nrow(param_df)) scenarios for RCP $rcp"

            # Switch RCPs so correct data is loaded
            target_rows = findall(scenarios_matrix("RCP") .== parse(Float64, rcp))
            if show_progress
                @showprogress run_msg 4 pmap(func, CachingPool(workers()), zip(target_rows, eachrow(scenarios_matrix[target_rows, :]), Iterators.repeated(domain, size(scenarios_matrix, 1))))
            else
                pmap(func, CachingPool(workers()), zip(target_rows, eachrow(scenarios_matrix[target_rows, :]), Iterators.repeated(domain, size(scenarios_matrix, 1))))
            end
        end
    else
        # Cache to reuse during scenario runs
        cache = setup_cache(domain)

        # Define local helper
        func = dfx -> run_scenario(dfx..., domain, data_store, cache)

        for rcp in RCP
            run_msg = "Running $(nrow(param_df)) scenarios for RCP $rcp"

            # Switch RCPs so correct data is loaded
            domain = switch_RCPs!(domain, rcp)
            target_rows = findall(scenarios_matrix("RCP") .== parse(Float64, rcp))
            if show_progress
                @showprogress run_msg 4 map(func, zip(target_rows, eachrow(scenarios_matrix[target_rows, :])))
            else
                map(func, zip(target_rows, eachrow(scenarios_matrix[target_rows, :])))
            end
        end
    end

    if parallel && remove_workers
        _remove_workers()
    end

    return load_results(_result_location(domain, RCP))
end

"""
    run_scenario(idx::Int64, param_set::Union{AbstractVector, DataFrameRow}, domain::Domain, data_store::NamedTuple, cache::NamedTuple)::Nothing
    run_scenario(idx::Int64, param_set::Union{AbstractVector, DataFrameRow}, domain::Domain, data_store::NamedTuple)::Nothing
    run_scenario(param_set::Union{AbstractVector, DataFrameRow}, domain::Domain, cache::NamedTuple)::NamedTuple
    run_scenario(param_set::NamedTuple, domain::Domain)::NamedTuple

Run individual scenarios for a given domain, saving results to a Zarr data store.
Results are stored in Zarr format at a pre-configured location.
Sets up a new `cache` if not provided.

# Notes
Logs of site ranks only store the mean site rankings over all environmental scenarios.
This is to reduce the volume of data stored.

# Returns
Nothing
"""
function run_scenario(idx::Int64, param_set::Union{AbstractVector,DataFrameRow}, domain::Domain,
    data_store::NamedTuple, cache::NamedTuple)::Nothing

    coral_params = to_coral_spec(param_set)
    if domain.RCP == ""
        local rcp
        try
            rcp = param_set("RCP")  # Try extracting from NamedDimsArray
        catch err
            if !(err isa MethodError)
                rethrow(err)
            end

            rcp = param_set.RCP  # Extract from dataframe
        end

        domain = switch_RCPs!(domain, string(Int64(rcp)))
    end

    result_set = run_model(domain, param_set, coral_params, cache)

    # Capture results to disk
    # Set values below threshold to 0 to save space
    threshold = parse(Float32, ENV["ADRIA_THRESHOLD"])

    rs_raw = result_set.raw
    vals = total_absolute_cover(rs_raw, site_area(domain))
    vals[vals.<threshold] .= 0.0
    data_store.total_absolute_cover[:, :, idx] .= vals

    vals .= absolute_shelter_volume(rs_raw, site_area(domain), param_set)
    vals[vals.<threshold] .= 0.0
    data_store.absolute_shelter_volume[:, :, idx] .= vals

    vals .= relative_shelter_volume(rs_raw, site_area(domain), site_k_area(domain), param_set)
    vals[vals.<threshold] .= 0.0
    data_store.relative_shelter_volume[:, :, idx] .= vals

    coral_spec::DataFrame = to_coral_spec(param_set)
    vals .= relative_juveniles(rs_raw, coral_spec)
    vals[vals.<threshold] .= 0.0
    data_store.relative_juveniles[:, :, idx] .= vals

    vals .= juvenile_indicator(rs_raw, coral_spec, site_area(domain), site_k_area(domain))
    vals[vals.<threshold] .= 0.0
    data_store.juvenile_indicator[:, :, idx] .= vals

    vals = relative_taxa_cover(rs_raw, site_k_area(domain), site_area(domain))
    vals[vals.<threshold] .= 0.0
    data_store.relative_taxa_cover[:, :, idx] .= vals

    vals = relative_loc_taxa_cover(rs_raw, site_k_area(domain), site_area(domain))
    vals = coral_evenness(NamedDims.unname(vals))
    vals[vals.<threshold] .= 0.0
    data_store.coral_evenness[:, :, idx] .= vals

    # Store raw results if no metrics specified
    # if length(metrics) == 0
    #     data_store.raw[:, :, :, idx] .= r.raw
    # end

    # Store logs
    tf = size(domain.dhw_scens, 1)  # time frame
    tmp_site_ranks = zeros(Float32, tf, nrow(domain.site_data), 2)
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
            tmp_site_ranks[:, :, :] .= vals
        elseif k == :coral_dhw_log
            getfield(data_store, k)[:, :, :, :, idx] .= vals
        else
            getfield(data_store, k)[:, :, idx] .= vals
        end
    end

    if !isnothing(data_store.site_ranks)
        # Squash site ranks down to average rankings over environmental repeats
        data_store.site_ranks[:, :, :, idx] .= tmp_site_ranks
    end

    # There seems to be an issue when using `pmap` for distributed tasks, where the
    # garbage collector is not run leading to excessive memory use (Julia v1.9.2).
    # Force garbage collection if > 95% of total memory is found to be in use.
    # if (Sys.free_memory() ./ Sys.total_memory()) .< 0.05
    #     GC.gc()
    # end

    return nothing
end
function run_scenario(idx::Int64, param_set::Union{AbstractVector,DataFrameRow}, domain::Domain, data_store::NamedTuple)::Nothing
    return run_scenario(idx, param_set, domain, data_store, setup_cache(domain))
end
function run_scenario(param_set::Union{AbstractVector,DataFrameRow}, domain::Domain)::NamedTuple
    coral_params = to_coral_spec(param_set)
    return run_model(domain, param_set, coral_params, setup_cache(domain))
end
function run_scenario(param_set::Union{AbstractVector,DataFrameRow}, domain::Domain, RCP::String)::NamedTuple
    domain = switch_RCPs!(domain, RCP)
    return run_scenario(param_set, domain)
end


"""
    run_model(domain::Domain, param_set::Union{NamedTuple,DataFrameRow}, corals::DataFrame, cache::NamedTuple)::NamedTuple

Core scenario running function.

# Notes
Only the mean site rankings are kept

# Returns
NamedTuple of collated results
"""
function run_model(domain::Domain, param_set::DataFrameRow, corals::DataFrame, cache::NamedTuple)::NamedTuple
    ps = NamedDimsArray(Vector(param_set), factors=names(param_set))
    return run_model(domain, ps, corals, cache)
end
function run_model(domain::Domain, param_set::NamedDimsArray, corals::DataFrame, cache::NamedTuple)::NamedTuple
    p = domain.coral_growth.ode_p

    # Set random seed using intervention values
    # TODO: More robust way of getting intervention/criteria values
    rnd_seed_val::Int64 = floor(Int64, sum(param_set(!=("RCP"))))  # select everything except RCP
    Random.seed!(rnd_seed_val)

    dhw_idx::Int64 = Int64(param_set("dhw_scenario"))
    wave_idx::Int64 = Int64(param_set("wave_scenario"))

    dhw_scen::Matrix{Float64} = copy(domain.dhw_scens[:, :, dhw_idx])

    tspan::Tuple = (0.0, 1.0)
    solver::Euler = Euler()

    MCDA_approach::Int64 = param_set("guided")

    # Sim constants
    sim_params = domain.sim_constants
    tf::Int64 = size(dhw_scen, 1)
    n_site_int::Int64 = sim_params.n_site_int
    n_sites::Int64 = domain.coral_growth.n_sites
    n_species::Int64 = domain.coral_growth.n_species
    n_groups::Int64 = domain.coral_growth.n_groups

    # years to start seeding/shading
    seed_start_year::Int64 = param_set("seed_year_start")
    shade_start_year::Int64 = param_set("shade_year_start")

    fogging::Real = param_set("fogging")  # percent reduction in bleaching mortality through fogging
    srm::Real = param_set("SRM")  # DHW equivalents reduced by some shading mechanism
    seed_years::Int64 = param_set("seed_years")  # number of years to seed
    shade_years::Int64 = param_set("shade_years")  # number of years to shade

    total_loc_area::Matrix{Float64} = cache.site_area
    fec_params_per_m²::Vector{Float64} = corals.fecundity  # number of larvae produced per m²

    # Caches
    TP_data = cache.TP_data
    sf = cache.sf
    fec_all = cache.fec_all
    fec_scope = cache.fec_scope
    dhw_t = cache.dhw_step
    Y_pstep = cache.cov_tmp
    depth_coeff = cache.depth_coeff

    site_data = domain.site_data
    depth_coeff .= depth_coefficient.(site_data.depth_med)

    Y_cover::Array{Float64,3} = zeros(tf, n_species, n_sites)  # Coral cover relative to total site area
    Y_cover[1, :, :] .= domain.init_coral_cover
    ode_u = zeros(n_species, n_sites)
    cover_tmp = p.cover  # pre-allocated matrix used to avoid memory allocations
    max_cover = site_k(domain)  # Max coral cover at each site (0 - 1).

    # Proportionally adjust initial cover (handles inappropriate initial conditions)
    proportional_adjustment!(@view(Y_cover[1, :, :]), cover_tmp, max_cover)

    site_ranks = SparseArray(zeros(tf, n_sites, 2))  # log seeding/fogging/shading ranks
    Yshade = SparseArray(spzeros(tf, n_sites))
    Yfog = SparseArray(spzeros(tf, n_sites))
    Yseed = SparseArray(zeros(tf, 3, n_sites))  # 3 = the number of seeded coral types

    # Prep scenario-specific flags/values
    # Intervention strategy: < 0 is no intervention, 0 is random location selection, > 0 is guided
    is_guided = param_set("guided") > 0

    # Decisions should place more weight on environmental conditions
    # closer to the decision point
    α = 0.95
    decay = α .^ (1:Int64(param_set("plan_horizon"))+1)

    # Years at which to reassess seeding site selection
    seed_decision_years = fill(false, tf)
    shade_decision_years = fill(false, tf)

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

    seed_locs::Vector{Int64} = zeros(Int, n_site_int)
    shade_locs::Vector{Int64} = zeros(Int, n_site_int)

    # Set other params for ODE
    p.r .= corals.growth_rate  # Assumed growth_rate
    p.mb .= corals.mb_rate  # background mortality
    @set! p.k = max_cover  # max coral cover

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

    # Defaults to considering all sites if depth cannot be considered.
    depth_priority = collect(1:n_sites)

    # Calculate total area to seed respecting tolerance for minimum available space to still
    # seed at a site
    area_to_seed = sum(seeded_area)

    # Filter out sites outside of desired depth range
    if .!all(site_data.depth_med .== 0)
        max_depth::Float64 = param_set("depth_min") + param_set("depth_offset")
        depth_criteria::BitArray{1} = (site_data.depth_med .>= param_set("depth_min")) .& (site_data.depth_med .<= max_depth)

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
        mcda_vars = DMCDA_vars(domain, param_set, depth_priority, sum(Y_cover[1, :, :], dims=1), area_to_seed)

        # Number of time steps in environmental layers to look ahead when making decisions
        plan_horizon::Int64 = Int64(param_set("plan_horizon"))
    end

    # Set up distributions for natural adaptation/heritability
    c_dist_t_1::Matrix{Distribution} = repeat(
        truncated.(
            Normal.(corals.dist_mean, corals.dist_std),
            0.0,
            corals.dist_mean .+ HEAT_UB
        ),
        1,
        n_sites
    )
    c_dist_t = copy(c_dist_t_1)
    dhw_tol_log = permutedims(repeat(c_dist_t_1, 1, 1, tf), [3, 1, 2])

    # Cache for proportional mortality and coral population increases
    bleaching_mort = zeros(tf, n_species, n_sites)

    #### End coral constants

    ## Update ecological parameters based on intervention option

    # Treat as enhancement from mean of "natural" DHW tolerance
    a_adapt[a_adapt.>0.0] .+= corals.dist_mean[a_adapt.>0.0]

    # TODO: Better conversion of Ub to wave mortality
    #       Currently scaling significant wave height by its max to non-dimensionalize values
    wave_scen::Matrix{Float64} = copy(domain.wave_scens[:, :, wave_idx])
    wave_scen .= wave_scen ./ maximum(wave_scen)

    # Pre-calculate proportion of survivers from wave stress
    # Sw_t = wave_damage!(cache.waves, wave_scen, corals.wavemort90, n_species)

    site_k_prop = max_cover'
    absolute_k_area = site_k_area(domain)'  # max possible coral area in m^2
    growth::ODEProblem = ODEProblem{true}(growthODE, ode_u, tspan, p)
    tmp::Matrix{Float64} = zeros(size(Y_cover[1, :, :]))  # temporary array to hold intermediate covers

    # basal_area_per_settler is the area in m^2 of a size class one coral
    basal_area_per_settler = colony_mean_area(corals.mean_colony_diameter_m[corals.class_id.==1])
    @inbounds for tstep::Int64 in 2:tf
        p_step::Int64 = tstep - 1
        Y_pstep[:, :] .= Y_cover[p_step, :, :]

        # Calculates scope for coral fedundity for each size class and at each site.
        fecundity_scope!(fec_scope, fec_all, fec_params_per_m², Y_pstep, total_loc_area)

        site_coral_cover = sum(Y_pstep, dims=1)  # dims: nsites * 1
        leftover_space_prop = relative_leftover_space(site_k_prop, site_coral_cover)
        leftover_space_m² = leftover_space_prop .* total_loc_area

        # Recruitment represents additional cover, relative to total site area
        # Gets used in ODE
        p.rec .= settler_cover(fec_scope, TP_data, leftover_space_prop,
            sim_params.max_settler_density, sim_params.max_larval_density, basal_area_per_settler)

        in_shade_years = (shade_start_year <= tstep) && (tstep <= (shade_start_year + shade_years - 1))
        in_seed_years = (seed_start_year <= tstep) && (tstep <= (seed_start_year + seed_years - 1))

        @views dhw_t .= dhw_scen[tstep, :]  # subset of DHW for given timestep

        # Apply regional cooling effect before selecting locations to seed
        if (srm > 0.0) && in_shade_years
            Yshade[tstep, :] .= srm

            # Apply reduction in DHW due to SRM
            dhw_t .= max.(0.0, dhw_t .- srm)
        end

        if is_guided && (in_seed_years || in_shade_years)
            # Update dMCDA values

            # Determine subset of data to select data for planning horizon
            horizon::UnitRange{Int64} = tstep:min(tstep+plan_horizon, tf)
            d_s::UnitRange{Int64} = 1:length(horizon)

            # Put more weight on projected conditions closer to the decision point
            @views env_horizon = decay[d_s] .* dhw_scen[horizon, :]
            mcda_vars.heat_stress_prob .= vec((mean(env_horizon, dims=1) .+ std(env_horizon, dims=1)) .* 0.5)

            @views env_horizon = decay[d_s] .* wave_scen[horizon, :]
            mcda_vars.dam_prob .= vec((mean(env_horizon, dims=1) .+ std(env_horizon, dims=1)) .* 0.5)
        end
        if is_guided && (in_seed_years || in_shade_years)
            mcda_vars.sum_cover .= site_coral_cover

            # Determine connectivity strength
            # Account for cases where there is no coral cover
            in_conn, out_conn, strong_pred = connectivity_strength(domain.TP_data .* site_k_area(domain), vec(site_coral_cover))
            (seed_locs, shade_locs, rankings) = guided_site_selection(mcda_vars, MCDA_approach,
                seed_decision_years[tstep], shade_decision_years[tstep],
                seed_locs, shade_locs, rankings, in_conn[mcda_vars.site_ids], out_conn[mcda_vars.site_ids], strong_pred[mcda_vars.site_ids])

            # Log site ranks
            # First col only holds site index ids so skip (with 2:end)
            site_ranks[tstep, rankings[:, 1], :] = rankings[:, 2:end]
        elseif seed_corals && (in_seed_years || in_shade_years)
            # Unguided deployment, seed/shade corals anywhere, so long as available space > 0
            seed_locs, shade_locs = unguided_site_selection(seed_locs, shade_locs,
                seed_decision_years[tstep], shade_decision_years[tstep],
                n_site_int, vec(leftover_space_m²), depth_priority)

            site_ranks[tstep, seed_locs, 1] .= 1.0
            site_ranks[tstep, shade_locs, 2] .= 1.0
        end

        has_shade_sites::Bool = !all(shade_locs .== 0)
        has_seed_sites::Bool = !all(seed_locs .== 0)

        # Fog selected locations
        if (fogging > 0.0) && in_shade_years && has_shade_sites
            fog_locations!(@view(Yfog[tstep, :]), shade_locs, dhw_t, fogging)
        end

        # Calculate and apply bleaching mortality
        #    This: `dhw_t .* (1.0 .- wave_scen[tstep, :])`
        #    attempts to account for the cooling effect of storms / high wave activity
        # `wave_scen` is normalized to the maximum value found for the given wave scenario
        # so what causes 100% mortality can differ between runs.
        bleaching_mortality!(Y_pstep, dhw_t .* (1.0 .- @view(wave_scen[tstep, :])), depth_coeff, c_dist_t, c_dist_t1, @view(bleaching_mort[tstep, :, :]))

        # Apply seeding
        if seed_corals && in_seed_years && has_seed_sites
            # Seed each selected site
            seed_corals!(Y_pstep, vec(total_loc_area), vec(leftover_space_m²),
                seed_locs, seeded_area, seed_sc, a_adapt, @view(Yseed[tstep, :, :]),
                corals.dist_std, c_dist_t)
        end

        # Note: ODE is run relative to `k` area, but values are otherwise recorded
        #       in relative to absolute area.
        # Update initial condition
        @. tmp = (Y_pstep * total_loc_area) / absolute_k_area
        replace!(tmp, Inf => 0.0, NaN => 0.0)
        growth.u0 .= tmp

        # X is cover relative to `k` (max. carrying capacity)
        # So we subtract from 1.0 to get leftover/available space, relative to `k`
        p.sXr .= max.(1.0 .- sum(tmp, dims=1), 0.0) .* tmp .* p.r  # leftover space * current cover * growth_rate
        p.X_mb .= tmp .* p.mb    # current cover * background mortality

        sol::ODESolution = solve(growth, solver, save_everystep=false, save_start=false,
            alg_hints=[:nonstiff], adaptive=false, dt=0.5)

        # Ensure values are ∈ [0, 1]
        @views Y_cover[tstep, :, :] .= clamp.(sol.u[end] .* absolute_k_area ./ total_loc_area, 0.0, 1.0)

        if tstep <= tf
            adjust_DHW_distribution!(@view(Y_cover[tstep-1:tstep, :, :]), n_groups, c_dist_t_1, c_dist_t,
                p.r, corals.dist_std, param_set("heritability"))

            dhw_tol_log[tstep, :, :] .= c_dist_t
            c_dist_t_1 .= c_dist_t
        end
    end

    # Could collate critical DHW threshold log for corals to reduce disk space...
    # dhw_tol_mean = dropdims(mean(mean.(dhw_tol_log), dims=3), dims=3)
    # dhw_tol_mean_std = dropdims(mean(std.(dhw_tol_log), dims=3), dims=3)
    # collated_dhw_tol_log = NamedDimsArray(cat(dhw_tol_mean, dhw_tol_mean_std, dims=3),
    #     timesteps=1:tf, species=corals.coral_id, stat=[:mean, :stdev])
    dhw_tol_mean_log = mean.(dhw_tol_log)
    dhw_tol_mean_std_log = std.(dhw_tol_log)

    collated_dhw_tol_log = NamedDimsArray(cat(dhw_tol_mean_log, dhw_tol_mean_std_log, dims=ndims(dhw_tol_log) + 1),
        timesteps=1:tf, species=corals.coral_id, sites=1:n_sites, stat=[:mean, :stdev])

    # Avoid placing importance on sites that were not considered
    # (lower values are higher importance)
    site_ranks[site_ranks.==0.0] .= n_sites + 1
    return (raw=Y_cover, seed_log=Yseed, fog_log=Yfog, shade_log=Yshade, site_ranks=site_ranks, bleaching_mortality=bleaching_mort, coral_dhw_log=collated_dhw_tol_log)
end
