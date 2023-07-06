"""Scenario running functions"""

using ADRIA.metrics: relative_cover, total_absolute_cover, absolute_shelter_volume, relative_shelter_volume
using ADRIA.metrics: relative_juveniles, relative_taxa_cover, juvenile_indicator


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
        prop_loss=zeros(n_species, n_sites),  # proportional loss
        Sbl=zeros(n_species, n_sites),   # bleaching survivors
        dhw_step=zeros(n_sites),  # DHW each time step
        cov_tmp=zeros(size(init_cov)...),  # Cover for previous timestep
        felt_dhw=zeros(size(init_cov)...),  # Store for felt DHW (DHW after reductions)
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

    # Setup cache to reuse for each scenario run
    cache = setup_cache(domain)

    parallel = (nrow(param_df) >= 4096) && (parse(Bool, ENV["ADRIA_DEBUG"]) == false)
    if parallel && nworkers() == 1
        @info "Setting up parallel processing..."
        spinup_time = @elapsed begin
            _setup_workers()
            sleep(2)  # wait a bit while workers spin-up

            # load ADRIA on workers and define helper function
            @everywhere @eval using ADRIA
            @everywhere @eval func = (dfx) -> run_scenario(dfx..., domain, data_store, cache)
        end

        @info "Time taken to spin up workers: $(spinup_time) seconds"

        # Define number of scenarios to run before returning results to main
        # https://discourse.julialang.org/t/parallelism-understanding-pmap-and-the-batch-size-parameter/15604/2
        # https://techytok.com/lesson-parallel-computing/
        b_size = 4
    else
        b_size = 1
    end

    # Define local helper
    func = (dfx) -> run_scenario(dfx..., domain, data_store, cache)

    for rcp in RCP
        run_msg = "Running $(nrow(param_df)) scenarios for RCP $rcp"

        # Switch RCPs so correct data is loaded
        domain = switch_RCPs!(domain, rcp)
        target_rows = findall(scenarios_matrix("RCP") .== parse(Float64, rcp))
        if show_progress
            @showprogress run_msg 4 pmap(func, zip(target_rows, eachrow(scenarios_matrix[target_rows, :])), batch_size=b_size)
        else
            pmap(func, zip(target_rows, eachrow(scenarios_matrix[target_rows, :])), batch_size=b_size)
        end
    end

    if parallel && remove_workers
        _remove_workers()
    end

    return load_results(_result_location(domain, RCP))
end

"""
    run_scenario(idx::Int64, param_set::Union{AbstractVector, DataFrameRow}, domain::Domain, data_store::NamedTuple, cache::NamedTuple)::NamedTuple
    run_scenario(idx::Int64, param_set::Union{AbstractVector, DataFrameRow}, domain::Domain, data_store::NamedTuple)::NamedTuple
    run_scenario(param_set::Union{AbstractVector, DataFrameRow}, domain::Domain, cache::NamedTuple)::NamedTuple
    run_scenario(param_set::NamedTuple, domain::Domain)::NamedTuple

Run individual scenarios for a given domain, saving results to a Zarr data store.
Results are stored in Zarr format at a pre-configured location.
Sets up a new `cache` if not provided.

# Notes
Logs of site ranks only store the mean site rankings over all environmental scenarios.
This is to reduce the volume of data stored.
"""
function run_scenario(idx::Int64, param_set::Union{AbstractVector,DataFrameRow}, domain::Domain,
    data_store::NamedTuple, cache::NamedTuple)

    result_set = run_scenario(param_set, domain, cache)

    # Capture results to disk
    # Set values below threshold to 0 to save space
    tf = size(domain.dhw_scens, 1)
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

    vals = relative_taxa_cover(rs_raw)
    vals[vals.<threshold] .= 0.0
    data_store.relative_taxa_cover[:, :, idx] .= vals

    # Store raw results if no metrics specified
    # if length(metrics) == 0
    #     data_store.raw[:, :, :, idx] .= r.raw
    # end

    # Store logs
    tmp_site_ranks = zeros(Float32, tf, nrow(domain.site_data), 2)
    c_dim = Base.ndims(result_set.raw) + 1
    log_stores = (:site_ranks, :seed_log, :fog_log, :shade_log)
    for k in log_stores
        if k == :seed_log || k == :site_ranks
            concat_dim = c_dim
        else
            concat_dim = c_dim - 1
        end

        vals = getfield(result_set, k)
        vals[vals.<threshold] .= 0.0

        if k == :seed_log
            getfield(data_store, k)[:, :, :, idx] .= vals
        elseif k == :site_ranks
            tmp_site_ranks[:, :, :] .= vals
        else
            getfield(data_store, k)[:, :, idx] .= vals
        end
    end

    if !isnothing(data_store.site_ranks)
        # Squash site ranks down to average rankings over environmental repeats
        data_store.site_ranks[:, :, :, idx] .= tmp_site_ranks
    end
end
function run_scenario(idx::Int64, param_set::Union{AbstractVector,DataFrameRow}, domain::Domain, data_store::NamedTuple)
    cache = setup_cache(domain)
    return run_scenario(idx, param_set, domain, data_store, cache)
end
function run_scenario(param_set::Union{AbstractVector,DataFrameRow}, domain::Domain, cache::NamedTuple)
    # Extract coral only parameters
    coral_params = to_coral_spec(param_set)
    if domain.RCP == ""
        throw("No RCP set for domain. Use `ADRIA.switch_RCPs!() to select an RCP to run scenarios with.")
    end

    return run_model(domain, param_set, coral_params, cache)
end
function run_scenario(param_set::Union{AbstractVector,DataFrameRow}, domain::Domain)::NamedTuple
    cache = setup_cache(domain)
    return run_scenario(param_set, domain, cache)
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
    sim_params = domain.sim_constants
    site_data = domain.site_data
    p = domain.coral_growth.ode_p

    # Set random seed using intervention values
    # TODO: More robust way of getting intervention/criteria values
    rnd_seed_val::Int64 = floor(Int64, sum(param_set(!=("RCP"))))  # select everything except RCP
    Random.seed!(rnd_seed_val)

    dhw_idx::Int64 = Int64(param_set("dhw_scenario"))
    wave_idx::Int64 = Int64(param_set("wave_scenario"))

    dhw_scen::Matrix{Float64} = @view domain.dhw_scens[:, :, dhw_idx]

    # TODO: Better conversion of Ub to wave mortality
    #       Currently scaling significant wave height by its max to non-dimensionalize values
    wave_scen::Matrix{Float64} = Matrix{Float64}(domain.wave_scens[:, :, wave_idx]) ./ maximum(domain.wave_scens[:, :, wave_idx])

    tspan::Tuple = (0.0, 1.0)
    solver::Euler = Euler()

    MCDA_approach::Int64 = param_set("guided")

    # sim constants
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

    total_site_area::Array{Float64,2} = cache.site_area
    fec_params_per_m²::Vector{Float64} = corals.fecundity  # number of larvae produced per m²

    # Caches
    TP_data = cache.TP_data
    sf = cache.sf
    fec_all = cache.fec_all
    fec_scope = cache.fec_scope
    prop_loss = cache.prop_loss
    Sbl = cache.Sbl
    dhw_t = cache.dhw_step
    Y_pstep = cache.cov_tmp
    felt_dhw = cache.felt_dhw
    depth_coeff = cache.depth_coeff

    depth_coeff .= depth_coefficient.(site_data.depth_med)

    Y_cover::Array{Float64,3} = zeros(tf, n_species, n_sites)  # Coral cover relative to total site area
    Y_cover[1, :, :] .= domain.init_coral_cover
    ode_u = zeros(n_species, n_sites)
    cover_tmp = p.cover  # pre-allocated matrix used to avoid memory allocations

    site_ranks = SparseArray(zeros(tf, n_sites, 2)) # log seeding/fogging/shading ranks
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
    seed_decision_years = repeat([false], tf)
    shade_decision_years = repeat([false], tf)

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

    prefseedsites::Vector{Int64} = zeros(Int, n_site_int)
    prefshadesites::Vector{Int64} = zeros(Int, n_site_int)

    # Max coral cover at each site (0 - 1).
    max_cover = site_k(domain)

    # Set other params for ODE
    p.r .= corals.growth_rate  # Assumed growth_rate

    p.mb .= corals.mb_rate  # background mortality
    @set! p.k = max_cover  # max coral cover

    # Proportionally adjust initial cover (handles inappropriate initial conditions)
    proportional_adjustment!(Y_cover[1, :, :], cover_tmp, max_cover)

    # Define taxa and size class to seed, and identify their factor names
    taxa_to_seed = [2, 3, 5]
    target_class_id::BitArray = corals.class_id .== 2  # seed second smallest size class
    taxa_names = param_set.factors[occursin.("N_seed_", param_set.factors)]

    # Identify taxa and size class to be seeded
    seed_sc = (corals.taxa_id .∈ [taxa_to_seed]) .& target_class_id

    # Extract colony areas for sites selected in m^2 and add adaptation values
    seeded_area = colony_mean_area(corals.mean_colony_diameter_m[seed_sc]) .* param_set(taxa_names)

    # Set up assisted adaptation values
    a_adapt = zeros(n_species)
    a_adapt[seed_sc] .= param_set("a_adapt")

    # Flag indicating whether to seed or not to seed
    seed_corals = any(param_set(taxa_names) .> 0.0)

    bleaching_sensitivity = corals.bleaching_sensitivity

    # Defaults to considering all sites if depth cannot be considered.
    depth_priority = collect(1:nrow(site_data))

    # Calculate total area to seed respecting tolerance for minimum available space to still
    # seed at a site
    area_to_seed = sum(seeded_area)

    # Filter out sites outside of desired depth range
    if .!all(site_data.depth_med .== 0)
        max_depth::Float64 = param_set("depth_min") + param_set("depth_offset")
        depth_criteria::BitArray{1} = (site_data.depth_med .>= param_set("depth_min")) .& (site_data.depth_med .<= max_depth)

        # TODO: Include this change in MATLAB version as well
        if any(depth_criteria .> 0)
            # If sites can be filtered based on depth, do so.
            depth_priority = depth_priority[depth_criteria]
        else
            # Otherwise if no sites can be filtered, remove depth as a criterion.
            @warn "No sites within provided depth range of $(param_set("depth_min")) - $(max_depth) meters. Considering all sites."
        end
    end

    if is_guided
        # pre-allocate rankings
        rankings = [depth_priority zeros(Int, length(depth_priority)) zeros(Int, length(depth_priority))]

        # Prep site selection
        mcda_vars = DMCDA_vars(domain, param_set, depth_priority, sum(Y_cover[1, :, :], dims=1), area_to_seed)
    end

    # Set up distributions for natural adaptation/heritability
    c_dist_t::Matrix{Distribution} = repeat(
        TruncatedNormal.(corals.dist_mean, corals.dist_std, 0.0, corals.dist_mean .+ HEAT_UB),
        1, n_sites
    )
    c_dist_t1 = copy(c_dist_t)

    # Identify juvenile classes
    # juveniles = corals.class_id .∈ [[1, 2]]
    # deployed_size_class = corals.class_id .∈ 2

    # Cache for proportional mortality and coral population increases
    bleaching_mort = zeros(tf, n_species, n_sites)
    c_increase = zeros(n_groups)

    #### End coral constants

    ## Update ecological parameters based on intervention option

    # Treat as enhancement from mean of "natural" DHW tolerance
    a_adapt[a_adapt.>0.0] .+= corals.dist_mean[a_adapt.>0.0]

    # Level of natural coral adaptation
    n_adapt = param_set("n_adapt")

    ## Extract other parameters
    LPdhwcoeff = sim_params.LPdhwcoeff # shape parameters relating dhw affecting cover to larval production
    DHWmaxtot = sim_params.DHWmaxtot # max assumed DHW for all scenarios.  Will be obsolete when we move to new, shared inputs for DHW projections
    LPDprm2 = sim_params.LPDprm2 # parameter offsetting LPD curve

    # Wave stress
    Sw_t::Array{Float64,3} = cache.waves
    wavemort90::Vector{Float64} = corals.wavemort90::Vector{Float64}  # 90th percentile wave mortality

    for sp::Int64 in 1:n_species
        @views Sw_t[:, sp, :] .= wavemort90[sp] .* wave_scen
    end

    clamp!(Sw_t, 0.0, 1.0)

    # Wave damage survival
    Sw_t .= 1.0 .- Sw_t

    site_k_prop = max_cover'
    absolute_k_area = site_k_area(domain)'  # max possible coral area in m^2
    growth::ODEProblem = ODEProblem{true}(growthODE, ode_u, tspan, p)
    tmp::Matrix{Float64} = zeros(size(Y_cover[1, :, :]))  # temporary array to hold intermediate covers

    env_horizon = zeros(Int64(param_set("plan_horizon") + 1), n_sites)  # temporary cache for planning horizon

    # basal_area_per_settler is the area in m^2 of a size class one coral
    basal_area_per_settler = colony_mean_area(corals.mean_colony_diameter_m[corals.class_id.==1])
    @inbounds for tstep::Int64 in 2:tf
        p_step::Int64 = tstep - 1
        Y_pstep[:, :] .= Y_cover[p_step, :, :]

        sf .= stressed_fecundity(tstep, a_adapt, n_adapt, dhw_scen[p_step, :],
            LPdhwcoeff, DHWmaxtot, LPDprm2, n_groups)

        # Calculates scope for coral fedundity for each size class and at each site.
        fecundity_scope!(fec_scope, fec_all, fec_params_per_m², Y_pstep, total_site_area)

        site_coral_cover = sum(Y_pstep, dims=1)  # dims: nsites * 1
        leftover_space_prop = relative_leftover_space(site_k_prop, site_coral_cover)
        leftover_space_m² = leftover_space_prop .* total_site_area

        # Recruitment represents additional cover, relative to total site area
        # Gets used in ODE
        p.rec .= settler_cover(fec_scope, sf, TP_data, leftover_space_prop,
            sim_params.max_settler_density, sim_params.max_larval_density, basal_area_per_settler)

        in_shade_years = (shade_start_year <= tstep) && (tstep <= (shade_start_year + shade_years - 1))
        in_seed_years = ((seed_start_year <= tstep) && (tstep <= (seed_start_year + seed_years - 1)))

        @views dhw_t .= dhw_scen[tstep, :]  # subset of DHW for given timestep

        # Apply regional cooling effect before selecting locations to seed
        if (srm > 0.0) && in_shade_years
            Yshade[tstep, :] .= srm

            # Apply reduction in DHW due to SRM
            dhw_t .= max.(0.0, dhw_t .- srm)
        end

        if is_guided && (in_seed_years || in_shade_years)
            # Update dMCDA values
            horizon = tstep:tstep+Int64(param_set("plan_horizon"))

            # Put more weight on projected conditions closer to the decision point
            env_horizon .= decay .* dhw_scen[horizon, :]
            mcda_vars.heat_stress_prob .= vec((mean(env_horizon, dims=1) .+ std(env_horizon, dims=1)) .* 0.5)

            env_horizon .= decay .* wave_scen[horizon, :]
            mcda_vars.dam_prob .= vec((mean(env_horizon, dims=1) .+ std(env_horizon, dims=1)) .* 0.5)
        end
        if is_guided && (in_seed_years || in_shade_years)
            mcda_vars.sum_cover .= site_coral_cover

            # Determine connectivity strength
            # Account for cases where no coral cover
            in_conn, out_conn, strong_pred = connectivity_strength(domain.TP_data .* site_k_area(domain), vec(site_coral_cover))
            (prefseedsites, prefshadesites, rankings) = guided_site_selection(mcda_vars, MCDA_approach,
                seed_decision_years[tstep], shade_decision_years[tstep],
                prefseedsites, prefshadesites, rankings, in_conn[mcda_vars.site_ids], out_conn[mcda_vars.site_ids], strong_pred[mcda_vars.site_ids])

            # Log site ranks
            # First col only holds site index ids so skip (with 2:end)
            site_ranks[tstep, rankings[:, 1], :] = rankings[:, 2:end]
        elseif seed_corals && (in_seed_years || in_shade_years)
            # Unguided deployment, seed/shade corals anywhere, so long as available space > 0
            prefseedsites, prefshadesites = unguided_site_selection(prefseedsites, prefshadesites,
                seed_decision_years[tstep], shade_decision_years[tstep],
                n_site_int, vec(leftover_space_m²), depth_priority)

            site_ranks[tstep, prefseedsites, 1] .= 1.0
            site_ranks[tstep, prefshadesites, 2] .= 1.0
        end

        has_shade_sites::Bool = !all(prefshadesites .== 0)
        has_seed_sites::Bool = !all(prefseedsites .== 0)

        # Fog selected locations
        if (fogging > 0.0) && in_shade_years && has_shade_sites
            site_locs = prefshadesites

            dhw_t[site_locs] .= dhw_t[site_locs] .* (1.0 .- fogging)
            Yfog[tstep, site_locs] .= fogging
        end

        # Calculate and apply bleaching mortality
        # bleaching_mortality!(Sbl, felt_dhw, depth_coeff, tstep, site_data.depth_med, bleaching_sensitivity, dhw_t, a_adapt, n_adapt)
        Sbl .= Y_pstep[:, :]
        bleaching_mortality!(Sbl, dhw_t, depth_coeff, c_dist_t, c_dist_t1, @view(bleaching_mort[tstep, :, :]))

        # Apply seeding
        if seed_corals && in_seed_years && has_seed_sites
            # Calculate proportion to seed based on current available space
            scaled_seed = distribute_seeded_corals(vec(total_site_area), prefseedsites, vec(leftover_space_m²), seeded_area)

            # Seed each site
            @views Y_pstep[seed_sc, prefseedsites] .+= scaled_seed
            Yseed[tstep, :, prefseedsites] .= scaled_seed
        end

        # Calculate survivors from bleaching and wave stress
        # @views @. prop_loss = Sbl * Sw_t[p_step, :, :]

        # Note: ODE is run relative to `k` area, but values are otherwise recorded
        #       in relative to absolute area.
        # Update initial condition
        # @. tmp = ((Y_pstep * prop_loss) * total_site_area) / absolute_k_area
        @. tmp = (Sbl * total_site_area) / absolute_k_area
        replace!(tmp, Inf => 0.0, NaN => 0.0)
        growth.u0 .= tmp

        # X is cover relative to `k` (max. carrying capacity)
        # So we subtract from 1.0 to get leftover/available space, relative to `k`
        p.sXr .= max.(1.0 .- sum(tmp, dims=1), 0.0) .* tmp .* p.r  # leftover space * current cover * growth_rate
        p.X_mb .= tmp .* p.mb    # current cover * background mortality

        sol::ODESolution = solve(growth, solver, save_everystep=false, save_start=false,
            alg_hints=[:nonstiff], adaptive=false, dt=0.5) 
        # Using the last step from ODE above, proportionally adjust site coral cover
        # if any are above the maximum possible (i.e., the site `k` value)
        @views Y_cover[tstep, :, :] .= clamp.(sol.u[end] .* absolute_k_area ./ total_site_area, 0.0, 1.0)

        if tstep < tf
            adjust_population_distribution!(Y_cover, n_groups, c_dist_t, c_dist_t1, tstep, c_increase)
        end
    end

    # Avoid placing importance on sites that were not considered
    # (lower values are higher importance)
    site_ranks[site_ranks.==0.0] .= n_sites + 1
    return (raw=Y_cover, seed_log=Yseed, fog_log=Yfog, shade_log=Yshade, site_ranks=site_ranks, bleaching_mortality=bleaching_mort)
end
