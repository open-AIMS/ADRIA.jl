"""Scenario running functions"""

using ADRIA.metrics: relative_cover, total_absolute_cover, absolute_shelter_volume, relative_shelter_volume
using ADRIA.metrics: relative_juveniles, relative_taxa_cover, juvenile_indicator


"""
    setup_cache(domain::Domain)::NamedTuple

Establish tuple of matrices/vectors for use as reusable data stores to avoid repeated memory allocations.
"""
function setup_cache(domain::Domain)::NamedTuple

    # sim constants
    n_locations::Int64 = domain.coral_growth.n_locations
    n_species::Int64 = domain.coral_growth.n_species
    n_groups::Int64 = domain.coral_growth.n_groups

    init_cov = domain.init_coral_cover
    cache = (
        sf=zeros(n_groups, n_locations),  # stressed fecundity
        fec_all=zeros(size(init_cov)...),  # all fecundity
        fec_scope=zeros(n_groups, n_locations),  # fecundity scope
        prop_loss=zeros(n_species, n_locations),  # proportional loss
        Sbl=zeros(n_species, n_locations),   # bleaching survivors
        dhw_step=zeros(n_locations),  # DHW each time step
        cov_tmp=zeros(size(init_cov)...),  # Cover for previous timestep
        felt_dhw=zeros(size(init_cov)...),  # Store for felt DHW (DHW after reductions)
        depth_coeff=zeros(n_locations),  # store for depth coefficient
        location_area=Matrix{Float64}(domain.location_data.area'),  # location areas
        TP_data=Matrix{Float64}(domain.TP_data),  # transition probabilities
        waves=zeros(length(timesteps(domain)), n_species, n_locations)
    )

    return cache
end


"""
    run_scenarios(param_df::DataFrame, domain::Domain; remove_workers=true)
    run_scenarios(param_df::DataFrame, domain::Domain, rcp::String; remove_workers=true)
    run_scenarios(param_df::DataFrame, domain::Domain, rcp::Array{String}; remove_workers=true)

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
julia> rs_all = ADRIA.run_scenarios(p_df, dom)
```

# Arguments
- param_df : DataFrame of scenarios to run
- domain : Domain, to run scenarios with
- rcp : ID of RCP(s) to run scenarios under.
- remove_workers : if running in parallel, removes workers after completion

# Returns
ResultSet
"""
function run_scenarios(param_df::DataFrame, domain::Domain; remove_workers=true)::ResultSet
    # Identify available data
    avail_data::Vector{String} = readdir(joinpath(domain.env_layer_md.dpkg_path, "DHWs"))
    RCP_ids = replace.(avail_data, "dhwRCP" => "", ".nc" => "")

    @info "Running scenarios for RCPs: $(RCP_ids)"
    return run_scenarios(param_df, domain, RCP_ids::Array{String}; remove_workers=remove_workers)
end
function run_scenarios(param_df::DataFrame, domain::Domain, RCP::String; show_progress=true, remove_workers=true)::ResultSet
    setup()
    parallel = (nrow(param_df) > 4096) && (parse(Bool, ENV["ADRIA_DEBUG"]) == false)
    if parallel
        _setup_workers()
        sleep(2)  # wait a bit while workers spin-up
        @eval @everywhere using ADRIA
    end

    domain = switch_RCPs!(domain, RCP)
    domain, data_store = ADRIA.setup_result_store!(domain, param_df)

    cache = setup_cache(domain)
    run_msg = "Running $(nrow(param_df)) scenarios for RCP $RCP"

    # Convert to named matrix for faster iteration
    scen_mat = NamedDimsArray(Matrix(param_df); scenarios=1:nrow(param_df), factors=names(param_df))

    # Batch run scenarios
    func = (dfx) -> run_scenario(dfx..., domain, data_store, cache)
    if parallel
        if show_progress
            @showprogress run_msg 4 pmap(func, enumerate(eachrow(scen_mat)))
        else
            pmap(func, enumerate(eachrow(scen_mat)))
        end

        if remove_workers
            _remove_workers()
        end
    else
        if show_progress
            @showprogress run_msg 4 map(func, enumerate(eachrow(scen_mat)))
        else
            map(func, enumerate(eachrow(scen_mat)))
        end
    end

    return load_results(domain)
end
function run_scenarios(param_df::DataFrame, domain::Domain, RCP_ids::Array{String}; show_progress=true, remove_workers=true)::ResultSet
    @info "Running $(nrow(param_df)) scenarios across $(length(RCP_ids)) RCPs"

    setup()
    output_dir = ENV["ADRIA_OUTPUT_DIR"]

    result_sets::Vector{ResultSet} = Vector{ResultSet}(undef, length(RCP_ids))
    for (i, RCP) in enumerate(RCP_ids)
        tmp_dir = mktempdir(prefix="ADRIA_")
        ENV["ADRIA_OUTPUT_DIR"] = tmp_dir

        result_sets[i] = run_scenarios(param_df, domain, RCP; show_progress=show_progress, remove_workers=false)
    end

    if remove_workers
        _remove_workers()
    end

    ENV["ADRIA_OUTPUT_DIR"] = output_dir

    rs = combine_results(result_sets...)

    # Remove temporary result dirs
    for t in result_sets
        rm(result_location(t); force=true, recursive=true)
    end

    return rs
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
Logs of location ranks only store the mean location rankings over all environmental scenarios.
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
    vals = total_absolute_cover(rs_raw, location_area(domain))
    vals[vals.<threshold] .= 0.0
    data_store.total_absolute_cover[:, :, idx] .= vals

    vals .= absolute_shelter_volume(rs_raw, location_area(domain), param_set)
    vals[vals.<threshold] .= 0.0
    data_store.absolute_shelter_volume[:, :, idx] .= vals

    vals .= relative_shelter_volume(rs_raw, location_area(domain), location_k_area(domain), param_set)
    vals[vals.<threshold] .= 0.0
    data_store.relative_shelter_volume[:, :, idx] .= vals

    coral_spec::DataFrame = to_coral_spec(param_set)
    vals .= relative_juveniles(rs_raw, coral_spec)
    vals[vals.<threshold] .= 0.0
    data_store.relative_juveniles[:, :, idx] .= vals

    vals .= juvenile_indicator(rs_raw, coral_spec, location_area(domain), location_k_area(domain))
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
    tmp_location_ranks = zeros(Float32, tf, nrow(domain.location_data), length(domain.interventions))
    c_dim = Base.ndims(result_set.raw) + 1
    log_stores = (:location_ranks, :seed_log, :fog_log, :shade_log)
    for k in log_stores
        if k == :seed_log || k == :location_ranks
            concat_dim = c_dim
        else
            concat_dim = c_dim - 1
        end

        vals = getfield(result_set, k)
        vals[vals.<threshold] .= 0.0

        if k == :seed_log
            getfield(data_store, k)[:, :, :, idx] .= vals
        elseif k == :location_ranks
            tmp_location_ranks[:, :, :] .= vals
        else
            getfield(data_store, k)[:, :, idx] .= vals
        end
    end

    if !isnothing(data_store.location_ranks)
        # Squash location ranks down to average rankings over environmental repeats
        data_store.location_ranks[:, :, :, idx] .= tmp_location_ranks
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
Only the mean location rankings are kept

# Returns
NamedTuple of collated results
"""
function run_model(domain::Domain, param_set::DataFrameRow, corals::DataFrame, cache::NamedTuple)::NamedTuple
    ps = NamedDimsArray(Vector(param_set), factors=names(param_set))
    return run_model(domain, ps, corals, cache)
end
function run_model(domain::Domain, param_set::NamedDimsArray, corals::DataFrame, cache::NamedTuple)::NamedTuple
    sim_params = domain.sim_constants
    location_data = domain.location_data
    p = domain.coral_growth.ode_p

    # Set random seed using intervention values
    # TODO: More robust way of getting intervention/criteria values
    rnd_seed_val::Int64 = floor(Int64, sum(param_set(!=("RCP"))))  # select everything except RCP
    Random.seed!(rnd_seed_val)

    ### TODO: All cached arrays/values to be moved to outer function and passed in
    # to reduce overall allocations (e.g., sim constants don't change across all scenarios)
    dhw_idx::Int64 = Int64(param_set("dhw_scenario"))
    wave_idx::Int64 = Int64(param_set("wave_scenario"))

    dhw_scen::Matrix{Float64} = @view domain.dhw_scens[:, :, dhw_idx]

    # TODO: Better conversion of Ub to wave mortality
    #       Currently scaling significant wave height by its max to non-dimensionalize values
    wave_scen::Matrix{Float64} = Matrix{Float64}(domain.wave_scens[:, :, wave_idx]) ./ maximum(domain.wave_scens[:, :, wave_idx])

    tspan::Tuple = (0.0, 1.0)
    solver::BS3 = BS3()

    MCDA_approach::Int64 = param_set("guided")

    # sim constants
    tf::Int64 = size(dhw_scen, 1)
    n_iv_locs::Int64 = sim_params.n_iv_locs
    n_locs::Int64 = domain.coral_growth.n_locations
    n_species::Int64 = domain.coral_growth.n_species
    n_groups::Int64 = domain.coral_growth.n_groups

    # years to start seeding/shading
    seed_start_year::Int64 = param_set("seed_year_start")
    shade_start_year::Int64 = param_set("shade_year_start")

    n_TA_to_seed::Int64 = param_set("seed_TA")  # tabular Acropora size class 2, per year per species per cluster
    n_CA_to_seed::Int64 = param_set("seed_CA")  # corymbose Acropora size class 2, per year per species per cluster
    fogging::Real = param_set("fogging")  # percent reduction in bleaching mortality through fogging
    srm::Real = param_set("SRM")  # DHW equivalents reduced by some shading mechanism
    seed_years::Int64 = param_set("seed_years")  # number of years to seed
    shade_years::Int64 = param_set("shade_years")  # number of years to shade

    ### END TODO

    total_location_area::Array{Float64,2} = cache.location_area
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

    Y_cover::Array{Float64,3} = zeros(tf, n_species, n_locs)  # Coral cover relative to total location area
    Y_cover[1, :, :] .= domain.init_coral_cover
    ode_u = zeros(n_species, n_locs)
    cover_tmp = p.cover  # pre-allocated matrix used to avoid memory allocations

    location_ranks = SparseArray(zeros(tf, n_locs, length(domain.interventions))) # log seeding/fogging/shading ranks
    Yshade = SparseArray(spzeros(tf, n_locs))
    Yfog = SparseArray(spzeros(tf, n_locs))
    Yseed = SparseArray(zeros(tf, 2, n_locs))  # 2 = the two enhanced coral types

    # Intervention strategy: < 0 is no intervention, 0 is random location selection, > 0 is guided
    is_guided = param_set("guided") > 0

    # Years at which to reassess seeding location selection
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

    int_log = NamedDimsArray(hcat(seed_decision_years, shade_decision_years, shade_decision_years), year=1:tf, log=[:seed, :fog, :shade])
    pref_locations = (seed=zeros(Int64, n_iv_locs), shade=zeros(Int64, n_iv_locs), fog=zeros(Int64, n_iv_locs))

    # Max coral cover at each location. Divided by 100 to convert to proportion
    max_cover = location_k(domain)

    # Set other params for ODE
    p.r .= corals.growth_rate  # Assumed growth_rate

    p.mb .= corals.mb_rate  # background mortality
    @set! p.k = max_cover  # max coral cover
    @set! p.comp = sim_params.comp  # competition rate

    # Proportionally adjust initial cover (handles inappropriate initial conditions)
    proportional_adjustment!(Y_cover[1, :, :], cover_tmp, max_cover)

    # Define constant table location for seed values
    tabular_enhanced::BitArray = corals.taxa_id .== 1
    corymbose_enhanced::BitArray = corals.taxa_id .== 3
    target_class_id::BitArray = corals.class_id .== 2  # seed second smallest size class
    seed_sc_TA::Int64 = first(findall(tabular_enhanced .& target_class_id))  # size class indices for TA and CA
    seed_sc_CA::Int64 = first(findall(corymbose_enhanced .& target_class_id))

    # Extract colony areas for locations selected and convert to m^2
    col_area_seed_TA = corals.colony_area_cm2[seed_sc_TA] / 10^4
    col_area_seed_CA = corals.colony_area_cm2[seed_sc_CA] / 10^4

    bleaching_sensitivity = corals.bleaching_sensitivity

    # Defaults to considering all locations if depth cannot be considered.
    depth_priority = collect(1:nrow(location_data))

    # calculate total area to seed respecting tolerance for minimum available space to still seed at a location
    area_to_seed = (col_area_seed_TA .* n_TA_to_seed) + (col_area_seed_CA .* n_CA_to_seed)

    # Filter out locations outside of desired depth range
    if .!all(location_data.depth_med .== 0)
        max_depth::Float64 = param_set("depth_min") + param_set("depth_offset")
        depth_criteria::BitArray{1} = (location_data.depth_med .>= param_set("depth_min")) .& (location_data.depth_med .<= max_depth)

        # TODO: Include this change in MATLAB version as well
        if any(depth_criteria .> 0)
            # If locations can be filtered based on depth, do so.
            depth_priority = depth_priority[depth_criteria]
        else
            # Otherwise if no locations can be filtered, remove depth as a criterion.
            @warn "No locations within provided depth range of $(param_set("depth_min")) - $(max_depth) meters. Considering all locations."
        end
    end

    if is_guided
        # pre-allocate rankings
        rankings, criteria_store, thresholds_store, min_distance, use_dist = initialize_mcda(domain, param_set, depth_priority, area_to_seed)
    end

    #### End coral constants

    ## Update ecological parameters based on intervention option
    # Set up assisted adaptation values
    a_adapt = zeros(n_species)
    a_adapt[tabular_enhanced] .= param_set("a_adapt")
    a_adapt[corymbose_enhanced] .= param_set("a_adapt")

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

    Sw_t .= 1.0 .- Sw_t

    # Flag indicating whether to seed or not to seed
    seed_corals::Bool = (n_TA_to_seed > 0) || (n_CA_to_seed > 0)

    location_k_prop = max_cover'
    absolute_k_area = total_location_area .* location_k_prop  # max possible coral area in m^2
    growth::ODEProblem = ODEProblem{true}(growthODE, ode_u, tspan, p)
    tmp::Matrix{Float64} = zeros(size(Y_cover[1, :, :]))  # temporary array to hold intermediate covers

    env_horizon = zeros(Int64(param_set("plan_horizon") + 1), n_sites)  # temporary cache for planning horizon

    # basal_area_per_settler is the area in m^2 of a size class one coral 
    basal_area_per_settler = corals.colony_area_cm2[corals.class_id.==1] ./ 100 .^ 2
    # debug_log = zeros(tf, n_sites)
    # rec_log = zeros(tf, 6, n_sites)
    # juv_log = zeros(tf, 12, n_sites)
    # juves = [1, 2, 7, 8, 13, 14, 19, 20, 25, 26, 31, 32]
    @inbounds for tstep::Int64 in 2:tf
        p_step::Int64 = tstep - 1
        Y_pstep[:, :] .= Y_cover[p_step, :, :]

        sf .= stressed_fecundity(tstep, a_adapt, n_adapt, dhw_scen[p_step, :],
            LPdhwcoeff, DHWmaxtot, LPDprm2, n_groups)

        # Calculates scope for coral fedundity for each size class and at each location.
        fecundity_scope!(fec_scope, fec_all, fec_params_per_m², Y_pstep, total_location_area)

        location_coral_cover = sum(Y_pstep, dims=1)  # dims: nlocations * 1
        leftover_space_prop = relative_leftover_space(location_k_prop, location_coral_cover)
        leftover_space_m² = leftover_space_prop .* total_location_area

        # Recruitment represents additional cover, relative to total location area
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
            update_criteria_store!(criteria_store, sum(Sw_t[tstep, :, :], dims=1)', dhw_t, domain.in_conn, domain.out_conn,
                total_location_area, location_coral_cover, domain.location_data, depth_priority)
        end
        if is_guided && (in_seed_years || in_shade_years)
            mcda_vars.sum_cover .= site_coral_cover

        if is_guided && in_seed_years
            (pref_locations, rankings) = guided_location_selection(criteria_store, domain.interventions,
                param_set, thresholds_store, n_iv_locs, domain.location_distances, min_distance, use_dist,
                int_log(year=tstep), pref_locations, rankings)

            # Log location ranks
            # First col only holds location index ids so skip (with 2:end)
            location_ranks[tstep, rankings[:seed][:, 1], 1] .= rankings[:seed][:, 2]
            location_ranks[tstep, rankings[:fog][:, 1], 2] .= rankings[:fog][:, 2]
            if int_log(year=tstep, log=:shade)
                for cluster in eachrow(rankings[:shade])
                    location_ranks[tstep, dropdims(cluster_ids(domain) .== cluster[1], dims=1), 3] .= cluster[2]
                end
            end
        elseif seed_corals && (in_seed_years || in_shade_years)
            # Unguided deployment, seed/shade corals anywhere, so long as available space > 0
            pref_locations = unguided_location_selection(pref_locations, int_log(year=tstep),
                n_iv_locs, vec(leftover_space_m²), depth_priority, unique(cluster_ids(domain)))

            location_ranks[tstep, pref_locations.seed[pref_locations.seed.>0], 1] .= 1.0
            location_ranks[tstep, pref_locations.fog[pref_locations.fog.>0], 2] .= 1.0
            locations_in_clusters = Bool.(dropdims(sum(in.(pref_locations.shade, cluster_ids(domain)), dims=1), dims=1))
            location_ranks[tstep, locations_in_clusters, 3] .= 1.0
        end

        has_fog_locations::Bool = !all(pref_locations.fog .== 0)
        has_seed_locations::Bool = !all(pref_locations.seed .== 0)

        if (srm > 0.0) && in_shade_years
            locations_in_clusters = Bool.(dropdims(sum(in.(pref_locations.shade, cluster_ids(domain)), dims=1), dims=1))

            Yshade[tstep, locations_in_clusters] .= srm

            # Apply reduction in DHW due to SRM
            dhw_t[locations_in_clusters] .= max.(0.0, dhw_t[locations_in_clusters] .- srm)
        end

        if (fogging > 0.0) && in_shade_years && (has_seed_locations || has_fog_locations)
            if has_seed_locations
                # Always fog seeded locations if they are selected
                location_locs::Vector{Int64} = pref_locations.seed[pref_locations.seed.>0]
            elseif has_fog_locations
                # Use locations selected for fogging otherwise
                location_locs = pref_locations.fog[pref_locations.fog.>0]
            end

            dhw_t[location_locs] .= dhw_t[location_locs] .* (1.0 .- fogging)
            Yfog[tstep, location_locs] .= fogging
        end

        # Calculate and apply bleaching mortality
        bleaching_mortality!(Sbl, felt_dhw, depth_coeff, tstep, location_data.depth_med, bleaching_sensitivity, dhw_t, a_adapt, n_adapt)

        # Apply seeding
        if seed_corals && in_seed_years && has_seed_locations
            # Calculate proportion to seed based on current available space
            seeded_area = (TA=n_TA_to_seed * col_area_seed_TA, CA=n_CA_to_seed * col_area_seed_CA)
            pref_location_seed = pref_locations.seed[pref_locations.seed.>0]
            scaled_seed = distribute_seeded_corals(vec(total_location_area), pref_location_seed, vec(leftover_space_m²), seeded_area)

            # Seed each location with TA or CA
            @views Y_pstep[seed_sc_TA, pref_location_seed] .= Y_pstep[seed_sc_TA, pref_location_seed] .+ scaled_seed.TA
            @views Y_pstep[seed_sc_CA, pref_location_seed] .= Y_pstep[seed_sc_CA, pref_location_seed] .+ scaled_seed.CA

            # Log seed values/locations (these values are relative to location area)
            Yseed[tstep, 1, pref_location_seed] .= scaled_seed.TA
            Yseed[tstep, 2, pref_location_seed] .= scaled_seed.CA
        end

        # Calculate survivors from bleaching and wave stress
        @views @. prop_loss = Sbl * Sw_t[p_step, :, :]

        # Note: ODE is run relative to `k` area, but values are otherwise recorded
        #       in relative to absolute area.
        # Update initial condition
        @. tmp = ((Y_pstep * prop_loss) * total_site_area) / absolute_k_area
        replace!(tmp, Inf => 0.0, NaN => 0.0)
        growth.u0 .= tmp

        # X is cover relative to `k` (max. carrying capacity)
        # So we subtract from 1.0 to get leftover/available space, relative to `k`
        p.sXr .= max.(1.0 .- sum(tmp, dims=1), 0.0) .* tmp .* p.r  # leftover space * current cover * growth_rate
        p.X_mb .= tmp .* p.mb    # current cover * background mortality

        # Background mortality of small massives (incorporates competition with larger tabular acroporas)
        p.M_sm .= tmp[p.small_massives, :] .* (p.mb[p.small_massives] .+ p.comp .* sum(tmp[p.acr_6_12, :], dims=1))

        # Space gained from small massives due to competition
        p.sm_comp .= p.comp .* sum(tmp[p.small_massives, :], dims=1)

        sol::ODESolution = solve(growth, solver, save_everystep=false, save_start=false,
            alg_hints=[:nonstiff], abstol=1e-6, reltol=1e-6)  # , adaptive=false, dt=1.0
        # Using the last step from ODE above, proportionally adjust site coral cover
        # if any are above the maximum possible (i.e., the site `k` value)
        @views Y_cover[tstep, :, :] .= clamp.(sol.u[end] .* absolute_k_area ./ total_location_area, 0.0, 1.0)
    end

    # Avoid placing importance on locations that were not considered
    # (lower values are higher importance)
    location_ranks[location_ranks.==0.0] .= n_locs + 1
    return (raw=Y_cover, seed_log=Yseed, fog_log=Yfog, shade_log=Yshade, location_ranks=location_ranks)
end
