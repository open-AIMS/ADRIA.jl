"""Scenario running functions"""

import ADRIA.metrics: relative_cover, total_absolute_cover, absolute_shelter_volume, relative_shelter_volume

"""
    setup_cache(domain::Domain)::NamedTuple

Establish tuple of matrices/vectors for use as reusable data stores to avoid repeated memory allocations.
"""
function setup_cache(domain::Domain)::NamedTuple

    # sim constants
    n_sites::Int64 = domain.coral_growth.n_sites
    n_species::Int64 = domain.coral_growth.n_species
    n_groups::Int64 = domain.coral_growth.n_groups

    # Strip names from NamedArrays
    init_cov = Matrix{Float64}(domain.init_coral_cover)

    cache = (
        sf=zeros(n_groups, n_sites),
        fec_all=zeros(size(init_cov)...),
        fec_scope=zeros(n_groups, n_sites),
        prop_loss=zeros(n_species, n_sites),
        Sbl=zeros(n_species, n_sites),
        dhw_step=zeros(n_sites),
        init_cov=init_cov,
        cov_tmp=zeros(size(init_cov)...),
        site_area=Array{Float64}(domain.site_data.area'),
        TP_data=Array{Float64,2}(domain.TP_data)
    )

    return cache
end




"""
    run_scenarios(param_df::DataFrame, domain::Domain)
    run_scenarios(param_df::DataFrame, domain::Domain, rcp::String)
    run_scenarios(param_df::DataFrame, domain::Domain, rcp::Array{String})

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

# Returns
ResultSet
"""
function run_scenarios(param_df::DataFrame, domain::Domain)::ResultSet

    # Identify available data
    avail_data::Vector{String} = readdir(joinpath(domain.env_layer_md.dpkg_path, "DHWs"))
    RCP_ids = replace.(avail_data, "dhwRCP" => "", ".mat" => "")

    @info "Running scenarios for RCPs: $(RCP_ids)"
    return run_scenarios(param_df, domain, RCP_ids::Array{String})
end
function run_scenarios(param_df::DataFrame, domain::Domain, RCP::String)::ResultSet
    setup()
    parallel = (nrow(param_df) > 256) && (parse(Bool, ENV["ADRIA_DEBUG"]) == false)

    switch_RCPs!(domain, RCP)
    domain, data_store = ADRIA.setup_result_store!(domain, param_df)
    cache = setup_cache(domain)
    run_msg = "Running $(nrow(param_df)) scenarios for RCP $RCP"

    # Batch run scenarios
    func = (dfx) -> run_scenario(dfx, domain, data_store, cache)
    if parallel
        @eval @everywhere using ADRIA

        @showprogress run_msg 4 pmap(func, enumerate(eachrow(param_df)))
    else
        @showprogress run_msg 1 map(func, enumerate(eachrow(param_df)))
    end

    return load_results(domain)
end
function run_scenarios(param_df::DataFrame, domain::Domain, RCP_ids::Array{String})::ResultSet
    setup()
    parallel = (nrow(param_df) > 256) && (parse(Bool, ENV["ADRIA_DEBUG"]) == false)

    if length(unique(RCP_ids)) != length(RCP_ids)
        # Disallow running duplicate RCP scenarios
        error("Duplicate RCP ids specified")
    end

    run_msg = "Running $(nrow(param_df)) scenarios "
    output_dir = ENV["ADRIA_OUTPUT_DIR"]
    tmp_result_dirs::Vector{String} = String[]

    for RCP in RCP_ids
        switch_RCPs!(domain, RCP)
        tmp_dir = mktempdir(prefix="ADRIA_")
        ENV["ADRIA_OUTPUT_DIR"] = tmp_dir

        domain, data_store = ADRIA.setup_result_store!(domain, param_df)
        cache = setup_cache(domain)

        push!(tmp_result_dirs, result_location(domain))

        # Batch run scenarios
        func = (dfx) -> run_scenario(dfx, domain, data_store, cache)
        msg = run_msg * "for RCP $RCP"
        if parallel
            @eval @everywhere using ADRIA

            @showprogress msg 4 pmap(func, enumerate(eachrow(param_df)))
        else
            @showprogress msg 1 map(func, enumerate(eachrow(param_df)))
        end
    end

    ENV["ADRIA_OUTPUT_DIR"] = output_dir
    rs = combine_results(tmp_result_dirs)

    # Remove temporary result dirs
    for t in tmp_result_dirs
        rm(t; force=true, recursive=true)
    end

    return rs
end


"""
    run_scenarios(scen::Tuple{Int, DataFrameRow}, domain::Domain, data_store::NamedTuple, cache::NamedTuple)
    run_scenario(scen::Tuple{Int, DataFrameRow}, domain::Domain, data_store::NamedTuple)

Run individual scenarios for a given domain.

Stores results on disk in Zarr format at pre-configured location.
Sets up a new `cache` if not provided.

# Notes
Logs of site ranks only store the mean site rankings over all environmental scenarios.
This is to reduce the volume of data stored.
"""
function run_scenario(scen::Tuple{Int,DataFrameRow}, domain::Domain, data_store::NamedTuple, cache::NamedTuple)
    # Update model with values in given DF row
    update_params!(domain, scen[2])

    dhw_scen = scen[2].dhw_scenario
    wave_scen = scen[2].wave_scenario

    # TODO: Modify all scenario constants here to avoid repeated allocations
    @set! domain.coral_growth.ode_p.k = (domain.site_data.k::Vector{Float64} / 100.0)  # Max possible cover at site
    @set! domain.coral_growth.ode_p.comp = domain.sim_constants.comp::Float64  # competition rate between two mature coral groups

    run_scenario(domain; idx=scen[1], dhw=dhw_scen, wave=wave_scen, data_store=data_store, cache=cache)
end
function run_scenario(scen::Tuple{Int,DataFrameRow}, domain::Domain, data_store::NamedTuple)
    run_scenario(scen, domain, data_store, setup_cache(domain))
end


"""
    run_scenario(domain::Domain; idx=1, dhw=1, wave=1, data_store::NamedTuple, cache::NamedTuple)::NamedTuple

Convenience function to directly run a scenario for a Domain with pre-set values.

Stores results on disk in Zarr format at pre-configured location.

# Notes
Logs of site ranks only store the mean site rankings over all environmental scenarios.
This is to reduce the volume of data stored.
"""
function run_scenario(domain::Domain; idx::Int=1, dhw::Int=1, wave::Int=1, data_store::NamedTuple, cache::NamedTuple)

    # Extract non-coral parameters
    df = DataFrame(domain.model)
    not_coral_params = df[!, :component] .!== Coral
    param_set = NamedTuple{tuple(df[not_coral_params, :fieldname]...)}(df[not_coral_params, :val])

    # Expand coral model to include its specifications across all taxa/species/groups
    coral_params = to_spec(component_params(domain.model, Coral))

    # Pass in environmental layer data stripped of named dimensions.
    all_dhws = Array{Float64}(domain.dhw_scens)
    all_waves = Array{Float64}(domain.wave_scens)

    result_set = run_scenario(domain, param_set, coral_params, domain.sim_constants, domain.site_data,
        domain.coral_growth.ode_p,
        all_dhws[:, :, dhw], all_waves[:, :, wave], cache)

    # Capture results to disk
    # Set values below threshold to 0 to save space
    tf = size(all_dhws, 1)
    threshold = parse(Float32, ENV["ADRIA_THRESHOLD"])

    tmp_site_ranks = zeros(Float32, tf, nrow(domain.site_data), 2)

    r_raw = result_set.raw
    vals = relative_cover(r_raw)
    vals[vals.<threshold] .= 0.0
    data_store.relative_cover[:, :, idx] .= vals

    p_tbl = param_table(domain)
    vals .= absolute_shelter_volume(r_raw, site_area(domain), p_tbl)
    vals[vals.<threshold] .= 0.0
    data_store.absolute_shelter_volume[:, :, idx] .= vals

    vals .= relative_shelter_volume(r_raw, site_area(domain), p_tbl)
    vals[vals.<threshold] .= 0.0
    data_store.relative_shelter_volume[:, :, idx] .= vals

    # Store raw results if no metrics specified
    # if length(metrics) == 0
    #     data_store.raw[:, :, :, idx] .= r.raw
    # end

    # Store logs
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


"""
    run_scenario(param_row::DataFrameRow, domain::Domain)::NamedTuple
    run_scenario(param_set::NamedTuple, domain::Domain)::NamedTuple

Run a single scenario and return results.
"""
function run_scenario(param_set::NamedTuple, domain::Domain)::NamedTuple

    # Expand coral model to include its specifications across all taxa/species/groups
    coral_params = to_spec(component_params(domain.model, Coral))

    dhw_rep_id = param_set.dhw_scenario
    wave_rep_id = param_set.wave_scenario

    cache = setup_cache(domain)
    return run_scenario(domain, param_set, coral_params, domain.sim_constants, domain.site_data,
        domain.coral_growth.ode_p,
        Matrix{Float64}(domain.dhw_scens[:, :, dhw_rep_id]),
        Matrix{Float64}(domain.wave_scens[:, :, wave_rep_id]), cache)
end
function run_scenario(param_row::DataFrameRow, domain::Domain)::NamedTuple
    # Update model with values in given DF row
    update_params!(domain, param_row)
    param_set = NamedTuple{domain.model[:fieldname]}(domain.model[:val])

    return run_scenario(param_set, domain)
end

"""
    run_scenario(domain, param_set, corals, sim_params, site_data, p::NamedTuple,
                 dhw_scen::Array, wave_scen::Array, cache::NamedTuple)::NamedTuple

Core scenario running function.

# Notes
Only the mean site rankings are kept
"""
function run_scenario(domain::Domain, param_set::NamedTuple, corals::DataFrame, sim_params::SimConstants, site_data::DataFrame,
    p::NamedTuple, dhw_scen::Matrix{Float64},
    wave_scen::Matrix{Float64}, cache::NamedTuple)::NamedTuple

    # Set random seed using intervention values
    # TODO: More robust way of getting intervention/criteria values
    rnd_seed_val = floor(Int, sum(values(param_set)))
    Random.seed!(rnd_seed_val)

    ### TODO: All cached arrays/values to be moved to outer function and passed in
    # to reduce overall allocations (e.g., sim constants don't change across all scenarios)

    tspan::Tuple = (0.0, 1.0)
    solver::RK4 = RK4()

    MCDA_approach::Int64 = param_set.guided

    # sim constants
    n_sites::Int64 = domain.coral_growth.n_sites
    nsiteint::Int64 = sim_params.nsiteint
    tf::Int64 = size(dhw_scen, 1)
    n_species::Int64 = domain.coral_growth.n_species
    n_groups::Int64 = domain.coral_growth.n_groups

    # years to start seeding/shading
    seed_start_year::Int64 = param_set.seed_year_start
    shade_start_year::Int64 = param_set.shade_year_start

    n_TA_to_seed::Int64 = param_set.seed_TA  # tabular Acropora size class 2, per year per species per cluster
    n_CA_to_seed::Int64 = param_set.seed_CA  # corymbose Acropora size class 2, per year per species per cluster
    fogging::Real = param_set.fogging  # percent reduction in bleaching mortality through fogging
    srm::Real = param_set.SRM  # DHW equivalents reduced by some shading mechanism
    seed_years::Int64 = param_set.seed_years  # number of years to seed
    shade_years::Int64 = param_set.shade_years  # number of years to shade

    ### END TODO

    total_site_area::Array{Float64,2} = cache.site_area
    fec_params_per_m²::Vector{Float64} = corals.fecundity  # number of larvae produced per m²

    # Caches
    TP_data = cache.TP_data
    # fec_all::Array{Float64,2} = zeros(size(init_cov)...)
    # fec_scope::Array{Float64,2} = zeros(n_groups, n_sites)
    # prop_loss::Array{Float64,2} = zeros(n_species, n_sites)
    # Sbl::Array{Float64,2} = zeros(n_species, n_sites)
    # dhw_step::Vector{Float64} = zeros(n_sites)
    sf = cache.sf[:, :]
    fec_all = cache.fec_all[:, :]
    fec_scope = cache.fec_scope[:, :]
    prop_loss = cache.prop_loss[:, :]
    Sbl = cache.Sbl[:, :]
    dhw_t = cache.dhw_step[:]
    Y_pstep = cache.cov_tmp[:, :]

    Y_cover::Array{Float64,3} = zeros(tf, n_species, n_sites)  # Coral cover relative to total site area
    Y_cover[1, :, :] .= cache.init_cov[:, :]
    ode_u = zeros(n_species, n_sites)
    cover_tmp = p.cover  # pre-allocated matrix used to avoid memory allocations

    site_ranks = SparseArray(zeros(tf, n_sites, 2)) # log seeding/fogging/shading ranks
    Yshade = SparseArray(spzeros(tf, n_sites))
    Yfog = SparseArray(spzeros(tf, n_sites))
    Yseed = SparseArray(zeros(tf, 2, n_sites))  # 2 = the two enhanced coral types

    # Intervention strategy: 0 is random, > 0 is guided
    is_guided = param_set.guided > 0

    # Years at which to reassess seeding site selection
    seed_decision_years = repeat([false], tf)
    shade_decision_years = repeat([false], tf)

    seed_start_year = max(seed_start_year, 2)
    if param_set.seed_freq > 0
        max_consider = min(seed_start_year + seed_years - 1, tf)
        seed_decision_years[seed_start_year:param_set.seed_freq:max_consider] .= true
    else
        # Start at year 2 or the given specified seed start year
        seed_decision_years[seed_start_year] = true
    end

    shade_start_year = max(shade_start_year, 2)
    if param_set.shade_freq > 0
        max_consider = min(shade_start_year + shade_years - 1, tf)
        shade_decision_years[shade_start_year:param_set.shade_freq:max_consider] .= true
    else
        # Start at year 2 or the given specified shade start year
        shade_decision_years[shade_start_year] = true
    end

    prefseedsites::Vector{Int64} = zeros(Int, nsiteint)
    prefshadesites::Vector{Int64} = zeros(Int, nsiteint)

    # Max coral cover at each site. Divided by 100 to convert to proportion
    max_cover = site_data.k / 100.0

    # Set other params for ODE
    p.r .= corals.growth_rate  # Assumed growth_rate
    p.mb .= corals.mb_rate  # background mortality

    # Proportionally adjust initial cover (handles inappropriate initial conditions)
    Y_cover[1, :, :] .= proportional_adjustment!(Y_cover[1, :, :], cover_tmp, max_cover)

    # Define constant table location for seed values
    tabular_enhanced::BitArray = corals.taxa_id .== 1
    corymbose_enhanced::BitArray = corals.taxa_id .== 3
    target_class_id::BitArray = corals.class_id .== 2  # seed second smallest size class
    seed_sc_TA::Int64 = first(findall(tabular_enhanced .& target_class_id))  # size class indices for TA and CA
    seed_sc_CA::Int64 = first(findall(corymbose_enhanced .& target_class_id))

    # extract colony areas for sites selected and convert to m^2
    col_area_seed_TA = corals.colony_area_cm2[seed_sc_TA] / 10^4
    col_area_seed_CA = corals.colony_area_cm2[seed_sc_CA] / 10^4

    if is_guided
        ## Weights for connectivity , waves (ww), high cover (whc) and low
        wtwaves = param_set.wave_stress # weight of wave damage in MCDA
        wtheat = param_set.heat_stress # weight of heat damage in MCDA
        wtconshade = param_set.shade_connectivity # weight of connectivity for shading in MCDA
        wtinconnseed = param_set.in_seed_connectivity # weight for seed sites with high number of incoming connections
        wtoutconnseed = param_set.out_seed_connectivity # weight for seed sites with high number of outgoing connections

        wthicover = param_set.coral_cover_high # weight of high coral cover in MCDA (high cover gives preference for seeding corals but high for SRM)
        wtlocover = param_set.coral_cover_low # weight of low coral cover in MCDA (low cover gives preference for seeding corals but high for SRM)
        wtpredecseed = param_set.seed_priority # weight for the importance of seeding sites that are predecessors of priority reefs
        wtpredecshade = param_set.shade_priority # weight for the importance of shading sites that are predecessors of priority reefs
        risktol = param_set.deployed_coral_risk_tol # risk tolerance
        covertol = param_set.coral_cover_tol

        # Defaults to considering all sites if depth cannot be considered.
        depth_priority = collect(1:nrow(site_data))

        # calculate total area to seed
        area_to_seed = (col_area_seed_TA .* n_TA_to_seed) + (col_area_seed_CA .* n_CA_to_seed)
        min_area = covertol * area_to_seed

        # Filter out sites outside of desired depth range
        if .!all(site_data.depth_med .== 0)
            max_depth::Float64 = param_set.depth_min + param_set.depth_offset
            depth_criteria::BitArray{1} = (site_data.depth_med .>= param_set.depth_min) .& (site_data.depth_med .<= max_depth)

            # TODO: Include this change in MATLAB version as well
            if any(depth_criteria .> 0)
                # If sites can be filtered based on depth, do so. Otherwise if no sites can be filtered, remove depth as a criterion.
                depth_priority = depth_priority[depth_criteria]
            else
                @warn "No sites within provided depth range of $(param_set.depth_min) - $(max_depth) meters. Considering all sites."
            end
        end

        # pre-allocate rankings
        rankings = [depth_priority zeros(Int, length(depth_priority)) zeros(Int, length(depth_priority))]

        # Prep site selection
        mcda_vars = DMCDA_vars(
            depth_priority,
            nsiteint,
            sim_params.prioritysites,
            domain.strongpred,
            domain.in_conn,
            domain.out_conn,
            zeros(n_species, n_sites),  # wave stress
            dhw_scen[1, :],  # heat stress
            sum(Y_cover[1, :, :], dims=1),  # sum coral cover
            max_cover,
            total_site_area,
            min_area,
            risktol,
            wtinconnseed,
            wtoutconnseed,
            wtconshade,
            wtwaves,
            wtheat,
            wthicover,
            wtlocover,
            wtpredecseed,
            wtpredecshade
        )
    end

    #### End coral constants

    ## Update ecological parameters based on intervention option
    # Set up assisted adaptation values
    a_adapt = zeros(n_species)
    a_adapt[tabular_enhanced] .= param_set.a_adapt
    a_adapt[corymbose_enhanced] .= param_set.a_adapt

    # Level of natural coral adaptation
    n_adapt = param_set.n_adapt
    bleach_resist = corals.bleach_resist
    bleaching_sensitivity = sim_params.bleaching_sensitivity

    ## Extract other parameters
    LPdhwcoeff = sim_params.LPdhwcoeff # shape parameters relating dhw affecting cover to larval production
    DHWmaxtot = sim_params.DHWmaxtot # max assumed DHW for all scenarios.  Will be obsolete when we move to new, shared inputs for DHW projections
    LPDprm2 = sim_params.LPDprm2 # parameter offsetting LPD curve

    # Wave stress
    mwaves::Array{Float64,3} = zeros(tf, n_species, n_sites)
    wavemort90::Vector{Float64} = corals.wavemort90::Vector{Float64}  # 90th percentile wave mortality

    @inbounds for sp::Int64 in 1:n_species
        @views mwaves[:, sp, :] .= wavemort90[sp] .* wave_scen[:, :, :]
    end

    clamp!(mwaves, 0.0, 1.0)

    Sw_t = 1.0 .- mwaves

    # Flag indicating whether to seed or not to seed
    seed_corals::Bool = (n_TA_to_seed > 0) || (n_CA_to_seed > 0)

    absolute_k_area = vec(total_site_area' .* max_cover)'  # max possible coral area in m^2
    growth::ODEProblem = ODEProblem{true}(growthODE, ode_u, tspan, p)
    @inbounds for tstep::Int64 in 2:tf
        p_step = tstep - 1
        Y_pstep[:, :] .= Y_cover[p_step, :, :]

        sf .= stressed_fecundity(tstep, a_adapt, n_adapt, dhw_scen[p_step, :],
            LPdhwcoeff, DHWmaxtot, LPDprm2, n_groups)

        # Calculates scope for coral fedundity for each size class and at each site.
        fecundity_scope!(fec_scope, fec_all, fec_params_per_m², Y_pstep, total_site_area)

        site_coral_cover = sum(Y_pstep, dims=1)  # dims: nsites * 1
        absolute_site_coral_cover = site_coral_cover .* total_site_area  # in m²
        leftover_space_m² = max.(absolute_k_area .- absolute_site_coral_cover, 0.0)

        area_settled = settler_cover(fec_scope, sf, TP_data, leftover_space_m²,
            sim_params.max_settler_density, sim_params.basal_area_per_settler)

        # Recruitment should represent additional cover, relative to total site area
        # Gets used in ODE
        # p.rec[:, :] .= (area_settled ./ total_site_area)
        p.rec[:, :] .= replace((area_settled ./ absolute_k_area), Inf => 0.0, NaN => 0.0)

        in_shade_years = (shade_start_year <= tstep) && (tstep <= (shade_start_year + shade_years - 1))
        in_seed_years = ((seed_start_year <= tstep) && (tstep <= (seed_start_year + seed_years - 1)))

        @views dhw_t .= dhw_scen[tstep, :]  # subset of DHW for given timestep
        if is_guided && (in_seed_years || in_shade_years)
            # Update dMCDA values
            mcda_vars.heatstressprob .= dhw_t
            mcda_vars.damprob .= @view mwaves[tstep, :, :]
        end

        if is_guided && in_seed_years
            mcda_vars.sumcover .= site_coral_cover
            (prefseedsites, prefshadesites, rankings) = dMCDA(mcda_vars, MCDA_approach,
                seed_decision_years[tstep], shade_decision_years[tstep],
                prefseedsites, prefshadesites, rankings)

            # Log site ranks
            # First col only holds site index ids so skip (with 2:end)
            site_ranks[tstep, rankings[:, 1], :] = rankings[:, 2:end]
        elseif seed_corals && in_seed_years
            # Unguided deployment, seed/shade corals anywhere, so long as available space > 0
            prefseedsites, prefshadesites = unguided_site_selection(prefseedsites, prefshadesites,
                seed_decision_years[tstep], shade_decision_years[tstep],
                nsiteint, vec(leftover_space_m²))

            site_ranks[tstep, prefseedsites, 1] .= 1.0
            site_ranks[tstep, prefshadesites, 2] .= 1.0
        end

        has_shade_sites::Bool = !all(prefshadesites .== 0)
        has_seed_sites::Bool = !all(prefseedsites .== 0)
        if (srm > 0.0) && in_shade_years
            Yshade[tstep, :] .= srm

            # Apply reduction in DHW due to SRM
            adjusted_dhw::Vector{Float64} = max.(0.0, dhw_t .- srm)
        else
            adjusted_dhw = dhw_t
        end

        if (fogging > 0.0) && in_shade_years && (has_seed_sites || has_shade_sites)
            if has_seed_sites
                # Always fog where sites are selected if possible
                site_locs::Vector{Int64} = prefseedsites
            elseif has_shade_sites
                # Otherwise, if no sites are selected, fog selected sites
                site_locs = prefshadesites
            end

            adjusted_dhw[site_locs] .= adjusted_dhw[site_locs] .* (1.0 .- fogging)
            Yfog[tstep, site_locs] .= fogging
        end

        # Calculate and apply bleaching mortality
        bleaching_mortality!(Sbl, tstep, site_data.depth_med, bleaching_sensitivity, adjusted_dhw,
            a_adapt, n_adapt, bleach_resist)

        # Apply seeding
        if seed_corals && in_seed_years && has_seed_sites
            # Calculate proportion to seed based on current available space
            seeded_area = (TA=n_TA_to_seed * col_area_seed_TA, CA=n_CA_to_seed * col_area_seed_CA)

            scaled_seed = distribute_seeded_corals(vec(total_site_area), prefseedsites, vec(leftover_space_m²), seeded_area)

            # Seed each site with TA or CA
            @views Y_pstep[seed_sc_TA, prefseedsites] .= Y_pstep[seed_sc_TA, prefseedsites] .+ scaled_seed.TA
            @views Y_pstep[seed_sc_CA, prefseedsites] .= Y_pstep[seed_sc_CA, prefseedsites] .+ scaled_seed.CA

            # Log seed values/sites (these values are relative to site area)
            Yseed[tstep, 1, prefseedsites] .= scaled_seed.TA
            Yseed[tstep, 2, prefseedsites] .= scaled_seed.CA
        end

        @views prop_loss = Sbl[:, :] .* Sw_t[p_step, :, :]

        # update initial condition
        tmp = ((Y_pstep[:, :] .* prop_loss[:, :]) .* total_site_area) ./ absolute_k_area
        growth.u0[:, :] .= replace(tmp, Inf => 0.0, NaN => 0.0)

        # growth.u0[:, :] .= Y_pstep[:, :] .* prop_loss[:, :]
        sol::ODESolution = solve(growth, solver, save_everystep=false, save_start=false,
            alg_hints=[:nonstiff], abstol=1e-9, reltol=1e-8)  # , adaptive=false, dt=1.0
        # Using the last step from ODE above, proportionally adjust site coral cover
        # if any are above the maximum possible (i.e., the site `k` value)
        Y_cover[tstep, :, :] .= proportional_adjustment!(sol.u[end] .* absolute_k_area ./ total_site_area, cover_tmp, max_cover)
        # Y_cover[tstep, :, :] .= (sol.u[end] .* absolute_k_area) ./ total_site_area

    end

    # Avoid placing importance on sites that were not considered
    # (lower values are higher importance)
    site_ranks[site_ranks.==0.0] .= n_sites + 1

    return (raw=Y_cover, seed_log=Yseed, fog_log=Yfog, shade_log=Yshade, site_ranks=site_ranks)
end
