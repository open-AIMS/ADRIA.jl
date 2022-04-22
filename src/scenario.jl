"""Scenario running functions"""


"""
    run_scenario(domain::Domain; reps=1, MCDA_approach=1)::NamedTuple

Convenience function to directly run a scenario for a Domain with pre-set values.
"""
function run_scenario(domain::Domain; reps=1, MCDA_approach=1)::NamedTuple

    # Set scenario constants here to avoid repeated allocations
    # TODO: Set all constants outside rep or scenario loop
    @set! domain.coral_growth.ode_p.P = domain.sim_constants.max_coral_cover::Float64
    @set! domain.coral_growth.ode_p.comp = domain.sim_constants.comp::Float64

    # TODO: Actually store all results
    for i in 1:reps
        x = run_scenario(domain.coral_growth, domain.intervention, domain.criteria,
            domain.coral_params, domain.sim_constants, domain.site_data,
            domain.init_coral_cover, domain.coral_growth.ode_p,
            domain.dhw_scens[:, :, i], domain.wave_scens[:, :, i], MCDA_approach)
    end

    return x
end


"""
    run_scenario(domain, interv, criteria, corals, sim_params, site_data, init_cov::Array{Float64, 2}, p::NamedTuple, 
                 dhw_scen::Array, wave_scen::Array, MCDA_approach::Int64)::NamedTuple

Core scenario running function.
"""
function run_scenario(domain, interv, criteria, corals, sim_params, site_data,
    init_cov::Array{Float64,2}, p::NamedTuple, dhw_scen::Array,
    wave_scen::Array, MCDA_approach::Int64)::NamedTuple

    # TODO: All cached arrays/values to be moved to outer function and passed in
    # to reduce overall allocations (e.g., sim constants don't change across all scenarios)

    tspan::Tuple = (0.0, 1.0)
    solver::BS3 = BS3()

    # sim constants
    n_sites = domain.n_sites
    nsiteint = sim_params.nsiteint::Int64
    tf = sim_params.tf::Int64
    n_species = domain.n_species
    n_groups = domain.n_groups

    # years to start seeding/shading
    seed_start_year = interv.seed_year_start.val
    shade_start_year = interv.shade_year_start.val

    seed_TA_vol = interv.seed_TA.val # tabular Acropora size class 2, per year per species per cluster
    seed_CA_vol = interv.seed_CA.val # corymbose Acropora size class 2, per year per species per cluster
    fogging = interv.fogging.val #  percent reduction in bleaching mortality through fogging
    srm = interv.SRM.val #  DHW equivalents reduced by some shading mechanism
    seed_years = interv.seed_years.val # number of years to seed
    shade_years = interv.shade_years.val # number of years to shade

    # Gompertz shape parameters for bleaching
    neg_e_p1::Real = -sim_params.gompertz_p1
    neg_e_p2::Real = -sim_params.gompertz_p2

    ## END TODO

    site_area::Array{Float64,2} = site_data.area'

    fec_params::Vector{Float64} = corals.fecundity
    potential_settler_cover::Float64 = (sim_params.max_settler_density *
                                        sim_params.basal_area_per_settler *
                                        sim_params.density_ratio_of_settlers_to_larvae)

    # Caches
    TP_data::Array{Float64,2} = rand(n_sites, n_sites)
    LPs::Array{Float64,2} = zeros(n_groups, n_sites)
    fec_all::Array{Float64,2} = similar(init_cov, Float64)
    fec_scope::Array{Float64,2} = zeros(n_groups, n_sites)
    prop_loss::Array{Float64,2} = zeros(n_species, n_sites)
    Sbl::Array{Float64,2} = zeros(n_species, n_sites)
    dhw_step::Vector{Float64} = zeros(n_sites)

    Yout::Array{Float64,3} = zeros(tf, n_species, n_sites)
    Yout[1, :, :] .= @view init_cov[:, :]
    cov_tmp::Array{Float64,2} = similar(init_cov, Float64)

    # Yshade = ndSparse(zeros(tf, nsites));
    # Yfog = ndSparse(zeros(tf, nsites));
    # Yseed = ndSparse(zeros(tf, 2, nsites)); % only log seedable corals to save memory

    site_rankings = SparseArray(zeros(tf, n_sites, 2)) # log seeding/shading ranks
    Yshade = SparseArray(spzeros(tf, n_sites))
    Yfog = SparseArray(spzeros(tf, n_sites))
    Yseed = SparseArray(zeros(tf, 2, n_sites))

    # Intervention strategy: 0 is random, 1 is guided
    strategy = interv.guided.val
    is_guided = strategy > 0

    # Years at which to reassess seeding site selection
    seed_decision_years = repeat([false], tf)
    shade_decision_years = repeat([false], tf)

    if interv.seed_freq.val > 0
        shade_decision_years[seed_start_year:interv.seed_freq.val:(seed_start_year+seed_years-1)] .= true
    else
        shade_decision_years[max(seed_start_year, 2)] = true
    end

    if interv.shade_freq.val > 0
        shade_decision_years[shade_start_year:interv.shade_freq.val:(shade_start_year+shade_years-1)] .= true
    else
        shade_decision_years[max(shade_start_year, 2)] = true
    end

    prefseedsites = false
    prefshadesites = false

    if is_guided
        ## Weights for connectivity , waves (ww), high cover (whc) and low
        wtwaves = criteria.wave_stress.val # weight of wave damage in MCDA
        wtheat = criteria.heat_stress.val # weight of heat damage in MCDA
        wtconshade = criteria.shade_connectivity.val # weight of connectivity for shading in MCDA
        wtconseed = criteria.seed_connectivity.val # weight of connectivity for seeding in MCDA
        wthicover = criteria.coral_cover_high.val # weight of high coral cover in MCDA (high cover gives preference for seeding corals but high for SRM)
        wtlocover = criteria.coral_cover_low.val # weight of low coral cover in MCDA (low cover gives preference for seeding corals but high for SRM)
        wtpredecseed = criteria.seed_priority.val # weight for the importance of seeding sites that are predecessors of priority reefs
        wtpredecshade = criteria.shade_priority.val # weight for the importance of shading sites that are predecessors of priority reefs
        risktol = criteria.deployed_coral_risk_tol.val # risk tolerance

        # Filter out sites outside of desired depth range
        if .!all(site_data.sitedepth .== 0)
            max_depth = criteria.depth_min.val + criteria.depth_offset.val
            depth_criteria = (site_data.sitedepth .> -max_depth) .& (site_data.sitedepth .< -criteria.depth_min)
            depth_priority = site_data[depth_criteria, domain.site_id_col]
        else
            # No depth data, so consider all sites
            depth_priority = site_data[:, domain.site_id_col]
        end

        # if any(isa.(depth_priority, String))
        #     # Catch edge case where IDs are interpreted as text/cells
        #     depth_priority = 1:length(depth_priority)';
        #     # depth_priority = depth_priority';
        # end

        max_cover = site_data.k / 100.0 # Max coral cover at each site. Divided by 100 to convert to proportion

        # pre-allocate prefseedsites, prefshadesites and rankings
        rankings = [depth_priority, zeros(length(depth_priority), 1), zeros(length(depth_priority), 1)]
        prefseedsites = zeros(1, nsiteint)
        prefshadesites = zeros(1, nsiteint)

        # Prep site selection
        mcda_vars = DMCDA_vars(
            depth_priority,
            nsiteint,
            sim_params.prioritysites,
            strongpred,
            domain.site_ranks.C1,
            0,  # dam prob
            0,  # heatstressprob
            0,  # sumcover
            max_cover,
            site_data.area,
            risktol,
            wtconseed,
            wtconshade,
            wtwaves,
            wtheat,
            wthicover,
            wtlocover,
            wtpredecseed,
            wtpredecshade
        )
    else
        Random.seed!(floor(Int, sum(Model(interv)[:val]) + sum(Model(criteria)[:val])))
    end

    # n_species

    # containers for seeding, shading and cooling
    nprefseed = zeros(tf)
    nprefshade = zeros(tf)

    # Define constant table location for seed values
    # Seed1 = Tabular Acropora Enhanced (taxa 1, size class 2)
    # Seed2 = Corymbose Acropora Enhanced (taxa 3, size class 2)
    tabular_enhanced::BitArray = corals.taxa_id .== 1
    corymbose_enhanced::BitArray = corals.taxa_id .== 3
    target_class_id::BitArray = corals.class_id .== 2
    seed_size_class1::Int64 = first(findall(tabular_enhanced .& target_class_id))
    seed_size_class2::Int64 = first(findall(corymbose_enhanced .& target_class_id))

    #### End coral constants

    ## Update ecological parameters based on intervention option
    # Set up assisted adaptation values
    a_adapt = zeros(n_species)

    # assign level of assisted coral adaptation
    a_adapt[tabular_enhanced] .= interv.a_adapt.val
    a_adapt[corymbose_enhanced] .= interv.a_adapt.val
    n_adapt = interv.n_adapt.val  # Level of natural coral adaptation

    # level of added natural coral adaptation
    n_adapt = interv.n_adapt.val  # natad = coral_params.natad + interv.Natad;
    bleach_resist = corals.bleach_resist

    ## Extract other parameters
    LPdhwcoeff = sim_params.LPdhwcoeff # shape parameters relating dhw affecting cover to larval production
    DHWmaxtot = sim_params.DHWmaxtot # max assumed DHW for all scenarios.  Will be obsolete when we move to new, shared inputs for DHW projections
    LPDprm2 = sim_params.LPDprm2 # parameter offsetting LPD curve

    # Wave stress
    mwaves::Array{Float64,3} = zeros(tf, n_species, n_sites)
    wavemort90::Vector{Float64} = corals.wavemort90::Vector{Float64}  # 90th percentile wave mortality
    @inbounds for sp::Int64 in 1:n_species
        mwaves[:, sp, :] .= wavemort90[sp] .* wave_scen
    end

    mwaves[mwaves.<0.0] .= 0.0
    mwaves[mwaves.>1.0] .= 1.0

    Sw_t = 1.0 .- mwaves

    # Flag indicating whether to seed or not to seed
    seed_corals = (seed_TA_vol > 0) || (seed_CA_vol > 0)

    # extract colony areas for sites selected and convert to m^2
    col_area_seed_TA = corals.colony_area_cm2[seed_size_class1] / 10^4
    col_area_seed_CA = corals.colony_area_cm2[seed_size_class2] / 10^4

    growth::ODEProblem = ODEProblem{true,false}(growthODE, cov_tmp, tspan, p)
    @inbounds for tstep::Int64 in 2:tf
        p_step = tstep - 1
        @views cov_tmp[:, :] .= Yout[p_step, :, :]

        LPs .= larval_production(tstep, a_adapt, n_adapt, dhw_scen[p_step, :],
            LPdhwcoeff, DHWmaxtot, LPDprm2, n_groups)

        # Calculates scope for coral fedundity for each size class and at
        # each site. Now using coral fecundity per m2 in 'coralSpec()'
        fecundity_scope!(fec_scope, fec_all, fec_params, cov_tmp, site_area)

        # adjusting absolute recruitment at each site by dividing by the area
        p.rec .= (potential_settler_cover * ((fec_scope .* LPs) * TP_data)) / site_area
        # mul!(p.rec, potential_settler_cover .* (fec_scope .* LPs), TP_data)
        # @. p.rec = p.rec / site_area

        @views dhw_step .= dhw_scen[tstep, :]  # subset of DHW for given timestep

        in_shade_years = (shade_start_year <= tstep) && (tstep <= (shade_start_year + shade_years - 1))
        in_seed_years = ((seed_start_year <= tstep) && (tstep <= (seed_start_year + seed_years - 1)))
        if is_guided
            # Update dMCDA values
            dMCDA_vars.damprob .= dropdims(mwaves(tstep, :, :), dims=1)'
            dMCDA_vars.heatstressprob .= dhw_step'

            dMCDA_vars.sumcover = sum(Y_pstep, dims=1)' # dims: nsites * 1

            (prefseedsites, prefshadesites, nprefseedsites, nprefshadesites, rankings) = dMCDA(mcda_vars, MCDA_approach,
                seed_decision_years[tstep], shade_decision_years[tstep],
                refseedsites, prefshadesites, rankings)

            nprefseed[tstep] = nprefseedsites    # number of preferred seeding sites
            nprefshade[tstep] = nprefshadesites  # number of preferred shading sites

            # Log site ranks
            # First col only holds site ids so skip (with 2:end)
            site_rankings[tstep, rankings[:, 1], :] = rankings[:, 2:end]
        else
            # Unguided
            if seed_decision_years[tstep]
                # Unguided deployment, seed/shade corals anywhere
                prefseedsites = rand(1:n_sites, nsiteint)
            end

            if shade_decision_years[tstep]
                prefshadesites = rand(1:n_sites, nsiteint)
            end
        end

        has_shade_sites = !all(prefshadesites .== 0)
        has_seed_sites = !all(prefseedsites .== 0)
        if (srm > 0) && in_shade_years && has_shade_sites
            Yshade[tstep, prefshadesites] = srm

            # Apply reduction in DHW due to shading
            adjusted_dhw::Array{Float64} = max(0.0, dhw_step - Yshade(tstep, :))
        else
            adjusted_dhw = dhw_step
        end

        if (fogging > 0.0) && in_shade_years && (has_seed_sites || has_shade_sites)
            if has_seed_sites
                # Always fog where sites are selected if possible
                adjusted_dhw[prefseedsites] = adjusted_dhw[prefseedsites] .* (1.0 - fogging)
                Yfog[tstep, prefseedsites] .= fogging
            elseif has_shade_sites
                # Otherwise, if no sites are selected, fog selected shade sites
                adjusted_dhw[prefshadesites] = adjusted_dhw[prefshadesites] .* (1.0 - fogging)
                Yfog[tstep, prefshadesites] .= fogging
            end
        end

        # Calculate and apply bleaching mortality
        bleaching_mortality!(Sbl, tstep, neg_e_p1,
            neg_e_p2, a_adapt, n_adapt,
            bleach_resist, adjusted_dhw)

        # proportional loss + proportional recruitment
        @views prop_loss = Sbl[:, :] .* Sw_t[p_step, :, :]
        @views @. cov_tmp = cov_tmp[:, :] * prop_loss[:, :]

        # Apply seeding
        if seed_corals && in_seed_years && has_seed_sites
            # extract site area for sites selected and scale by available space for populations (k/100)
            site_area_seed = site_area[prefseedsites] .* max_cover[prefseedsites]

            scaled_seed_TA = (((seed_TA_vol / nsiteint) * col_area_seed_TA) ./ site_area_seed)'
            scaled_seed_CA = (((seed_CA_vol / nsiteint) * col_area_seed_CA) ./ site_area_seed)'

            # Seed each site with the value indicated with seed1/seed2
            cov_tmp[seed_size_class1, prefseedsites] .= cov_tmp[seed_size_class1, prefseedsites] .+ scaled_seed_TA  # seed Enhanced Tabular Acropora
            cov_tmp[seed_size_class2, prefseedsites] .= cov_tmp[seed_size_class2, prefseedsites] .+ scaled_seed_CA  # seed Enhanced Corymbose Acropora

            # Log seed values/sites
            Yseed[tstep, 1, prefseedsites] = scaled_seed1  # log site as seeded with Enhanced Tabular Acropora
            Yseed[tstep, 2, prefseedsites] = scaled_seed2  # log site as seeded with Enhanced Corymbose Acropora
        end

        growth.u0[:, :] .= @views cov_tmp[:, :] .* prop_loss[:, :]  # update initial condition
        sol::ODESolution = solve(growth, solver, save_everystep=false, abstol=1e-7, reltol=1e-4)
        @views Yout[tstep, :, :] .= sol.u[end]

        # growth::ODEProblem = ODEProblem{true, false}(growthODE, cov_tmp, tspan, p)
        # sol::ODESolution = solve(growth, solver, save_everystep=false, abstol=1e-6, reltol=1e-7)
        # Yout[tstep, :, :] .= sol.u[end]
    end

    results::NamedTuple = (Y=Yout, seed_log=Yseed, fog_log=Yfog, shade_log=Yshade, site_rankings=site_rankings)

    return results
end