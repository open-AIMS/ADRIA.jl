using Distributed


"""Scenario running functions"""


"""
    run_scenario(domain::Domain)::NamedTuple

Convenience function to run a scenario for a Domain directly.
"""
function run_scenario(domain::Domain; rep=1)::NamedTuple
    x = (x=nothing, )
    for i in 1:rep
        x = run_scenario(domain.coral_domain, domain.intervention, domain.criteria,
                         domain.coral_params, domain.sim_constants, domain.site_data, 
                         domain.init_coral_cover, domain.coral_domain.ode_p,
                         domain.dhw_scens[:, :, i], domain.wave_scens[:, :, i])
    end

    return x
end


"""
    run_scenario(domain, interv, criteria, corals, sim_params, site_data, init_cov::Array{Float64, 2})::NamedTuple


"""
function run_scenario(domain, interv, criteria, corals, sim_params, site_data, 
                      init_cov::Array{Float64, 2}, p::NamedTuple, dhw_scen::Array, wave_scen::Array)::NamedTuple
    
    # TODO: All cached arrays/values to be moved to outer function and passed in
    # to reduce overall allocations

    tspan::Tuple = (0.0, 1.0)
    solver::BS3 = BS3()
    # solver::Tsit5 = Tsit5()

    # sim constants
    n_sites = domain.n_sites
    nsiteint = sim_params.nsiteint::Int64
    tf = sim_params.tf::Int64
    n_species = domain.n_species
    n_groups = domain.n_groups

    # Gompertz shape parameters for bleaching
    neg_e_p1::Real = -sim_params.gompertz_p1;
    neg_e_p2::Real = -sim_params.gompertz_p2;

    ## TODO: constants for ode_p (should be set outside of run_scenario...)
    @set! p.P = sim_params.max_coral_cover::Float64; # max total coral cover

    # competition factor between Small Massives and Acropora
    @set! p.comp = sim_params.comp::Float64;

    ## END TODO

    # Wave stress
    wavemort90::Vector{Float64} = corals.wavemort90::Vector{Float64};  # 90th percentile wave mortality
    mwaves::Array{Float64, 3} = zeros(tf, n_species, n_sites);
    for sp::Int64 in 1:n_species
        mwaves[:, sp, :] = wavemort90[sp] .* wave_scen;
    end

    mwaves[mwaves .< 0.0] .= 0.0;
    mwaves[mwaves .> 1.0] .= 1.0;

    Sw_t = 1.0 .- wave_scen;

    ## TODO: Define Coral constants
    
    # Define constant table location for seed values
    # Seed1 = Tabular Acropora Enhanced (taxa 1, size class 2)
    # Seed2 = Corymbose Acropora Enhanced (taxa 3, size class 2)
    tabular_enhanced::BitArray = corals.taxa_id .== 1
    corymbose_enhanced::BitArray = corals.taxa_id .== 3
    target_class_id::BitArray = corals.class_id .== 2
    seed_size_class1::Int64 = first(findall(tabular_enhanced .&& target_class_id))
    seed_size_class2::Int64 = first(findall(corymbose_enhanced .&& target_class_id))

    # Update ecological parameters based on intervention option

    # Set up assisted adaptation values
    a_adapt = zeros(n_species);

    # assign level of assisted coral adaptation
    a_adapt[tabular_enhanced] .= interv.a_adapt.val
    a_adapt[corymbose_enhanced] .= interv.a_adapt.val

    # Level of added natural coral adaptation
    n_adapt = corals.n_adapt .+ interv.n_adapt.val;

    #### End coral constants

    site_area::Array{Float64, 2} = site_data.area'
    fec_params::Vector{Float64} = corals.fecundity

    # Data that should be inputs
    potential_settler_cover::Float64 = 9.8175e-08

    # Caches
    TP_data::Array{Float64, 2} = rand(n_sites, n_sites)
    LPs::Array{Float64, 2} = zeros(n_groups, n_sites)
    fec_all::Array{Float64, 2} = similar(init_cov, Float64)
    fec_scope::Array{Float64, 2} = zeros(n_groups, n_sites)
    prop_loss::Array{Float64, 2} = zeros(n_species, n_sites)
    Sbl::Array{Float64, 2} = zeros(n_species, n_sites)
    dhw_step::Vector{Float64} = zeros(n_sites)

    Yout::Array{Float64, 3} = zeros(tf, n_species, n_sites)
    Yout[1, :, :] .= @view init_cov[:, :]
    cov_tmp::Array{Float64, 2} = similar(init_cov, Float64)

    growth::ODEProblem = ODEProblem{true, false}(growthODE, cov_tmp, tspan, p)
    @inbounds for tstep::Int64 in 2:tf
        # TODO: Larval production
        LPs .= rand(n_groups, n_sites)

        p_step = tstep - 1;
        @views cov_tmp[:, :] .= Yout[p_step, :, :];

        # Calculates scope for coral fedundity for each size class and at
        # each site. Now using coral fecundity per m2 in 'coralSpec()'
        fecundity_scope!(fec_scope, fec_all, fec_params, cov_tmp, site_area);

        # adjusting absolute recruitment at each site by dividing by the area
        p.rec .= (potential_settler_cover * ((fec_scope .* LPs) * TP_data)) ./ site_area
        # mul!(p.rec, potential_settler_cover .* (fec_scope .* LPs), TP_data)
        # @. p.rec = p.rec / site_area

        @views dhw_step .= dhw_scen[tstep, :];  # subset of DHW for given timestep

        # TODO: Apply fogging
        # adjusted_dhw::Array{Float64} = dhw_step;

        # Calculate and apply bleaching mortality
        bleaching_mortality!(Sbl, tstep, neg_e_p1,
                             neg_e_p2, a_adapt, n_adapt,
                             corals.bleach_resist, dhw_step);
        # Sbl .= 1.0 .- Sbl

        # proportional loss + proportional recruitment
        @views @. prop_loss = Sbl[:, :] * Sw_t[p_step, :]';

        # @views cov_tmp .= cov_tmp[:, :] .* prop_loss[:, :];

        # Apply seeding
        # Yin1[seed_species1, int_sites] = ...
        # Yin2[seed_species2, int_sites] = ...

        @views growth.u0[:, :] .= cov_tmp[:, :] .* prop_loss[:, :]  # update initial condition
        sol::ODESolution = solve(growth, solver, save_everystep=false, abstol=1e-7, reltol=1e-4)
        @views Yout[tstep, :, :] .= sol.u[end]

        # growth::ODEProblem = ODEProblem{true, false}(growthODE, cov_tmp, tspan, p)
        # sol::ODESolution = solve(growth, solver, save_everystep=false, abstol=1e-6, reltol=1e-7)
        # Yout[tstep, :, :] .= sol.u[end]
    end

    results::NamedTuple = (Y=Yout, )

    return results
end