using Setfield

import ADRIA: to_spec
import ModelParameters: Model


function run_scenario(domain, interv, criteria, corals, sim_params, site_data, init_cov::Array{Float64, 2})::NamedTuple
    p::NamedTuple = domain.ode_p
    tspan::Tuple = (0.0, 1.0)
    solver::BS3 = BS3()

    # sim constants
    nsiteint = sim_params.nsiteint
    tf = sim_params.tf

    # Gompertz shape parameters for bleaching
    # neg_e_p1 = -sim_params.gompertz_p1;
    # neg_e_p2 = -sim_params.gompertz_p2;

    ## TODO: constants for ode_p (should be set outside of run_scenario...)
    @set! p.P = sim_params.max_coral_cover; # max total coral cover

    # competition factor between Small Massives and Acropora
    @set! p.comp = sim_params.comp;

    ## END TODO

    # Data that should be inputs
    site_area::Array{Float64, 2} = rand(1, 561)
    # coral_params = rand(36, 13)
    fec_params::Array{Float64, 1} = corals.fecundity
    potential_settler_cover::Float64 = 9.8175e-08

    # Caches
    TP_data::Array{Float64, 2} = rand(561, 561)
    LPs::Array{Float64, 2} = zeros(6, 561)
    fec_all::Array{Float64, 2} = similar(init_cov, Float64)
    fec_scope::Array{Float64, 2} = zeros(6, 561)

    Yout::Array{Float64, 3} = zeros(75, 36, 561)
    Yout[1, :, :] .= init_cov
    cov_tmp::Array{Float64, 2} = similar(init_cov, Float64)

    # growth::ODEProblem = ODEProblem{true, false}(growthODE, cov_tmp, tspan, p)
    @inbounds for tstep in 2:75
        LPs .= rand(6, 561)

        @views cov_tmp .= Yout[tstep-1, :, :];

        # Calculates scope for coral fedundity for each size class and at
        # each site. Now using coral fecundity per m2 in 'coralSpec()'
        fecundity_scope!(fec_scope, fec_all, fec_params, cov_tmp, site_area);

        # adjusting absolute recruitment at each site by dividing by the area
        # p.rec .= (potential_settler_cover * ((fec_scope .* LPs) * TP_data)) ./ site_area
        mul!(p.rec, potential_settler_cover .* (fec_scope .* LPs), TP_data)
        p.rec .= p.rec ./ site_area

        # Apply fogging
        # adjusted_dhw = ...

        # Calculate and apply bleaching mortality
        # Sbl = 1.0 - bleaching_mortality(tstep, neg_e_p1,
        #                 neg_e_p2, assistadapt, natad,
        #                 bleach_resist, adjusted_dhw);

        # # proportional loss + proportional recruitment
        # prop_loss = Sbl .* squeeze(Sw_t(p_step, :, :));

        # cov_tmp .= cov_tmp .* prop_loss;

        # Apply seeding
        # Yin1[seed1, int_sites] = ...
        # Yin2[seed2, int_sites] = ...
        
        # @views growth.u0 .= cov_tmp[:, :]  # update initial condition
        # sol::ODESolution = solve(growth, solver, save_everystep=false, abstol=1e-6, reltol=1e-7)
        # Yout[tstep, :, :] .= sol.u[end]

        growth::ODEProblem = ODEProblem{true, false}(growthODE, cov_tmp, tspan, p)
        sol::ODESolution = solve(growth, solver, save_everystep=false, abstol=1e-6, reltol=1e-7)
        Yout[tstep, :, :] .= sol.u[end]
    end

    results::NamedTuple = (Y=Yout, )

    return results
end