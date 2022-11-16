using Statistics, Distributed, Logging
import BlackBoxOptim: bboptimize, ParetoFitnessScheme


function robust_optim(dom::Domain; N=16, rcps=["45", "60", "85"])
    interv_names = ADRIA.component_params(dom.model, Intervention).fieldname
    crit_names = ADRIA.component_params(dom.model, Criteria).fieldname
    lever_names = vcat(interv_names, crit_names)

    # TODO FIX: interventions currently include natural adaptation
    # (hang over from MATLAB).
    scens = ADRIA.sample(dom, N)

    # Add RCP column
    scens = repeat(scens, length(rcps))
    insertcols!(scens, 1, :RCP => 45)
    for (subset, rcp) in zip(Iterators.partition(1:nrow(scens), Int(nrow(scens) / length(rcps))), rcps)
        scens[subset, :RCP] .= parse(Int, rcp)
    end

    cache = setup_cache(dom)

    function interv_optim_run(x; dom=dom, scens=scens, rcps=rcps, cache=cache)
        # For a given intervention set, run across sampled combinations 
        # of environmental and biological conditions.

        y = zeros(nrow(scens))
        spec = model_spec(dom)

        # Update scenario spec with given intervention set, ignoring the RCP column
        scens[!, lever_names] .= x'
        scens[!, Not("RCP")] = adjust_samples(dom, spec, scens[!, Not("RCP")])
        with_logger(NullLogger()) do
            @inbounds for i in 1:nrow(scens)
                # update model parameter table with intervention set
                if dom.RCP != string(scens[i, :RCP])
                    dom = switch_RCPs!(dom, string(scens[i, :RCP]))
                end
                update_params!(dom, scens[i, :])

                raw_set = ADRIA.run_direct(dom, cache)

                y[i] = mean(ADRIA.metrics.relative_cover(raw_set.raw))
                # y[i, 2] = mean(ADRIA.metrics.relative_shelter_volume(raw_set.raw))
            end
        end

        # @info y

        # m_rc = mean(ADRIA.metrics.scenario_relative_cover(rs), dims=(:scenarios, :timesteps))[1]
        # m_rsv = mean(ADRIA.metrics.scenario_rsv(rs), dims=(:scenarios, :timesteps))[1]

        return (mean(y), ADRIA.performance.environmental_diversity(spec, scens))
    end

    mod = ADRIA.model_spec(dom)
    bounds = mod[in.(mod[!, :fieldname], [lever_names]), [:lower_bound, :upper_bound]]
    bounds[lever_names.==:guided, :lower_bound] .= 1

    res = bboptimize(interv_optim_run;
        SearchRange=Tuple.(eachrow(bounds)),
        FitnessScheme=ParetoFitnessScheme{2}(is_minimizing=false),
        TraceInterval=60.0,
        Method=:borg_moea,
        # NThreads=Threads.nthreads() - 1,
        MaxTime=60
        # Workers=workers()
    )
    return res
end
