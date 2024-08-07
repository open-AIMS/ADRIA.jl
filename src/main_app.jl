using Plots, Plots.Measures

const valid_methods = ["run", "load", "help"]

"""

Main entry point for ADRIA application.
"""
function julia_main()::Cint
    method = ARGS[1]
    if method == "run"
        adria_cmd_run()
    elseif method == "load"
        adria_cmd_load()
    elseif method == "help"
        msg = """
        Available commands:
        run :  Run scenarios defined in a CSV file, storing results in the output directory specified in config.toml
               ADRIA run [data package location] [RCP scenario] [scenario file]

        load : Load and display results of previously run scenarios.
               ADRIA load [result location]

        help : Display this help list.
        """
        println(msg)

        return 0
    else
        error("Unknown command $(method), valid commands are: $(valid_methods)")
        return 1
    end

    return 0  # if things finished successfully
end

function adria_cmd_run()
    config = TOML.parsefile("config.toml")
    ENV["ADRIA_OUTPUT_DIR"] = config["results"]["output_dir"]

    data_pkg_loc = ARGS[2]
    rcp = ARGS[3]
    scenario_file = ARGS[4]
    scenarios = CSV.read(scenario_file, DataFrame; comment="#")

    # If number of scenarios <= 4, not worth multiprocessing...
    if nrow(scenarios) > 4
        # Use specified number of cores if specified, or all available if given value is <= 0 (the default)
        active_cores = config["operation"]["num_cores"]
        use_cores = (active_cores <= 0) ? num_physical_cores() : active_cores
        addprocs(use_cores; exeflags="--project")

        @eval @everywhere using ADRIA
    end

    domain = load_domain(data_pkg_loc, rcp)

    n_scenarios = nrow(scenarios)
    println("ADRIA.jl - Prototype App")
    println("Running $(n_scenarios) scenarios under RCP $(rcp) for $(domain.name)")

    d = run_scenarios(domain, scenarios)

    res = ADRIA.load_results(d)

    println("Results stored in: $(ADRIA.result_location(res))")

    _indicative_result_display(res)
    return nothing
end

function adria_cmd_load()
    res_loc = ARGS[2]
    println("Loading results stored in: $(res_loc)")

    res = ADRIA.load_results(res_loc)
    _indicative_result_display(res)
    return nothing
end

"""
Display results for indicative purposes, just to demonstrate things are working.
Not intended for production.
"""
function _indicative_result_display(res)
    nodeploy_scens = findall(select(res, "guided .== 0") .&& select(res, "seed_TA .== 0"))
    unguided_scens = findall(select(res, "guided .== 0") .&& select(res, "seed_TA .> 0"))
    guided_scens = findall(select(res, "guided .> 0") .&& select(res, "seed_TA .> 0"))

    Y_no = ADRIA.metrics.summarize_relative_cover(selectdim(res.raw, 5, nodeploy_scens))
    Y_ung = ADRIA.metrics.summarize_relative_cover(selectdim(res.raw, 5, unguided_scens))
    Y_g = ADRIA.metrics.summarize_relative_cover(selectdim(res.raw, 5, guided_scens))

    year_axis = [t % 5 == 0 || t == 2099 ? string(t) : "" for t in 2025:2099]

    # No deployment
    upper = Y_no[:upper_95]
    lower = Y_no[:lower_95]

    p = plot(upper; fillrange=lower, color=:lightsalmon1, alpha=0.8, label="")
    plot!(Y_no[:median]; label="No Deployment median", linecolor=:red, alpha=0.8)

    # Unguided Deployment
    upper = Y_ung[:upper_95]
    lower = Y_ung[:lower_95]
    p = plot!(upper; fillrange=lower, color=:lightblue2, alpha=0.8, label="")
    plot!(Y_ung[:median]; label="Unguided median", linecolor=:blue, alpha=0.5)

    # Guided
    upper = Y_g[:upper_95]
    lower = Y_g[:lower_95]

    plot!(upper; fillrange=lower, color=:lightseagreen, alpha=0.4, label="")
    plot!(Y_g[:median];
        label="Guided median", linecolor=:green, alpha=0.4,
        xlabel="Year", ylabel="Relative Cover",
        xticks=(1:75, year_axis))

    p2 = plot(Y_ung[:mean] - Y_no[:mean]; label="Guided - No Deployment (μ)",
        xlabel="Year", ylabel="δ Relative Cover",
        xticks=(1:75, year_axis), color=:red
    )
    plot!(Y_g[:mean] - Y_ung[:mean]; label="Guided - Unguided (μ)", color=:blue)

    fig = plot(p, p2; size=(1000, 500), layout=(1, 2), left_margin=5mm, bottom_margin=5mm,
        xrotation=45,
        legend=:best, fg_legend=:transparent, bg_legend=:transparent)
    # display(fig)
    # gui(fig)

    savefig(joinpath(ENV["ADRIA_OUTPUT_DIR"], "$(ADRIA.store_name(res)).png"))

    # TODO: Force display from commandline
    # https://discourse.julialang.org/t/how-to-display-the-plots-by-executing-the-file-from-command-line/13822/2
    return nothing
end
