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

    reps = config["operation"]["reps"]

    data_pkg_loc = ARGS[2]
    rcp = ARGS[3]
    scenario_file = ARGS[4]
    scenarios = CSV.read(scenario_file, DataFrame, comment="#")

    # If number of scenarios <= 4, not worth multiprocessing...
    if nrow(scenarios) > 4
        # Use specified number of cores if specified, or all available if given value is <= 0 (the default)
        active_cores = config["operation"]["num_cores"]
        use_cores = (active_cores <= 0) ? num_physical_cores() : active_cores
        addprocs(use_cores; exeflags="--project")

        @eval @everywhere begin
            using Statistics
            using ADRIA
        end
    end

    domain = load_domain(data_pkg_loc, rcp)

    n_scenarios = nrow(scenarios)
    println("ADRIA.jl - Prototype App")
    println("Running $(n_scenarios) scenarios under RCP $(rcp) for $(domain.name)")

    d = run_scenarios(scenarios, domain; reps=reps)

    res = ADRIA.load_results(d)

    println("Results stored in: $(ADRIA.store_location(res))")

    _indicative_result_display(res)
end


function adria_cmd_load()
    res_loc = ARGS[2]
    println("Loading results stored in: $(res_loc)")

    res = ADRIA.load_results(res_loc)
    _indicative_result_display(res)
end


"""
Display results for indicative purposes, just to demonstrate things are working.
Not intended for production.
"""
function _indicative_result_display(res)
    unguided_scens = select(res, "guided .== 0.0")
    guided_scens = select(res, "guided .> 0.0")

    Y_o = ADRIA.metrics.summarize_total_cover(res)

    year_axis = [t % 5 == 0 || t == 2099 ? string(t) : "" for t in 2025:2099]

    # Unguided
    upper = maximum(Y_o.maximum[:, :, :, unguided_scens], dims=(3,4))
    u_d = Array{Float32}(dropdims(upper, dims=(3,4)))
    lower = minimum(Y_o.minimum[:, :, :, unguided_scens], dims=(3,4))
    l_d = Array{Float32}(dropdims(lower, dims=(3,4)))

    p = plot(u_d, fillrange=l_d, color=:lightsalmon1, alpha=0.7, label="")
    mean_unguided = dropdims(mean(Y_o.mean[:, :, :, unguided_scens], dims=(3,4)), dims=(3,4))
    plot!(mean_unguided, label="Unguided mean", linecolor=:red, alpha=0.7)

    # Guided
    upper = maximum(Y_o.maximum[:, :, :, guided_scens], dims=(3,4))
    u_d = Array{Float32}(dropdims(upper, dims=(3,4)))
    lower = minimum(Y_o.minimum[:, :, :, guided_scens], dims=(3,4))
    l_d = Array{Float32}(dropdims(lower, dims=(3,4)))

    plot!(u_d, fillrange=l_d, color=:lightblue2, alpha=0.7, label="")
    mean_guided = dropdims(mean(Y_o.mean[:, :, :, guided_scens], dims=(3,4)), dims=(3,4))
    plot!(mean_guided, 
        label="Guided mean", linecolor=:blue, alpha=0.7, 
        xlabel="Year", ylabel="Relative Cover", 
        xticks=(1:75, year_axis))

    p2 = plot(mean_guided - mean_unguided, label="Difference (Mean Guided - Unguided)", 
            xlabel="Year", ylabel="δ Relative Cover",
            xticks=(1:75, year_axis)
    )

    fig = plot(p, p2, size=(1000, 500), layout=(1,2), left_margin=5mm, bottom_margin=5mm, xrotation=45, 
            legend=:best, fg_legend=:transparent, bg_legend=:transparent)
    # display(fig)
    # gui(fig)

    savefig("$(ADRIA.store_name(res)).png")

    # TODO: Force display from commandline
    # https://discourse.julialang.org/t/how-to-display-the-plots-by-executing-the-file-from-command-line/13822/2
end
