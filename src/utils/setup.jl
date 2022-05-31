function setup()
    try
        # Load in configuration settings
        config = TOML.parsefile(joinpath(pwd(), "config.toml"))
        ENV["ADRIA_OUTPUT_DIR"] = config["results"]["output_dir"]
        ENV["ADRIA_NUM_CORES"] = config["operation"]["num_cores"]
        ENV["ADRIA_reps"] = config["operation"]["reps"]
        ENV["ADRIA_THRESHOLD"] = config["operation"]["threshold"]
    catch
        @warn "Could not find config.toml file.\nApplying default configuration and saving results to 'Outputs' in current directory."

        # Note: anything stored in ENV will be stored as String.
        ENV["ADRIA_OUTPUT_DIR"] = "./Outputs"
        ENV["ADRIA_NUM_CORES"] = -1
        ENV["ADRIA_reps"] = 20
        ENV["ADRIA_THRESHOLD"] = Float32(1e-8)
    end


    # Spin up workers if needed
    if nprocs() == 1
        active_cores = parse(Int, ENV["ADRIA_NUM_CORES"])
        if active_cores <= 0
            active_cores = cpucores()
        end

        addprocs(active_cores, exeflags="--project")
    end
end


"""Check to ensure setup has been run."""
function has_setup()
    try
        ENV["ADRIA_OUTPUT_DIR"]
    catch err
        if isa(err, KeyError)
            error("Setup has not been run.")
        else
            rethrow(err)
        end
    end
end