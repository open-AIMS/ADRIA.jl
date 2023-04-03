"""
    setup()::Nothing

Initialize ADRIA configuration options from `config.toml` 
or load defaults if not found.
"""
function setup()::Nothing
    if has_setup()
        return
    end

    try
        # Load in configuration settings
        config = TOML.parsefile(joinpath(pwd(), "config.toml"))
        ENV["ADRIA_OUTPUT_DIR"] = config["results"]["output_dir"]
        ENV["ADRIA_NUM_CORES"] = config["operation"]["num_cores"]
        ENV["ADRIA_THRESHOLD"] = config["operation"]["threshold"]
        ENV["ADRIA_DEBUG"] = haskey(config["operation"], "debug") ? config["operation"]["debug"] : false
    catch
        @warn "Could not find config.toml file.\nApplying default configuration and saving results to 'Outputs' in current directory."

        # Note: anything stored in ENV will be stored as String.
        ENV["ADRIA_OUTPUT_DIR"] = "./Outputs"
        ENV["ADRIA_NUM_CORES"] = 1
        ENV["ADRIA_THRESHOLD"] = Float32(1e-8)
        ENV["ADRIA_DEBUG"] = false
    end

    return
end


"""Check to ensure setup has been run."""
function has_setup()
    try
        ENV["ADRIA_OUTPUT_DIR"]
    catch err
        if !isa(err, KeyError)
            rethrow(err)
        end

        return false
    end

    return true
end

"""Spin up workers if needed."""
function _setup_workers()::Nothing
    if nprocs() == 1 && (parse(Bool, ENV["ADRIA_DEBUG"]) == false)
        active_cores = parse(Int, ENV["ADRIA_NUM_CORES"])
        if active_cores <= 0
            active_cores = cpucores()
        end

        if active_cores > 1
            addprocs(active_cores; exeflags="--project=$(Base.active_project())")
        end
    end

    return
end

"""Remove workers and free up memory."""
function _remove_workers()::Nothing
    if nprocs() > 1
        rmprocs(workers()...)
    end

    return
end

