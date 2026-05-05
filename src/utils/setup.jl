"""
    setup()::Nothing

Initialize ADRIA configuration options from `config.toml`
or load defaults if not found.
"""
function setup()::Nothing
    if has_setup()
        return nothing
    end

    try
        # Load in configuration settings
        config = TOML.parsefile(joinpath(pwd(), "config.toml"))
        config_operation = config["operation"]
        ENV["ADRIA_OUTPUT_DIR"] = config["results"]["output_dir"]
        ENV["ADRIA_NUM_CORES"] = config_operation["num_cores"]
        ENV["ADRIA_THRESHOLD"] = config_operation["threshold"]
        ENV["ADRIA_DEBUG"] =
            haskey(config_operation, "debug") ? config_operation["debug"] : false

        ENV["ADRIA_LOG_DHW_TOLS"] =
            haskey(config_operation, "log_dhw_tols") ?
            config_operation["log_dhw_tols"] :
            false

        ENV["ADRIA_LOG_COVER"] =
            haskey(config_operation, "log_cover") ? config_operation["log_cover"] : false

        ENV["ADRIA_RNG_SEED"] =
            haskey(config_operation, "rng_seed") ? config_operation["rng_seed"] : false
    catch
        @warn "Could not find config.toml file.\nApplying default configuration and saving results to 'Outputs' in current directory."

        # Note: anything stored in ENV will be stored as String.
        ENV["ADRIA_OUTPUT_DIR"] = "./Outputs"
        ENV["ADRIA_NUM_CORES"] = 1
        ENV["ADRIA_THRESHOLD"] = Float32(1e-8)
        ENV["ADRIA_DEBUG"] = false
        ENV["ADRIA_LOG_DHW_TOLS"] = false
        ENV["ADRIA_LOG_COVER"] = false
    end

    return nothing
end

"""Check to ensure setup has been run."""
function has_setup()::Bool
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
        active_cores::Int64 = parse(Int64, ENV["ADRIA_NUM_CORES"])
        if active_cores <= 0
            active_cores = cpucores()
        end

        if active_cores > 1
            addprocs(active_cores; exeflags="--project=$(Base.active_project())")
        end
    end

    return nothing
end

"""Remove workers and free up memory."""
function _remove_workers()::Nothing
    if nworkers() > 1 || nprocs() > 1
        rmprocs(workers())
    end

    return nothing
end

function is_test_env()::Bool
    return get(ENV, "ADRIA_TEST", "false") == "true"
end

"""
    set_random_seed(param_set)::Nothing

Set the RNG seed from `ENV["ADRIA_RNG_SEED"]`. If unset or `"false"`, derives the seed
from the sum of all non-RCP parameters in `param_set`.

# Arguments
- `param_set` : YAXArray of scenario parameter values for a single run.
"""
function set_random_seed(param_set::YAXArray)::AbstractRNG
    rnd_seed_val::Int64 =
        get(ENV, "ADRIA_RNG_SEED", "false") == "false" ?
        floor(Int64, sum(param_set[Where(x -> x != "RCP")])) : # select everything except RCP
        parse(Int64, ENV["ADRIA_RNG_SEED"])
    return Xoshiro(rnd_seed_val)
end
