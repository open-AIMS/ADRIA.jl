using Dates
using SHA
using TOML

Base.@kwdef struct CalibrationRunConfig
    mode::String
    seed::Int
    max_steps::Int
    max_evals::Int
    max_time::Float64
    method::Symbol
    run_id::String
    resume::Bool
    output_dir::String
end

function calibration_run_config(repo_root::AbstractString)::CalibrationRunConfig
    mode = uppercase(get(ENV, "BBO_MODE", "FOCUSED"))
    seed = parse(Int, get(ENV, "BBO_SEED", "20260928"))
    max_steps = parse(Int, get(ENV, "BBO_MAX_STEPS", "20"))
    max_evals = parse(Int, get(ENV, "BBO_MAX_EVALS", "0"))
    max_time = parse(Float64, get(ENV, "BBO_MAX_TIME", "0.0"))
    method = Symbol(get(ENV, "BBO_METHOD", "adaptive_de_rand_1_bin_radiuslimited"))
    resume = lowercase(get(ENV, "BBO_RESUME", "false")) == "true"
    default_id = Dates.format(now(), "yyyymmddTHHMMSS") * "_" * lowercase(mode) * "_seed" * string(seed)
    run_id = get(ENV, "BBO_RUN_ID", default_id)
    occursin(r"^[A-Za-z0-9_.-]+$", run_id) || error("BBO_RUN_ID contains unsupported characters")
    output_root = get(ENV, "BBO_OUTPUT_ROOT", joinpath(repo_root, "sandbox", "calibration", "runs"))
    output_dir = normpath(joinpath(output_root, run_id))
    mkpath(output_dir)
    return CalibrationRunConfig(; mode, seed, max_steps, max_evals, max_time, method, run_id, resume, output_dir)
end

function git_revision(path::AbstractString)::String
    try
        return readchomp(`git -C $path rev-parse HEAD`)
    catch
        return "unknown"
    end
end

file_sha256(path::AbstractString)::String = isfile(path) ? bytes2hex(sha256(read(path))) : "missing"

function write_run_metadata(
    path::AbstractString,
    config::CalibrationRunConfig,
    repo_root::AbstractString,
    param_names,
    search_space,
    input_paths::Vector{String}
)
    metadata = Dict(
        "run" => Dict(
            "run_id" => config.run_id,
            "mode" => config.mode,
            "seed" => config.seed,
            "max_steps" => config.max_steps,
            "max_evals" => config.max_evals,
            "max_time_seconds" => config.max_time,
            "method" => string(config.method),
            "resume" => config.resume,
            "started_at" => string(now()),
            "julia_version" => string(VERSION),
        ),
        "revision" => Dict(
            "adria_repository" => git_revision(repo_root),
            "cotsmod_repository" => git_revision(normpath(joinpath(repo_root, "..", "COTSMod.jl"))),
        ),
        "search" => Dict(
            "parameter_names" => collect(param_names),
            "lower_bounds" => [bounds[1] for bounds in search_space],
            "upper_bounds" => [bounds[2] for bounds in search_space],
        ),
        "inputs" => Dict(path => file_sha256(path) for path in input_paths),
        "environment" => Dict(
            key => get(ENV, key, "") for key in [
                "COTS_EXTERNAL_PULSE", "COTS_PULSE_START", "COTS_PULSE_DURATION",
                "COTS_PULSE_REPEAT_INTERVAL", "COTS_PULSE_RELATIVE_MAGNITUDE",
                "COTS_SEED_FIRST_N", "COTS_INITIAL_MULTIPLIER",
                "ADRIA_COTS_CONNECTIVITY_MODE", "BBO_EXCLUDE_REEFS",
                "COTS_TEMPORAL_SPLIT_YEAR"
            ]
        ),
    )
    open(path, "w") do io
        TOML.print(io, metadata; sorted=true)
    end
    return metadata
end

function complete_run_metadata(path::AbstractString; best_loss::Float64, evaluations::Int)
    metadata = TOML.parsefile(path)
    metadata["run"]["completed_at"] = string(now())
    metadata["result"] = Dict("best_loss" => best_loss, "logged_evaluations" => evaluations)
    open(path, "w") do io
        TOML.print(io, metadata; sorted=true)
    end
    return metadata
end
