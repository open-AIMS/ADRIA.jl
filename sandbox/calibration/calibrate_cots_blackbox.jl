using Pkg
const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
Pkg.activate(joinpath(REPO_ROOT, "sandbox"))
cd(REPO_ROOT)

using ADRIA
using CSV, DataFrames, Statistics, Dates
using BlackBoxOptim

include(joinpath(@__DIR__, "cots_cycle_metrics.jl"))

println("=== ADRIA COTS Submodel BlackBoxOptim Calibration Driver ===")

# 1. Load historical domain & observation mapping
println("Loading Lizard Island historical domain...")
dom = ADRIA.load_domain(ADRIA.LizardDomain, "sandbox/data/Lizard_Historical_v0.1", "historical")
site_to_reef = CSV.read("sandbox/data/Lizard_Historical_v0.1/site_to_reef.csv", DataFrame)
emp_df = CSV.read("sandbox/data/reef_cots.csv", DataFrame)

site_to_reef.reef_name_clean = [split(r, " (")[1] for r in site_to_reef.reef_name]
target_reefs = ["Lizard Island Reef", "MacGillivray Reef", "North Direction Reef", "Eyrie Reef"]
sim_to_emp_map = Dict(
    "Lizard Island Reef" => "Lizard Isles",
    "MacGillivray Reef" => "Macgillivray Reef",
    "North Direction Reef" => "North Direction Island",
    "Eyrie Reef" => "Eyrie Reef"
)

# Set base environment defaults
ENV["COTS_EXTERNAL_PULSE"] = get(ENV, "COTS_EXTERNAL_PULSE", "false")
ENV["COTS_SEED_FIRST_N"] = get(ENV, "COTS_SEED_FIRST_N", "10")
ENV["ADRIA_DEBUG_SEED_FIRST_N"] = ENV["COTS_SEED_FIRST_N"]

# Build base scenario parameter template
scen_template = ADRIA.sample(dom, 2)[1:1, :]
p_df = ADRIA.param_table(dom)
for col in names(scen_template)
    scen_template[!, col] .= p_df[1, col]
end

# 2. Configure Search Mode & Space
mode = uppercase(get(ENV, "BBO_MODE", "FOCUSED"))

param_names = String[]
search_space = Tuple{Float64, Float64}[]

if mode == "FOCUSED"
    param_names = ["a_F", "a_S", "IMM", "seed_mult"]
    search_space = [
        (0.1, 2.0),   # a_F
        (0.01, 0.9),  # a_S
        (0.0, 0.05),  # IMM
        (0.5, 4.0)    # seed_mult
    ]
elseif mode == "EXPANDED"
    param_names = [
        "a_F", "a_S", "IMM", "seed_mult",
        "a_ricker", "b_ricker", "m1", "m2", "m3",
        "p_tilde", "C_max", "tau_condition", "allee_threshold",
        "imm_threshold", "eta_imm"
    ]
    search_space = [
        (0.1, 2.0),   # a_F
        (0.01, 0.9),  # a_S
        (0.0, 0.05),  # IMM
        (0.5, 4.0),   # seed_mult
        (2.0, 10.0),  # a_ricker
        (0.01, 0.5),  # b_ricker
        (0.1, 0.9),   # m1
        (0.05, 0.5),  # m2
        (0.05, 0.3),  # m3
        (0.8, 1.0),   # p_tilde
        (0.4, 1.0),   # C_max
        (1.0, 10.0),  # tau_condition
        (0.1, 5.0),   # allee_threshold
        (0.1, 0.8),   # imm_threshold
        (1.0, 5.0)    # eta_imm
    ]
elseif mode == "EXPANDED_PULSE"
    param_names = [
        "a_F", "a_S", "IMM", "seed_mult",
        "a_ricker", "b_ricker", "m1", "m2", "m3",
        "p_tilde", "C_max", "tau_condition", "allee_threshold",
        "imm_threshold", "eta_imm",
        "pulse_start", "pulse_duration", "pulse_relative_magnitude"
    ]
    search_space = [
        (0.1, 2.0),   # a_F
        (0.01, 0.9),  # a_S
        (0.0, 0.05),  # IMM
        (0.5, 4.0),   # seed_mult
        (2.0, 10.0),  # a_ricker
        (0.01, 0.5),  # b_ricker
        (0.1, 0.9),   # m1
        (0.05, 0.5),  # m2
        (0.05, 0.3),  # m3
        (0.8, 1.0),   # p_tilde
        (0.4, 1.0),   # C_max
        (1.0, 10.0),  # tau_condition
        (0.1, 5.0),   # allee_threshold
        (0.1, 0.8),   # imm_threshold
        (1.0, 5.0),   # eta_imm
        (15.0, 30.0), # pulse_start
        (1.0, 4.0),   # pulse_duration
        (0.0, 1.0)    # pulse_relative_magnitude
    ]
else
    error("Unknown BBO_MODE: '$mode'. Supported options: FOCUSED, EXPANDED, EXPANDED_PULSE")
end

println("Search Mode: $mode")
println("Parameters ($(length(param_names))): ", join(param_names, ", "))

# 3. Setup Detailed Candidate Logging
candidates_log_file = joinpath(REPO_ROOT, "sandbox", "data", "bbo_evaluated_candidates.csv")
summary_best_file = joinpath(REPO_ROOT, "sandbox", "data", "bbo_best_summary.csv")

log_columns = vcat(
    ["eval_id", "timestamp", "loss", "cycle_loss", "legacy_loss"],
    param_names,
    [
        "mean_rmse", "mean_abs_percent_bias", "mean_pearson", "mean_spearman",
        "mean_peak_count_penalty", "mean_peak_timing_penalty", "mean_period_penalty",
        "mean_amplitude_penalty", "mean_flatline_penalty", "mean_lag_correlation_penalty",
        "mean_best_lag_years", "mean_best_lag_pearson", "mean_best_lag_spearman"
    ]
)

eval_log_df = DataFrame([col => Any[] for col in log_columns])
CSV.write(candidates_log_file, eval_log_df) # Initialize CSV with header

eval_counter = 0

# Helper metric functions
function pearson_corr(x::Vector{Float64}, y::Vector{Float64})::Float64
    length(x) < 2 && return NaN
    std(x) == 0.0 && return NaN
    std(y) == 0.0 && return NaN
    return cor(x, y)
end

function rank_avg(v::Vector{Float64})::Vector{Float64}
    order = sortperm(v)
    ranks = zeros(Float64, length(v))
    i = 1
    while i <= length(v)
        j = i
        while j < length(v) && v[order[j + 1]] == v[order[i]]
            j += 1
        end
        avg_rank = (i + j) / 2
        for k in i:j
            ranks[order[k]] = avg_rank
        end
        i = j + 1
    end
    return ranks
end

# 4. Objective Function Driver
function evaluate_candidate(candidate_vec::Vector{Float64})::Float64
    global eval_counter += 1
    curr_id = eval_counter
    ts_str = Dates.format(now(), "yyyy-mm-dd HH:MM:SS")

    # Map candidate values to dictionary
    cand_dict = Dict(name => candidate_vec[i] for (i, name) in enumerate(param_names))

    # Configure scenario parameter row
    scen_df = copy(scen_template)

    # Apply factor parameter updates
    for name in param_names
        if hasproperty(scen_df, Symbol(name))
            scen_df[1, Symbol(name)] = cand_dict[name]
        end
    end


    # Handle ENV-controlled parameters
    if haskey(cand_dict, "seed_mult")
        ENV["COTS_INITIAL_MULTIPLIER"] = string(cand_dict["seed_mult"])
    end

    if haskey(cand_dict, "pulse_relative_magnitude")
        rel_mag = cand_dict["pulse_relative_magnitude"]
        if rel_mag > 0.05
            ENV["COTS_EXTERNAL_PULSE"] = "true"
            ENV["COTS_PULSE_START"] = string(round(Int, cand_dict["pulse_start"]))
            ENV["COTS_PULSE_DURATION"] = string(round(Int, cand_dict["pulse_duration"]))
            ENV["COTS_PULSE_RELATIVE_MAGNITUDE"] = string(rel_mag)
        else
            ENV["COTS_EXTERNAL_PULSE"] = "false"
        end
    end

    # Run scenario
    local rs
    try
        rs = ADRIA.run_scenario(dom, scen_df[1, :])
    catch err
        @warn "Scenario evaluation failed for candidate $curr_id: $err"
        return 9999.0
    end

    adult_cots_site = rs.cots_log[:, 3, :]
    n_timesteps = size(adult_cots_site, 1)

    # Calculate simulated adult COTS normalized per reef
    sim_df = DataFrame(reef_name=String[], year=Int[], sim_cots_adult=Float64[], sim_cots_norm=Float64[])
    for reef in unique(site_to_reef.reef_name_clean)
        site_indices = findall(site_to_reef.reef_name_clean .== reef)
        isempty(site_indices) && continue
        reef_sim_cots = [mean(adult_cots_site[t, site_indices]) for t in 1:n_timesteps]
        for t in 1:n_timesteps
            push!(sim_df, (reef, 1984 + t, reef_sim_cots[t], 0.0))
        end
    end

    for reef in unique(sim_df.reef_name)
        reef_rows = sim_df.reef_name .== reef
        max_sim = maximum(sim_df[reef_rows, :sim_cots_adult])
        if max_sim > 0
            sim_df[reef_rows, :sim_cots_norm] .= sim_df[reef_rows, :sim_cots_adult] ./ max_sim
        end
    end

    # Evaluate cycle metrics per target reef
    reef_rmses = Float64[]
    reef_biases = Float64[]
    reef_pearsons = Float64[]
    reef_spearmans = Float64[]

    cycle_losses = Float64[]
    peak_count_penalties = Float64[]
    peak_timing_penalties = Float64[]
    period_penalties = Float64[]
    amplitude_penalties = Float64[]
    flatline_penalties = Float64[]
    lag_corr_penalties = Float64[]

    best_lag_years_vec = Float64[]
    best_lag_pearsons = Float64[]
    best_lag_spearmans = Float64[]

    for reef in target_reefs
        reef_sim = sim_df[sim_df.reef_name .== reef, :]
        isempty(reef_sim) && continue

        emp_name = sim_to_emp_map[reef]
        emp_reef = emp_df[emp_df.reef_name .== emp_name, :]
        isempty(emp_reef) && continue
        max_obs = maximum(emp_reef.cotsptow)
        max_obs <= 0 && continue

        obs_by_year = combine(groupby(transform(emp_reef, :cotsptow => ByRow(x -> x / max_obs) => :obs_norm), :year), :obs_norm => mean => :obs)
        sim_by_year = combine(groupby(reef_sim, :year), :sim_cots_norm => median => :sim)
        joined = innerjoin(sim_by_year, obs_by_year; on=:year)
        isempty(joined) && continue

        sim_vec = Vector{Float64}(joined.sim)
        obs_vec = Vector{Float64}(joined.obs)

        # Baseline metrics
        r = pearson_corr(sim_vec, obs_vec)
        rho = pearson_corr(rank_avg(sim_vec), rank_avg(obs_vec))
        rmse = sqrt(mean((sim_vec .- obs_vec) .^ 2))
        pbias = sum(obs_vec) == 0.0 ? 0.0 : 100.0 * (sum(sim_vec) - sum(obs_vec)) / sum(obs_vec)

        push!(reef_rmses, rmse)
        push!(reef_biases, abs(pbias))
        push!(reef_pearsons, isnan(r) ? 0.0 : r)
        push!(reef_spearmans, isnan(rho) ? 0.0 : rho)

        # Detailed COTS cycle score
        cycle = cots_cycle_score(
            Vector{Int}(sim_by_year.year),
            Vector{Float64}(sim_by_year.sim),
            Vector{Int}(obs_by_year.year),
            Vector{Float64}(obs_by_year.obs)
        )

        push!(cycle_losses, cycle.total_loss)
        push!(peak_count_penalties, cycle.peak_count_penalty)
        push!(peak_timing_penalties, cycle.peak_timing_penalty)
        push!(period_penalties, cycle.period_penalty)
        push!(amplitude_penalties, cycle.amplitude_penalty)
        push!(flatline_penalties, cycle.flatline_penalty)
        push!(lag_corr_penalties, cycle.lag_correlation_penalty)
        push!(best_lag_years_vec, Float64(cycle.best_lag_years))
        push!(best_lag_pearsons, cycle.best_lag_pearson)
        push!(best_lag_spearmans, cycle.best_lag_spearman)
    end

    if isempty(cycle_losses)
        return 9999.0
    end

    mean_cycle_loss = mean(cycle_losses)
    mean_rmse = mean(reef_rmses)
    mean_abs_pbias = mean(reef_biases)
    mean_pearson = mean(reef_pearsons)
    mean_spearman = mean(reef_spearmans)

    legacy_loss = mean_rmse + 0.002 * mean_abs_pbias - 0.05 * mean_pearson
    total_loss = mean_cycle_loss + 0.25 * legacy_loss

    # Construct evaluation log row
    log_row = Dict{String, Any}(
        "eval_id" => curr_id,
        "timestamp" => ts_str,
        "loss" => total_loss,
        "cycle_loss" => mean_cycle_loss,
        "legacy_loss" => legacy_loss,
        "mean_rmse" => mean_rmse,
        "mean_abs_percent_bias" => mean_abs_pbias,
        "mean_pearson" => mean_pearson,
        "mean_spearman" => mean_spearman,
        "mean_peak_count_penalty" => mean(peak_count_penalties),
        "mean_peak_timing_penalty" => mean(peak_timing_penalties),
        "mean_period_penalty" => mean(period_penalties),
        "mean_amplitude_penalty" => mean(amplitude_penalties),
        "mean_flatline_penalty" => mean(flatline_penalties),
        "mean_lag_correlation_penalty" => mean(lag_corr_penalties),
        "mean_best_lag_years" => mean(best_lag_years_vec),
        "mean_best_lag_pearson" => mean(best_lag_pearsons),
        "mean_best_lag_spearman" => mean(best_lag_spearmans)
    )

    for (name, val) in cand_dict
        log_row[name] = val
    end

    row_data = [get(log_row, col, NaN) for col in log_columns]
    push!(eval_log_df, row_data)

    # Append to CSV on disk
    CSV.write(candidates_log_file, eval_log_df[end:end, :], append=true)

    println("Eval #$(curr_id) | Loss: $(round(total_loss; digits=4)) | Cycle Loss: $(round(mean_cycle_loss; digits=4)) | Lag: $(round(mean(best_lag_years_vec); digits=1)) yrs")
    return total_loss
end

# 5. Execute BlackBoxOptim Setup & Run
max_steps = parse(Int, get(ENV, "BBO_MAX_STEPS", "20"))
max_time = parse(Float64, get(ENV, "BBO_MAX_TIME", "0.0"))

println("\nStarting BlackBoxOptim (Max Steps: $max_steps, Max Time: $max_time s)...")

opt_args = Dict{Symbol, Any}(
    :SearchSpace => search_space,
    :NumDimensions => length(param_names),
    :MaxSteps => max_steps,
    :TraceMode => :silent
)
if max_time > 0.0
    opt_args[:MaxTime] = max_time
end

opt_res = bboptimize(evaluate_candidate; opt_args...)

best_cand = best_candidate(opt_res)
best_fit = best_fitness(opt_res)

println("\n================ Optimization Completed ================")
println("Best Loss Score: ", round(best_fit; digits=4))
println("Best Candidate Vector:")
for (i, name) in enumerate(param_names)
    println("  $name = ", round(best_cand[i]; digits=6))
end

# Save Best Candidates Summary
if nrow(eval_log_df) > 0
    best_summary_df = sort(eval_log_df, :loss)
    CSV.write(summary_best_file, best_summary_df)
    println("\nEvaluated candidates logged to: $candidates_log_file")
    println("Top candidates summary saved to: $summary_best_file")
end
