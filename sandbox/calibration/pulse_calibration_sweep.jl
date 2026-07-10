using Pkg
const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
Pkg.activate(joinpath(REPO_ROOT, "sandbox"))
cd(REPO_ROOT)
using ADRIA
using CSV, DataFrames, Statistics

function parse_int_grid(env_name::String, default::String)::Vector{Int}
    return parse.(Int, split(get(ENV, env_name, default), ","))
end

function parse_float_grid(env_name::String, default::String)::Vector{Float64}
    return parse.(Float64, split(get(ENV, env_name, default), ","))
end

function pearson(x::Vector{Float64}, y::Vector{Float64})::Float64
    length(x) < 2 && return NaN
    std(x) == 0.0 && return NaN
    std(y) == 0.0 && return NaN
    return cor(x, y)
end

function rank_average(v::Vector{Float64})::Vector{Float64}
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

function validation_metrics(sim::Vector{Float64}, obs::Vector{Float64})::NamedTuple
    n = length(sim)
    n == 0 && return (n=0, pearson=NaN, spearman=NaN, rmse=NaN, percent_bias=NaN)
    r = pearson(sim, obs)
    rho = pearson(rank_average(sim), rank_average(obs))
    rmse = sqrt(mean((sim .- obs) .^ 2))
    percent_bias = sum(obs) == 0.0 ? NaN : 100.0 * (sum(sim) - sum(obs)) / sum(obs)
    return (n=n, pearson=r, spearman=rho, rmse=rmse, percent_bias=percent_bias)
end

function reef_validation_rows(
    sim_df::DataFrame,
    emp_df::DataFrame,
    sim_to_emp_map::Dict{String,String},
    reef_names::Vector{String},
    variant::NamedTuple
)::DataFrame
    rows = DataFrame(
        pulse_start=Int[],
        pulse_duration=Int[],
        pulse_repeat_interval=Int[],
        pulse_relative_magnitude=Float64[],
        reef_name=String[],
        n_obs=Int[],
        pearson=Float64[],
        spearman=Float64[],
        rmse=Float64[],
        percent_bias=Float64[],
        post_2010_peak_year=Int[]
    )

    for reef in reef_names
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

        m = validation_metrics(Vector{Float64}(joined.sim), Vector{Float64}(joined.obs))
        median_by_year = combine(groupby(reef_sim, :year), :sim_cots_norm => median => :median)
        post_2010 = median_by_year[median_by_year.year .>= 2010, :]
        peak_year = isempty(post_2010) ? missing : post_2010.year[argmax(post_2010.median)]

        push!(
            rows,
            (
                variant.start,
                variant.duration,
                variant.repeat_interval,
                variant.relative_magnitude,
                reef,
                m.n,
                m.pearson,
                m.spearman,
                m.rmse,
                m.percent_bias,
                peak_year
            )
        )
    end

    return rows
end

println("Loading Lizard Island historical domain...")
dom = ADRIA.load_domain(ADRIA.LizardDomain, "sandbox/data/Lizard_Historical_v0.1", "historical")
site_to_reef = CSV.read("sandbox/data/Lizard_Historical_v0.1/site_to_reef.csv", DataFrame)
emp_df = CSV.read("sandbox/data/reef_cots.csv", DataFrame)

site_to_reef.reef_name_clean = [split(r, " (")[1] for r in site_to_reef.reef_name]
reef_names = ["Lizard Island Reef", "MacGillivray Reef", "North Direction Reef", "Eyrie Reef"]
sim_to_emp_map = Dict(
    "Lizard Island Reef" => "Lizard Isles",
    "MacGillivray Reef" => "Macgillivray Reef",
    "North Direction Reef" => "North Direction Island",
    "Eyrie Reef" => "Eyrie Reef"
)

best_a_F = 1.34378
best_a_S = 0.0957831
best_IMM = 0.0947791
best_seed_mult = 2.23695

starts = parse_int_grid("PULSE_SWEEP_STARTS", "16,18,20,22,24,26,28,30")
durations = parse_int_grid("PULSE_SWEEP_DURATIONS", "1,2,3")
repeat_intervals = parse_int_grid("PULSE_SWEEP_REPEAT_INTERVALS", "0")
relative_magnitudes = parse_float_grid("PULSE_SWEEP_RELATIVE_MAGNITUDES", "0.0,0.25,0.5,0.75,1.0")

ENV["COTS_INITIAL_MULTIPLIER"] = string(best_seed_mult)
ENV["COTS_EXTERNAL_PULSE"] = "false"
ENV["COTS_SEED_FIRST_N"] = get(ENV, "COTS_SEED_FIRST_N", "10")
ENV["ADRIA_DEBUG_SEED_FIRST_N"] = ENV["COTS_SEED_FIRST_N"]

scen = ADRIA.sample(dom, 2)[1:1, :]
p_df = ADRIA.param_table(dom)
for col in names(scen)
    scen[!, col] .= p_df[1, col]
end
scen.a_F .= best_a_F
scen.a_S .= best_a_S
scen.IMM .= best_IMM

all_rows = DataFrame()

for start in starts, duration in durations, repeat_interval in repeat_intervals, relative_magnitude in relative_magnitudes
    variant = (
        start=start,
        duration=duration,
        repeat_interval=repeat_interval,
        relative_magnitude=relative_magnitude
    )

    println("Pulse sweep: start=$start duration=$duration repeat=$repeat_interval magnitude=$relative_magnitude")
    ENV["COTS_EXTERNAL_PULSE"] = relative_magnitude > 0.0 ? "true" : "false"
    ENV["COTS_PULSE_START"] = string(start)
    ENV["COTS_PULSE_DURATION"] = string(duration)
    ENV["COTS_PULSE_REPEAT_INTERVAL"] = string(repeat_interval)
    ENV["COTS_PULSE_RELATIVE_MAGNITUDE"] = string(relative_magnitude)

    rs = ADRIA.run_scenario(dom, scen[1, :])
    adult_cots_site = rs.cots_log[:, 3, :]
    n_timesteps = size(adult_cots_site, 1)

    sim_df = DataFrame(sim_id=Int[], year=Int[], reef_name=String[], sim_cots_adult=Float64[], sim_cots_norm=Float64[])
    for reef in unique(site_to_reef.reef_name_clean)
        site_indices = findall(site_to_reef.reef_name_clean .== reef)
        isempty(site_indices) && continue

        reef_sim_cots = [mean(adult_cots_site[t, site_indices]) for t in 1:n_timesteps]
        for t in 1:n_timesteps
            push!(sim_df, (1, 1984 + t, reef, reef_sim_cots[t], 0.0))
        end
    end

    for reef in unique(sim_df.reef_name)
        reef_rows = sim_df.reef_name .== reef
        max_sim = maximum(sim_df[reef_rows, :sim_cots_adult])
        if max_sim > 0
            sim_df[reef_rows, :sim_cots_norm] .= sim_df[reef_rows, :sim_cots_adult] ./ max_sim
        end
    end

    rows = reef_validation_rows(sim_df, emp_df, sim_to_emp_map, reef_names, variant)
    global all_rows = vcat(all_rows, rows)
end

summary = combine(
    groupby(all_rows, [:pulse_start, :pulse_duration, :pulse_repeat_interval, :pulse_relative_magnitude]),
    :rmse => mean => :mean_rmse,
    :pearson => mean => :mean_pearson,
    :spearman => mean => :mean_spearman,
    :percent_bias => (x -> mean(abs.(x))) => :mean_abs_percent_bias
)
summary.loss = summary.mean_rmse .+ 0.002 .* summary.mean_abs_percent_bias .- 0.05 .* coalesce.(summary.mean_pearson, 0.0)
sort!(summary, :loss)

CSV.write("sandbox/data/pulse_calibration_sweep_by_reef.csv", all_rows)
CSV.write("sandbox/data/pulse_calibration_sweep_summary.csv", summary)
println("Saved pulse sweep results to sandbox/data/pulse_calibration_sweep_summary.csv")
println(first(summary, min(10, nrow(summary))))

