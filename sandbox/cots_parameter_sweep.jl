using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
using ADRIA
using Plots
using Statistics
using DataFrames
using CSV
using Base.Iterators

DOMAIN_PATH = joinpath(@__DIR__, "data", "Moore_2025-08-22_v080")
dom = ADRIA.load_domain(DOMAIN_PATH, "26")

# Find the reef grouping column
reef_col = :Reef
if !(:Reef in propertynames(dom.loc_data))
    reef_col = :reef_siteid
end
println("Using reef column: ", reef_col)

# Setup seeding specifically on Elford_16073
target_reef = "Elford_16073"
locs_in_reef = findall(dom.loc_data[!, reef_col] .== target_reef)
n_seed = ceil(Int, length(locs_in_reef) * 0.6)
seed_locs = locs_in_reef[1:n_seed]
println("Seeding locations: ", seed_locs)

ENV["ADRIA_COTS_ENABLED"] = "true"
ENV["ADRIA_DEBUG_SEED_LOCATIONS"] = join(seed_locs, ",")
ENV["ADRIA_DEBUG_BLEACHING_SCALAR"] = "0.0" 
ENV["ADRIA_DEBUG_CYCLONE_SCALAR"] = "0.0"   

# Boost initial coral cover to give COTS something to eat initially
dom.init_coral_cover.data .*= 2.0 

tf = length(ADRIA.timesteps(dom))
scens = ADRIA.sample(dom, 2) # Sobol needs at least 2

ENV["COTS_allee_threshold"] = "1.0"
ENV["ADRIA_DEBUG_INIT_DENSITY"] = "2.0"

# Parameter grids
a_rickers = [2.0, 4.0, 6.0]
b_rickers = [0.1, 0.5, 0.9]
p_tildes = [0.90, 0.95]
tau_conditions = [1.0, 3.0, 5.0]

results_df = DataFrame(
    a_ricker = Float64[],
    b_ricker = Float64[],
    p_tilde = Float64[],
    tau_condition = Float64[],
    num_waves = Int[],
    avg_cycle_length = Float64[],
    avg_min_cots = Float64[]
)

out_dir = joinpath(@__DIR__, "sweep_results")
mkpath(out_dir)

function count_waves(cots_timeseries, min_height=0.5)
    # Find peaks in COTS density
    peaks = Int[]
    for i in 2:(length(cots_timeseries)-1)
        if cots_timeseries[i] > cots_timeseries[i-1] && cots_timeseries[i] > cots_timeseries[i+1]
            if cots_timeseries[i] > min_height
                push!(peaks, i)
            end
        end
    end
    
    # Filter out peaks that are too close together (within 5 years) to avoid noise
    filtered_peaks = Int[]
    last_peak = -10
    for p in peaks
        if (p - last_peak) >= 5
            push!(filtered_peaks, p)
            last_peak = p
        end
    end
    
    num_waves = length(filtered_peaks)
    avg_cycle = 0.0
    avg_min = 0.0
    
    if num_waves > 1
        avg_cycle = mean(diff(filtered_peaks))
        mins = Float64[]
        for i in 1:(num_waves-1)
            push!(mins, minimum(cots_timeseries[filtered_peaks[i]:filtered_peaks[i+1]]))
        end
        avg_min = mean(mins)
    elseif num_waves == 1
        avg_min = minimum(cots_timeseries[filtered_peaks[1]:end])
    end
    
    return num_waves, avg_cycle, avg_min, filtered_peaks
end

Plots.default(size=(800, 400), margin=5Plots.mm)

run_idx = 1
for (a, b, p_t, tau) in product(a_rickers, b_rickers, p_tildes, tau_conditions)
    println("Run $run_idx: a=$a, b=$b, p_tilde=$p_t, tau=$tau")
    
    ENV["COTS_a_ricker"] = string(a)
    ENV["COTS_b_ricker"] = string(b)
    ENV["COTS_p_tilde"] = string(p_t)
    ENV["COTS_tau_condition"] = string(tau)
    
    # Run simulation
    res = ADRIA.run_scenario(dom, scens[1, :])
    
    raw_data = hasproperty(res.raw, :data) ? res.raw.data : res.raw
    total_coral_per_loc = dropdims(sum(raw_data, dims=(2, 3)), dims=(2, 3)) # [time, loc]
    
    cots_data = hasproperty(res.cots_log, :data) ? res.cots_log.data : res.cots_log
    adult_cots_per_loc = cots_data[:, 3, :] # [time, loc]
    
    reef_coral = vec(mean(total_coral_per_loc[:, locs_in_reef], dims=2))
    reef_cots = vec(mean(adult_cots_per_loc[:, locs_in_reef], dims=2))
    
    num_waves, avg_cycle, avg_min, peaks = count_waves(reef_cots, 0.5)
    
    push!(results_df, (a, b, p_t, tau, num_waves, avg_cycle, avg_min))
    
    if num_waves > 1
        # Save plot
        p = plot(1:tf, reef_coral .* 100, label="Coral %", color=:blue, linewidth=2, 
                 title="a=$a b=$b p_tilde=$p_t tau=$tau\nWaves: $num_waves | Min COTS: $(round(avg_min, digits=4))", 
                 ylabel="Percentage / Scaled Density", legend=:topleft, ylims=(0, maximum(reef_coral)*120))
                 
        p = plot!(1:tf, reef_cots .* 20, label="Adult COTS (x20)", color=:red, linewidth=2)
        
        # Mark peaks
        scatter!(peaks, reef_cots[peaks] .* 20, color=:red, markersize=5, label="")
        
        savefig(p, joinpath(out_dir, "plot_run_$run_idx.png"))
    end
    
    global run_idx += 1
end

CSV.write(joinpath(out_dir, "sweep_summary.csv"), results_df)
println("Sweep complete. Saved to $out_dir")
