using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
using ADRIA
using Plots
using Statistics
using DataFrames

DOMAIN_PATH = joinpath(@__DIR__, "data", "Moore_2025-08-22_v080")
dom = ADRIA.load_domain(DOMAIN_PATH, "26")

# Find the reef grouping column
reef_col = :Reef
if !(:Reef in propertynames(dom.loc_data))
    reef_col = :reef_siteid
end

println("Using reef column: ", reef_col)
reefs = unique(dom.loc_data[!, reef_col])
println("Found ", length(reefs), " reefs")

n_reefs_to_seed = get(ENV, "ADRIA_DEBUG_SEED_N_REEFS", "all")
if n_reefs_to_seed == "all"
    target_reefs = reefs
else
    n = parse(Int, n_reefs_to_seed)
    target_reefs = reefs[1:min(n, length(reefs))]
end

println("Target reefs for seeding: ", target_reefs)

seed_locs = Int[]
for r in target_reefs
    locs_in_reef = findall(dom.loc_data[!, reef_col] .== r)
    # Seed >50% (60%) of the locations in this reef
    n_seed = ceil(Int, length(locs_in_reef) * 0.6)
    append!(seed_locs, locs_in_reef[1:n_seed])
end
println("Seeding locations: ", seed_locs)

ENV["ADRIA_COTS_ENABLED"] = "true"
ENV["COTS_EXTERNAL_PULSE"] = "true"
ENV["ADRIA_DEBUG_INIT_DENSITY"] = "2.0"
ENV["ADRIA_DEBUG_SEED_LOCATIONS"] = join(seed_locs, ",")
dom.init_coral_cover.data .*= 2.0 

tf = length(ADRIA.timesteps(dom))
scens = ADRIA.sample(dom, 2) # Sobol needs at least 2
res = ADRIA.run_scenario(dom, scens[1, :])

ENV["ADRIA_COTS_ENABLED"] = "false"
res_no_cots = ADRIA.run_scenario(dom, scens[1, :])
ENV["ADRIA_COTS_ENABLED"] = "true" # Reset for consistency

# Aggregate by reef
plot_list = []

# Extract raw data robustly to handle YAXArrays
raw_data = hasproperty(res.raw, :data) ? res.raw.data : res.raw
# Sum over functional groups (dim 2) and size classes (dim 3)
total_coral_per_loc = dropdims(sum(raw_data, dims=(2, 3)), dims=(2, 3)) # Shape: [time, loc]

raw_data_no_cots = hasproperty(res_no_cots.raw, :data) ? res_no_cots.raw.data : res_no_cots.raw
total_coral_per_loc_no_cots = dropdims(sum(raw_data_no_cots, dims=(2, 3)), dims=(2, 3)) # Shape: [time, loc]

cots_data = hasproperty(res.cots_log, :data) ? res.cots_log.data : res.cots_log
# adult COTS is age class 3 (dim 2)
adult_cots_per_loc = cots_data[:, 3, :] # Shape: [time, loc]

for (i, r) in enumerate(reefs)
    # Get locs
    locs = findall(dom.loc_data[!, reef_col] .== r)
    if length(locs) == 0
        continue
    end
    
    # Mean across locations in this reef
    reef_coral = vec(mean(total_coral_per_loc[:, locs], dims=2))
    reef_cots = vec(mean(adult_cots_per_loc[:, locs], dims=2))
    
    # Base plot
    p = plot(title="Reef: $r ($(length(locs)) locs)", ylabel="Percentage / Scaled Density", legend=:topleft, ylims=(0, maximum(reef_coral)*120))
    
    # Plot individual location trajectories as semi-transparent lines
    for loc in locs
        loc_coral = total_coral_per_loc[:, loc] .* 100
        loc_cots = adult_cots_per_loc[:, loc] .* 20
        plot!(p, 1:tf, loc_coral, color=:blue, alpha=0.15, label="")
        plot!(p, 1:tf, loc_cots, color=:red, alpha=0.15, label="")
    end
    
    # Overlay the reef mean lines (bold)
    plot!(p, 1:tf, reef_coral .* 100, label="Coral % Mean", color=:blue, linewidth=3)
    plot!(p, 1:tf, reef_cots .* 20, label="Adult COTS (x20)", color=:red, linewidth=3)
              
    push!(plot_list, p)
end

# Extract DHW Data
dhw_scen_id = Int(scens[1, :dhw_scenario])
sz = size(dom.dhw_scens)
dhw_matrix = if sz[2] == size(dom.loc_data, 1)
    dom.dhw_scens[:, :, dhw_scen_id]
else
    dom.dhw_scens[:, dhw_scen_id, :]
end
dhw_data = hasproperty(dhw_matrix, :data) ? dhw_matrix.data : dhw_matrix

dhw_plot_list = []
for (i, r) in enumerate(reefs)
    locs = findall(dom.loc_data[!, reef_col] .== r)
    if length(locs) == 0
        continue
    end
    
    p2 = plot(title="Reef: $r DHW", ylabel="Degree Heating Weeks", legend=false, ylims=(0, maximum(dhw_data) * 1.2))
    
    for loc in locs
        plot!(p2, 1:tf, dhw_data[:, loc], color=:orange, alpha=0.3)
    end
    
    reef_dhw = vec(mean(dhw_data[:, locs], dims=2))
    plot!(p2, 1:tf, reef_dhw, color=:red, linewidth=3)
    
    push!(dhw_plot_list, p2)
end

n_plots = length(plot_list)
cols = 2
rows = ceil(Int, n_plots / cols)
final_plot = plot(plot_list..., layout=(rows, cols), size=(1200, rows * 350), margin=5Plots.mm)
savefig(final_plot, joinpath(@__DIR__, "faceted_reef_dynamics_with_disturbances.png"))
println("Saved faceted plot to faceted_reef_dynamics_with_disturbances.png")

final_dhw_plot = plot(dhw_plot_list..., layout=(rows, cols), size=(1200, rows * 350), margin=5Plots.mm)
savefig(final_dhw_plot, joinpath(@__DIR__, "faceted_reef_dhw.png"))
println("Saved faceted DHW plot to faceted_reef_dhw.png")

comp_plot_list = []
for (i, r) in enumerate(reefs)
    locs = findall(dom.loc_data[!, reef_col] .== r)
    if length(locs) == 0
        continue
    end
    
    reef_coral = vec(mean(total_coral_per_loc[:, locs], dims=2))
    reef_coral_no_cots = vec(mean(total_coral_per_loc_no_cots[:, locs], dims=2))
    
    max_val = max(maximum(reef_coral), maximum(reef_coral_no_cots)) * 100
    p3 = plot(title="Reef: $r (w/ vs w/o COTS)", ylabel="Coral %", legend=:topleft, ylims=(0, max_val*1.2))
    
    for loc in locs
        loc_coral = total_coral_per_loc[:, loc] .* 100
        loc_coral_no_cots = total_coral_per_loc_no_cots[:, loc] .* 100
        plot!(p3, 1:tf, loc_coral_no_cots, color=:green, alpha=0.15, label="")
        plot!(p3, 1:tf, loc_coral, color=:blue, alpha=0.15, label="")
    end
    
    plot!(p3, 1:tf, reef_coral_no_cots .* 100, label="No COTS Mean", color=:green, linewidth=3)
    plot!(p3, 1:tf, reef_coral .* 100, label="With COTS Mean", color=:blue, linewidth=3)
    
    push!(comp_plot_list, p3)
end

final_comp_plot = plot(comp_plot_list..., layout=(rows, cols), size=(1200, rows * 350), margin=5Plots.mm)
savefig(final_comp_plot, joinpath(@__DIR__, "faceted_reef_coral_comparison.png"))
println("Saved faceted comparison plot to faceted_reef_coral_comparison.png")
