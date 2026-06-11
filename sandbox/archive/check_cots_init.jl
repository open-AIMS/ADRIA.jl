using Pkg; Pkg.activate(".")
using ADRIA
using Plots

dom = ADRIA.load_domain(ADRIA.LizardDomain, "sandbox/data/Lizard_Island_Cluster_v0.1", "45")

cots_init = dom.cots_init_density
n_sites = length(cots_init)
println("Total sites: ", n_sites)

println("Sites > 0.9: ", sum(cots_init .> 0.9))
println("Sites > 0.7: ", sum(cots_init .> 0.7))
println("Sites > 0.5: ", sum(cots_init .> 0.5))
println("Sites > 0.3: ", sum(cots_init .> 0.3))
println("Sites > 0.1: ", sum(cots_init .> 0.1))
println("Sites > 0.0: ", sum(cots_init .> 0.0))

# Also let's run a quick simulation and save a plot of total COTS population over time
p_cots = ADRIA.param_table(dom)
rs_cots = ADRIA.run_scenario(dom, p_cots[1, :])

cots_pop = rs_cots.cots_log

# cots_log has dimensions (timesteps, stages, locations)
# Let's sum across locations
# The array is actually (timesteps, stages, locations) or (timesteps, locations, scenarios) ?
# In scenario.jl it's Ycots which has shape (tf, 3, n_locs)
total_cots_over_time = [sum(cots_pop[t, 3, :]) for t in 1:size(cots_pop, 1)] # sum adults

println("Initial total adult COTS: ", total_cots_over_time[1])
println("Max total adult COTS: ", maximum(total_cots_over_time))

plt = plot(1:length(total_cots_over_time), total_cots_over_time, 
           xlabel="Timestep (Years)", ylabel="Total Adult COTS Density", 
           title="Adult COTS Density Over Time (Threshold > 0.7)",
           label="Adult COTS", lw=2)

savefig(plt, "sandbox/cots_pop_test.png")
println("Saved plot to sandbox/cots_pop_test.png")
