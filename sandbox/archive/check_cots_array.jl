using Pkg; Pkg.activate(".")
using ADRIA

dom = ADRIA.load_domain(ADRIA.LizardDomain, "sandbox/data/Lizard_Island_Cluster_v0.1", "45")
p_cots = ADRIA.param_table(dom)
rs_cots = ADRIA.run_scenario(dom, p_cots[1, :])

cots_pop = rs_cots.cots_log
println("T=1 to 10 adults: ", [sum(cots_pop[t, 3, :]) for t in 1:10])
println("T=1 to 10 juveniles: ", [sum(cots_pop[t, 2, :]) for t in 1:10])
