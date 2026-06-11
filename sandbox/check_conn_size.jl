using Pkg; Pkg.activate(".")
using ADRIA

println("Loading domain...")
dom = ADRIA.load_domain(ADRIA.LizardDomain, "sandbox/data/Lizard_Island_Cluster_v0.1", "45")

println("conn size: ", size(dom.conn))
println("loc_data size: ", size(dom.loc_data, 1))
